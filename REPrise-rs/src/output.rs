//! Production-grade output formatting and export system
//!
//! Provides multiple output formats (BED, TSV, JSON, GFF3) with compression support,
//! metadata headers, and checksums for data integrity verification.

use crate::config::{CompressionType, OutputFormat, OutputSettings};
use crate::error::{RepriseError, Result};
use crate::genome::Genome;
use crate::pipeline::DetectedRepeat;
use bzip2::write::BzEncoder;
use chrono::{DateTime, Utc};
use flate2::{write::GzEncoder, Compression};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;
use uuid::Uuid;

/// Metadata information for output files
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutputMetadata {
    /// Unique run identifier
    pub run_id: Uuid,
    /// Timestamp when the run started
    pub timestamp: DateTime<Utc>,
    /// REPrise version
    pub version: String,
    /// Command line used to generate this output
    pub command_line: String,
    /// Input file path
    pub input_file: PathBuf,
    /// Processing parameters
    pub parameters: ProcessingParameters,
    /// Statistics from the run
    pub statistics: OutputStatistics,
    /// Genome information
    pub genome_info: GenomeInfo,
}

/// Processing parameters used for the analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProcessingParameters {
    pub k_mer_length: usize,
    pub min_frequency: u32,
    pub max_frequency: Option<u32>,
    pub region_extension: u64,
    pub min_alignment_score: i32,
    pub min_identity: f64,
    pub num_workers: usize,
}

/// Statistics from the processing run
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutputStatistics {
    pub total_repeats_found: u64,
    pub candidates_processed: u64,
    pub processing_time_seconds: f64,
    pub genome_coverage: f64,
    pub average_repeat_length: f64,
    pub average_identity: f64,
    pub average_score: f64,
}

/// Genome information for output metadata
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GenomeInfo {
    pub total_length: u64,
    pub num_contigs: usize,
    pub n50: u64,
    pub gc_content: f64,
}

/// A single repeat entry for export
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RepeatEntry {
    /// First region chromosome/contig
    pub chr1: String,
    /// First region start position (0-based)
    pub start1: u64,
    /// First region end position (exclusive)
    pub end1: u64,
    /// Second region chromosome/contig
    pub chr2: String,
    /// Second region start position (0-based)
    pub start2: u64,
    /// Second region end position (exclusive)
    pub end2: u64,
    /// Repeat name/identifier
    pub name: String,
    /// Alignment score
    pub score: i32,
    /// Percent identity
    pub identity: f64,
    /// Alignment length
    pub length: u64,
    /// Strand information
    pub strand1: char,
    pub strand2: char,
    /// Repeat type/class (optional)
    pub repeat_type: Option<String>,
    /// Additional attributes
    pub attributes: HashMap<String, String>,
}

/// Checksum information for output files
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FileChecksum {
    pub file_path: PathBuf,
    pub size_bytes: u64,
    pub sha256: String,
    pub created: DateTime<Utc>,
}

/// Writer that can compress output based on configuration
pub enum CompressedWriter {
    Plain(BufWriter<File>),
    Gzip(BufWriter<GzEncoder<File>>),
    Bzip2(BufWriter<BzEncoder<File>>),
}

impl CompressedWriter {
    /// Create a new compressed writer
    pub fn new(path: &Path, compression: &CompressionType, compression_level: u32) -> Result<Self> {
        let file = OpenOptions::new()
            .create(true)
            .truncate(true)
            .write(true)
            .open(path)
            .map_err(|e| RepriseError::io_error(format!("Failed to create output file: {}", e)))?;

        let writer = match compression {
            CompressionType::None => Self::Plain(BufWriter::new(file)),
            CompressionType::Gzip => {
                let compression = Compression::new(compression_level);
                let encoder = GzEncoder::new(file, compression);
                Self::Gzip(BufWriter::new(encoder))
            }
            CompressionType::Bzip2 => {
                let encoder = BzEncoder::new(file, bzip2::Compression::new(compression_level));
                Self::Bzip2(BufWriter::new(encoder))
            }
        };

        Ok(writer)
    }

    /// Write data to the compressed writer
    pub fn write_all(&mut self, data: &[u8]) -> Result<()> {
        match self {
            Self::Plain(writer) => writer.write_all(data),
            Self::Gzip(writer) => writer.write_all(data),
            Self::Bzip2(writer) => writer.write_all(data),
        }
        .map_err(|e| RepriseError::io_error(format!("Write error: {}", e)))
    }

    /// Flush the writer
    pub fn flush(&mut self) -> Result<()> {
        match self {
            Self::Plain(writer) => writer.flush(),
            Self::Gzip(writer) => writer.flush(),
            Self::Bzip2(writer) => writer.flush(),
        }
        .map_err(|e| RepriseError::io_error(format!("Flush error: {}", e)))
    }
}

impl Write for CompressedWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Self::Plain(writer) => writer.write(buf),
            Self::Gzip(writer) => writer.write(buf),
            Self::Bzip2(writer) => writer.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Self::Plain(writer) => writer.flush(),
            Self::Gzip(writer) => writer.flush(),
            Self::Bzip2(writer) => writer.flush(),
        }
    }
}

/// Output format trait for different export formats
pub trait OutputFormatter {
    /// Write file header with metadata
    fn write_header(&mut self, metadata: &OutputMetadata) -> Result<()>;
    
    /// Write a single repeat entry
    fn write_entry(&mut self, entry: &RepeatEntry) -> Result<()>;
    
    /// Write file footer (if needed)
    fn write_footer(&mut self, total_entries: u64) -> Result<()>;
    
    /// Get file extension for this format
    fn file_extension(&self) -> &str;
    
    /// Check if this format supports metadata headers
    fn supports_metadata(&self) -> bool {
        true
    }
}

/// TSV (Tab-Separated Values) formatter
pub struct TsvFormatter {
    writer: CompressedWriter,
    entry_count: AtomicU64,
}

impl TsvFormatter {
    pub fn new(path: &Path, compression: &CompressionType, compression_level: u32) -> Result<Self> {
        let writer = CompressedWriter::new(path, compression, compression_level)?;
        Ok(Self {
            writer,
            entry_count: AtomicU64::new(0),
        })
    }
}

impl OutputFormatter for TsvFormatter {
    fn write_header(&mut self, metadata: &OutputMetadata) -> Result<()> {
        if self.supports_metadata() {
            writeln!(self.writer, "# REPrise repeat detection results")?;
            writeln!(self.writer, "# Generated: {}", metadata.timestamp)?;
            writeln!(self.writer, "# Version: {}", metadata.version)?;
            writeln!(self.writer, "# Input: {}", metadata.input_file.display())?;
            writeln!(self.writer, "# Parameters: k={}, min_freq={}, min_score={}, min_identity={:.2}%",
                metadata.parameters.k_mer_length,
                metadata.parameters.min_frequency,
                metadata.parameters.min_alignment_score,
                metadata.parameters.min_identity * 100.0)?;
            writeln!(self.writer, "# Statistics: {} repeats, {:.1}s processing time",
                metadata.statistics.total_repeats_found,
                metadata.statistics.processing_time_seconds)?;
        }
        
        writeln!(self.writer, "chr1\tstart1\tend1\tchr2\tstart2\tend2\tname\tscore\tidentity\tlength\tstrand1\tstrand2\ttype")?;
        Ok(())
    }

    fn write_entry(&mut self, entry: &RepeatEntry) -> Result<()> {
        writeln!(self.writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}",
            entry.chr1,
            entry.start1,
            entry.end1,
            entry.chr2,
            entry.start2,
            entry.end2,
            entry.name,
            entry.score,
            entry.identity,
            entry.length,
            entry.strand1,
            entry.strand2,
            entry.repeat_type.as_deref().unwrap_or("unknown"))?;
        
        self.entry_count.fetch_add(1, Ordering::Relaxed);
        Ok(())
    }

    fn write_footer(&mut self, total_entries: u64) -> Result<()> {
        writeln!(self.writer, "# Total entries: {}", total_entries)?;
        self.writer.flush()
    }

    fn file_extension(&self) -> &str {
        "tsv"
    }
}

/// BED (Browser Extensible Data) formatter
pub struct BedFormatter {
    writer: CompressedWriter,
    entry_count: AtomicU64,
}

impl BedFormatter {
    pub fn new(path: &Path, compression: &CompressionType, compression_level: u32) -> Result<Self> {
        let writer = CompressedWriter::new(path, compression, compression_level)?;
        Ok(Self {
            writer,
            entry_count: AtomicU64::new(0),
        })
    }
}

impl OutputFormatter for BedFormatter {
    fn write_header(&mut self, metadata: &OutputMetadata) -> Result<()> {
        if self.supports_metadata() {
            writeln!(self.writer, "track name=\"REPrise Repeats\" description=\"Interspersed repeats detected by REPrise\" visibility=2 itemRgb=\"On\"")?;
            writeln!(self.writer, "# Generated: {}", metadata.timestamp)?;
            writeln!(self.writer, "# REPrise version: {}", metadata.version)?;
        }
        Ok(())
    }

    fn write_entry(&mut self, entry: &RepeatEntry) -> Result<()> {
        // BED format: chr, start, end, name, score, strand, ...
        // Write both regions as separate BED entries
        let color = if entry.identity > 0.8 { "255,0,0" } else if entry.identity > 0.6 { "255,165,0" } else { "255,255,0" };
        
        // First region
        writeln!(self.writer, "{}\t{}\t{}\t{}_1\t{}\t{}\t{}\t{}\t{}\t1\t{}\t0",
            entry.chr1,
            entry.start1,
            entry.end1,
            entry.name,
            entry.score.min(1000), // BED score range 0-1000
            entry.strand1,
            entry.start1, // thickStart
            entry.end1,   // thickEnd
            color,
            entry.end1 - entry.start1)?;

        // Second region
        writeln!(self.writer, "{}\t{}\t{}\t{}_2\t{}\t{}\t{}\t{}\t{}\t1\t{}\t0",
            entry.chr2,
            entry.start2,
            entry.end2,
            entry.name,
            entry.score.min(1000),
            entry.strand2,
            entry.start2,
            entry.end2,
            color,
            entry.end2 - entry.start2)?;

        self.entry_count.fetch_add(2, Ordering::Relaxed); // Two entries per repeat
        Ok(())
    }

    fn write_footer(&mut self, _total_entries: u64) -> Result<()> {
        self.writer.flush()
    }

    fn file_extension(&self) -> &str {
        "bed"
    }
}

/// JSON formatter for structured data export
pub struct JsonFormatter {
    writer: CompressedWriter,
    entries: Vec<RepeatEntry>,
    metadata: Option<OutputMetadata>,
}

impl JsonFormatter {
    pub fn new(path: &Path, compression: &CompressionType, compression_level: u32) -> Result<Self> {
        let writer = CompressedWriter::new(path, compression, compression_level)?;
        Ok(Self {
            writer,
            entries: Vec::new(),
            metadata: None,
        })
    }
}

impl OutputFormatter for JsonFormatter {
    fn write_header(&mut self, metadata: &OutputMetadata) -> Result<()> {
        self.metadata = Some(metadata.clone());
        Ok(())
    }

    fn write_entry(&mut self, entry: &RepeatEntry) -> Result<()> {
        self.entries.push(entry.clone());
        Ok(())
    }

    fn write_footer(&mut self, _total_entries: u64) -> Result<()> {
        let output = serde_json::json!({
            "metadata": self.metadata,
            "repeats": self.entries
        });
        
        let json_str = serde_json::to_string_pretty(&output)
            .map_err(|e| RepriseError::config(format!("JSON serialization error: {}", e)))?;
        
        self.writer.write_all(json_str.as_bytes())?;
        self.writer.flush()
    }

    fn file_extension(&self) -> &str {
        "json"
    }
}

/// GFF3 (General Feature Format) formatter
pub struct Gff3Formatter {
    writer: CompressedWriter,
    entry_count: AtomicU64,
}

impl Gff3Formatter {
    pub fn new(path: &Path, compression: &CompressionType, compression_level: u32) -> Result<Self> {
        let writer = CompressedWriter::new(path, compression, compression_level)?;
        Ok(Self {
            writer,
            entry_count: AtomicU64::new(0),
        })
    }
}

impl OutputFormatter for Gff3Formatter {
    fn write_header(&mut self, metadata: &OutputMetadata) -> Result<()> {
        writeln!(self.writer, "##gff-version 3")?;
        if self.supports_metadata() {
            writeln!(self.writer, "##date {}", metadata.timestamp.format("%Y-%m-%d"))?;
            writeln!(self.writer, "##source REPrise {}", metadata.version)?;
        }
        Ok(())
    }

    fn write_entry(&mut self, entry: &RepeatEntry) -> Result<()> {
        // GFF3 format: seqid, source, type, start, end, score, strand, phase, attributes
        let attributes = format!("ID={};Name={};Target={} {} {};Identity={:.4};Length={}",
            entry.name,
            entry.name,
            entry.chr2,
            entry.start2 + 1, // GFF3 is 1-based
            entry.end2,
            entry.identity,
            entry.length);

        writeln!(self.writer, "{}\tREPrise\trepeat_region\t{}\t{}\t{}\t{}\t.\t{}",
            entry.chr1,
            entry.start1 + 1, // GFF3 is 1-based
            entry.end1,
            entry.score,
            entry.strand1,
            attributes)?;

        self.entry_count.fetch_add(1, Ordering::Relaxed);
        Ok(())
    }

    fn write_footer(&mut self, _total_entries: u64) -> Result<()> {
        writeln!(self.writer, "##FASTA")?; // GFF3 can include sequences
        self.writer.flush()
    }

    fn file_extension(&self) -> &str {
        "gff3"
    }
}

/// Main output manager for handling multiple formats
pub struct OutputManager {
    settings: OutputSettings,
    metadata: OutputMetadata,
    formatters: Vec<Box<dyn OutputFormatter + Send>>,
    output_paths: Vec<PathBuf>,
    checksums: Vec<FileChecksum>,
}

impl OutputManager {
    /// Create a new output manager
    pub fn new(
        settings: OutputSettings,
        metadata: OutputMetadata,
    ) -> Result<Self> {
        let mut output_manager = Self {
            settings,
            metadata,
            formatters: Vec::new(),
            output_paths: Vec::new(),
            checksums: Vec::new(),
        };

        output_manager.initialize_formatters()?;
        Ok(output_manager)
    }

    /// Initialize all configured formatters
    fn initialize_formatters(&mut self) -> Result<()> {
        for format in &self.settings.formats {
            let (path, formatter) = self.create_formatter(format)?;
            self.output_paths.push(path);
            self.formatters.push(formatter);
        }
        Ok(())
    }

    /// Create a formatter for the specified format
    fn create_formatter(
        &self,
        format: &OutputFormat,
    ) -> Result<(PathBuf, Box<dyn OutputFormatter + Send>)> {
        let base_name = &self.settings.file_prefix;
        let extension = format.to_string();
        
        let mut filename = format!("{}.{}", base_name, extension);
        
        // Add compression extension if enabled
        if self.settings.enable_compression {
            filename = match self.settings.compression_type {
                CompressionType::Gzip => format!("{}.gz", filename),
                CompressionType::Bzip2 => format!("{}.bz2", filename),
                CompressionType::None => filename,
            };
        }

        let path = if let Some(dir) = &self.settings.output_dir {
            dir.join(filename)
        } else {
            PathBuf::from(filename)
        };

        let formatter: Box<dyn OutputFormatter + Send> = match format {
            OutputFormat::Tsv => Box::new(TsvFormatter::new(
                &path,
                &self.settings.compression_type,
                self.settings.compression_level,
            )?),
            OutputFormat::Bed => Box::new(BedFormatter::new(
                &path,
                &self.settings.compression_type,
                self.settings.compression_level,
            )?),
            OutputFormat::Json => Box::new(JsonFormatter::new(
                &path,
                &self.settings.compression_type,
                self.settings.compression_level,
            )?),
            OutputFormat::Gff3 => Box::new(Gff3Formatter::new(
                &path,
                &self.settings.compression_type,
                self.settings.compression_level,
            )?),
            OutputFormat::Reprise => {
                return Err(RepriseError::config(
                    "REPrise native format not yet implemented".to_string(),
                ))
            }
        };

        Ok((path, formatter))
    }

    /// Write output files with detected repeats
    pub fn write_results(
        mut self,
        repeats: &[DetectedRepeat],
        genome: &Genome,
    ) -> Result<Vec<FileChecksum>> {
        // Convert DetectedRepeats to RepeatEntries
        let entries: Vec<RepeatEntry> = repeats
            .iter()
            .enumerate()
            .map(|(i, repeat)| self.convert_repeat_to_entry(repeat, i, genome))
            .collect();

        // Write headers
        for formatter in &mut self.formatters {
            if self.settings.include_metadata {
                formatter.write_header(&self.metadata)?;
            }
        }

        // Write entries
        for entry in &entries {
            for formatter in &mut self.formatters {
                formatter.write_entry(entry)?;
            }
        }

        // Write footers
        let total_entries = entries.len() as u64;
        for formatter in &mut self.formatters {
            formatter.write_footer(total_entries)?;
        }

        // Generate checksums if requested
        if self.settings.generate_checksums {
            for path in &self.output_paths {
                let checksum = self.generate_checksum(path)?;
                self.checksums.push(checksum);
            }
            
            // Write checksum file
            self.write_checksum_file()?;
        }

        Ok(self.checksums)
    }

    /// Convert DetectedRepeat to RepeatEntry
    fn convert_repeat_to_entry(
        &self,
        repeat: &DetectedRepeat,
        index: usize,
        genome: &Genome,
    ) -> RepeatEntry {
        let name = format!("repeat_{:06}", index);
        
        // Get contig names (fallback to IDs if names not available)
        let chr1 = format!("contig_{}", repeat.region1.contig);
        let chr2 = format!("contig_{}", repeat.region2.contig);

        // Determine repeat type based on characteristics
        let repeat_type = if repeat.identity > 0.95 {
            Some("exact_repeat".to_string())
        } else if repeat.identity > 0.8 {
            Some("tandem_repeat".to_string())
        } else if repeat.identity > 0.6 {
            Some("interspersed_repeat".to_string())
        } else {
            Some("divergent_repeat".to_string())
        };

        // Additional attributes
        let mut attributes = HashMap::new();
        attributes.insert("seed_kmer".to_string(), format!("{:016x}", repeat.seed_kmer.canonical()));
        attributes.insert("score".to_string(), repeat.score.to_string());
        attributes.insert("length".to_string(), repeat.length.to_string());

        RepeatEntry {
            chr1,
            start1: repeat.region1.start,
            end1: repeat.region1.end,
            chr2,
            start2: repeat.region2.start,
            end2: repeat.region2.end,
            name,
            score: repeat.score,
            identity: repeat.identity,
            length: repeat.length,
            strand1: '+', // TODO: Determine actual strand from alignment
            strand2: '+', // TODO: Determine actual strand from alignment
            repeat_type,
            attributes,
        }
    }

    /// Generate SHA-256 checksum for a file
    fn generate_checksum(&self, path: &Path) -> Result<FileChecksum> {
        use std::fs::File;
        use std::io::Read;
        
        let mut file = File::open(path)
            .map_err(|e| RepriseError::io_error(format!("Failed to open file for checksum: {}", e)))?;
        
        let metadata = file.metadata()
            .map_err(|e| RepriseError::io_error(format!("Failed to get file metadata: {}", e)))?;
        
        let mut hasher = sha2::Sha256::new();
        let mut buffer = vec![0; 8192];
        
        loop {
            let bytes_read = file.read(&mut buffer)
                .map_err(|e| RepriseError::io_error(format!("Failed to read file for checksum: {}", e)))?;
            
            if bytes_read == 0 {
                break;
            }
            
            hasher.update(&buffer[..bytes_read]);
        }
        
        let hash = hasher.finalize();
        let hash_string = format!("{:x}", hash);

        Ok(FileChecksum {
            file_path: path.to_path_buf(),
            size_bytes: metadata.len(),
            sha256: hash_string,
            created: Utc::now(),
        })
    }

    /// Write checksum file
    fn write_checksum_file(&self) -> Result<()> {
        let checksum_filename = format!("{}.checksums", self.settings.file_prefix);
        let checksum_path = if let Some(dir) = &self.settings.output_dir {
            dir.join(checksum_filename)
        } else {
            PathBuf::from(checksum_filename)
        };

        let mut writer = BufWriter::new(File::create(&checksum_path)
            .map_err(|e| RepriseError::io_error(format!("Failed to create checksum file: {}", e)))?);

        writeln!(writer, "# REPrise output file checksums")?;
        writeln!(writer, "# Generated: {}", Utc::now())?;
        writeln!(writer, "# Format: SHA256 *filename")?;
        writeln!(writer)?;

        for checksum in &self.checksums {
            writeln!(writer, "{} *{}", 
                checksum.sha256, 
                checksum.file_path.file_name().unwrap().to_string_lossy())?;
        }

        writer.flush()
            .map_err(|e| RepriseError::io_error(format!("Failed to flush checksum file: {}", e)))?;

        Ok(())
    }
}

// Import sha2 for checksums
use sha2::{Digest, Sha256};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::Kmer;
    use crate::pipeline::{CandidatePair, GenomicRegion};
    use std::io::Write;
    use tempfile::{NamedTempFile, TempDir};

    fn create_test_metadata() -> OutputMetadata {
        OutputMetadata {
            run_id: Uuid::new_v4(),
            timestamp: Utc::now(),
            version: "1.0.0".to_string(),
            command_line: "reprise --input test.fa --output test".to_string(),
            input_file: PathBuf::from("test.fa"),
            parameters: ProcessingParameters {
                k_mer_length: 13,
                min_frequency: 3,
                max_frequency: Some(10000),
                region_extension: 100,
                min_alignment_score: 10,
                min_identity: 0.50,
                num_workers: 4,
            },
            statistics: OutputStatistics {
                total_repeats_found: 100,
                candidates_processed: 1000,
                processing_time_seconds: 60.0,
                genome_coverage: 0.25,
                average_repeat_length: 500.0,
                average_identity: 0.75,
                average_score: 50.0,
            },
            genome_info: GenomeInfo {
                total_length: 1000000,
                num_contigs: 10,
                n50: 100000,
                gc_content: 0.42,
            },
        }
    }

    fn create_test_entry() -> RepeatEntry {
        let mut attributes = HashMap::new();
        attributes.insert("test_attr".to_string(), "test_value".to_string());

        RepeatEntry {
            chr1: "chr1".to_string(),
            start1: 1000,
            end1: 1500,
            chr2: "chr2".to_string(),
            start2: 5000,
            end2: 5500,
            name: "repeat_001".to_string(),
            score: 75,
            identity: 0.85,
            length: 500,
            strand1: '+',
            strand2: '+',
            repeat_type: Some("interspersed".to_string()),
            attributes,
        }
    }

    #[test]
    fn test_tsv_formatter() {
        let temp_file = NamedTempFile::new().unwrap();
        let mut formatter = TsvFormatter::new(
            temp_file.path(),
            &CompressionType::None,
            6,
        ).unwrap();

        let metadata = create_test_metadata();
        let entry = create_test_entry();

        formatter.write_header(&metadata).unwrap();
        formatter.write_entry(&entry).unwrap();
        formatter.write_footer(1).unwrap();

        // Verify the output contains expected content
        let content = std::fs::read_to_string(temp_file.path()).unwrap();
        assert!(content.contains("chr1\tstart1\tend1"));
        assert!(content.contains("chr1\t1000\t1500"));
        assert!(content.contains("REPrise repeat detection results"));
    }

    #[test]
    fn test_json_formatter() {
        let temp_file = NamedTempFile::new().unwrap();
        let mut formatter = JsonFormatter::new(
            temp_file.path(),
            &CompressionType::None,
            6,
        ).unwrap();

        let metadata = create_test_metadata();
        let entry = create_test_entry();

        formatter.write_header(&metadata).unwrap();
        formatter.write_entry(&entry).unwrap();
        formatter.write_footer(1).unwrap();

        // Verify the output is valid JSON
        let content = std::fs::read_to_string(temp_file.path()).unwrap();
        let parsed: serde_json::Value = serde_json::from_str(&content).unwrap();
        assert!(parsed["metadata"].is_object());
        assert!(parsed["repeats"].is_array());
    }

    #[test]
    fn test_bed_formatter() {
        let temp_file = NamedTempFile::new().unwrap();
        let mut formatter = BedFormatter::new(
            temp_file.path(),
            &CompressionType::None,
            6,
        ).unwrap();

        let metadata = create_test_metadata();
        let entry = create_test_entry();

        formatter.write_header(&metadata).unwrap();
        formatter.write_entry(&entry).unwrap();
        formatter.write_footer(2).unwrap(); // BED creates 2 entries per repeat

        let content = std::fs::read_to_string(temp_file.path()).unwrap();
        assert!(content.contains("track name"));
        assert!(content.contains("chr1\t1000\t1500"));
        assert!(content.contains("chr2\t5000\t5500"));
    }

    #[test]
    fn test_gff3_formatter() {
        let temp_file = NamedTempFile::new().unwrap();
        let mut formatter = Gff3Formatter::new(
            temp_file.path(),
            &CompressionType::None,
            6,
        ).unwrap();

        let metadata = create_test_metadata();
        let entry = create_test_entry();

        formatter.write_header(&metadata).unwrap();
        formatter.write_entry(&entry).unwrap();
        formatter.write_footer(1).unwrap();

        let content = std::fs::read_to_string(temp_file.path()).unwrap();
        assert!(content.contains("##gff-version 3"));
        assert!(content.contains("chr1\tREPrise\trepeat_region"));
        assert!(content.contains("ID=repeat_001"));
    }

    #[test]
    fn test_compressed_writer_gzip() {
        let temp_file = NamedTempFile::new().unwrap();
        let mut writer = CompressedWriter::new(
            temp_file.path(),
            &CompressionType::Gzip,
            6,
        ).unwrap();

        let test_data = b"Hello, compressed world!";
        writer.write_all(test_data).unwrap();
        writer.flush().unwrap();

        // The file should be smaller than uncompressed and contain gzip magic bytes
        let compressed_content = std::fs::read(temp_file.path()).unwrap();
        assert!(compressed_content.len() < test_data.len() + 100); // Some overhead is expected
        assert_eq!(&compressed_content[0..2], &[0x1f, 0x8b]); // Gzip magic bytes
    }

    #[test]
    fn test_output_format_display() {
        assert_eq!(OutputFormat::Tsv.to_string(), "tsv");
        assert_eq!(OutputFormat::Bed.to_string(), "bed");
        assert_eq!(OutputFormat::Json.to_string(), "json");
        assert_eq!(OutputFormat::Gff3.to_string(), "gff3");
        assert_eq!(OutputFormat::Reprise.to_string(), "reprof");
    }

    #[test]
    fn test_metadata_serialization() {
        let metadata = create_test_metadata();
        
        // Test JSON serialization
        let json_str = serde_json::to_string(&metadata).unwrap();
        let deserialized: OutputMetadata = serde_json::from_str(&json_str).unwrap();
        assert_eq!(metadata.version, deserialized.version);
        assert_eq!(metadata.parameters.k_mer_length, deserialized.parameters.k_mer_length);
    }
}