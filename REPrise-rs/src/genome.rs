//! High-performance genome loading and management
//!
//! This module provides efficient loading and access to large genomes with
//! support for highly fragmented assemblies (250k+ contigs). It uses the
//! needletail crate for fast FASTA parsing and maintains contig boundaries
//! for safe concurrent processing.

use crate::error::{REPriseError, Result};
use needletail::{parse_fastx_file, Sequence};
use std::collections::HashMap;
use std::ops::Range;
use std::path::Path;

/// Contig identifier (supports up to 4 billion contigs)
pub type ContigId = u32;

/// Maximum supported genome size (16 TiB)
pub const MAX_GENOME_SIZE: u64 = 16 * 1024 * 1024 * 1024 * 1024;

/// Information about a single contig
#[derive(Debug, Clone)]
pub struct ContigInfo {
    /// Human-readable contig name
    pub name: String,
    /// Global start position in concatenated genome
    pub start: u64,
    /// Global end position (exclusive) in concatenated genome
    pub end: u64,
    /// Length of this contig
    pub length: u64,
}

impl ContigInfo {
    /// Returns the range of this contig in the global genome
    pub fn range(&self) -> Range<u64> {
        self.start..self.end
    }
    
    /// Converts a local position (0-based, within contig) to global position
    pub fn local_to_global(&self, local_pos: u64) -> Option<u64> {
        if local_pos < self.length {
            Some(self.start + local_pos)
        } else {
            None
        }
    }
    
    /// Converts a global position to local position within this contig
    pub fn global_to_local(&self, global_pos: u64) -> Option<u64> {
        if global_pos >= self.start && global_pos < self.end {
            Some(global_pos - self.start)
        } else {
            None
        }
    }
}

/// High-performance genome representation with contig awareness
///
/// The Genome struct stores the entire genome as a concatenated sequence while
/// maintaining contig boundary information for safe parallel processing. It
/// supports efficient random access and contig-aware operations.
pub struct Genome {
    /// Concatenated genome sequence (all contigs joined)
    sequence: Vec<u8>,
    /// Contig information indexed by ContigId
    contigs: Vec<ContigInfo>,
    /// Fast lookup from contig name to ContigId
    name_to_id: HashMap<String, ContigId>,
    /// Total genome size in bases
    total_size: u64,
}

impl Genome {
    /// Builds a Genome from a FASTA file using high-performance needletail parser
    ///
    /// # Arguments
    /// * `path` - Path to the FASTA file
    ///
    /// # Returns
    /// * `Ok(Genome)` - Successfully loaded genome
    /// * `Err(REPriseError)` - File I/O error, invalid FASTA, or genome too large
    ///
    /// # Performance
    /// This function is optimized for large genomes and can handle files with
    /// hundreds of thousands of contigs efficiently.
    pub fn from_fasta<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        
        // Pre-allocate vectors with reasonable capacity
        let mut sequence = Vec::with_capacity(100_000_000); // Start with 100MB
        let mut contigs = Vec::with_capacity(10_000); // Support many contigs
        let mut name_to_id = HashMap::with_capacity(10_000);
        
        // Use needletail for high-performance FASTA parsing
        let mut reader = parse_fastx_file(path)
            .map_err(|e| REPriseError::Io(std::io::Error::new(std::io::ErrorKind::Other, e)))?;
        
        let mut current_pos = 0u64;
        let mut contig_id = 0u32;
        
        while let Some(record) = reader.next() {
            let record = record
                .map_err(|e| REPriseError::Io(std::io::Error::new(std::io::ErrorKind::InvalidData, e)))?;
            
            // Extract contig name and clean it
            let raw_name = String::from_utf8_lossy(record.id());
            let name = clean_contig_name(&raw_name);
            
            // Get sequence and validate
            let seq = record.normalize(false); // Don't reverse complement
            let seq_len = seq.len() as u64;
            
            if seq_len == 0 {
                continue; // Skip empty sequences
            }
            
            // Check for genome size overflow
            if current_pos.saturating_add(seq_len) > MAX_GENOME_SIZE {
                return Err(REPriseError::GenomeTooLarge(
                    current_pos + seq_len,
                    MAX_GENOME_SIZE,
                ));
            }
            
            // Convert sequence to numeric format (A=0, C=1, G=2, T=3, others=99)
            let numeric_seq: Vec<u8> = seq.iter()
                .map(|&b| nucleotide_to_numeric(b))
                .collect();
            
            // Store contig information
            let contig_info = ContigInfo {
                name: name.clone(),
                start: current_pos,
                end: current_pos + seq_len,
                length: seq_len,
            };
            
            // Check for duplicate contig names
            if name_to_id.contains_key(&name) {
                return Err(REPriseError::invalid_fasta(
                    contig_id as usize + 1,
                    format!("Duplicate contig name: {}", name),
                ));
            }
            
            // Add to collections
            contigs.push(contig_info);
            name_to_id.insert(name, contig_id);
            sequence.extend(numeric_seq);
            
            current_pos += seq_len;
            contig_id += 1;
            
            // Safety check for maximum number of contigs
            if contig_id == u32::MAX {
                return Err(REPriseError::Config(
                    "Maximum number of contigs (4 billion) exceeded".to_string(),
                ));
            }
        }
        
        if contigs.is_empty() {
            return Err(REPriseError::invalid_fasta(0, "No valid sequences found in FASTA file"));
        }
        
        // Shrink collections to actual size to free unused memory
        sequence.shrink_to_fit();
        contigs.shrink_to_fit();
        name_to_id.shrink_to_fit();
        
        Ok(Genome {
            sequence,
            contigs,
            name_to_id,
            total_size: current_pos,
        })
    }
    
    /// Get a slice of the genome sequence
    ///
    /// # Arguments
    /// * `range` - Range of positions to extract (global coordinates)
    ///
    /// # Returns
    /// Reference to the sequence slice, or empty slice if range is invalid
    ///
    /// # Performance
    /// This is a zero-copy operation that returns a reference to the underlying data.
    pub fn slice(&self, range: Range<u64>) -> &[u8] {
        let start = range.start as usize;
        let end = range.end.min(self.total_size) as usize;
        
        if start >= self.sequence.len() || start >= end {
            return &[];
        }
        
        &self.sequence[start..end.min(self.sequence.len())]
    }
    
    /// Find which contig contains a given global position
    ///
    /// # Arguments
    /// * `pos` - Global position in the genome
    ///
    /// # Returns
    /// * `Some(ContigId)` - ID of the contig containing this position
    /// * `None` - Position is out of bounds
    pub fn contig_of(&self, pos: u64) -> Option<ContigId> {
        // Binary search for the contig containing this position
        match self.contigs.binary_search_by(|contig| {
            if pos < contig.start {
                std::cmp::Ordering::Greater
            } else if pos >= contig.end {
                std::cmp::Ordering::Less
            } else {
                std::cmp::Ordering::Equal
            }
        }) {
            Ok(index) => Some(index as ContigId),
            Err(_) => None,
        }
    }
    
    /// Get the range of positions occupied by a contig
    ///
    /// # Arguments
    /// * `id` - ContigId to query
    ///
    /// # Returns
    /// * `Some(Range<u64>)` - Global range of the contig
    /// * `None` - Invalid contig ID
    pub fn contig_range(&self, id: ContigId) -> Option<Range<u64>> {
        self.contigs.get(id as usize).map(|info| info.range())
    }
    
    /// Check if a region is entirely within a single contig
    ///
    /// # Arguments
    /// * `start` - Start position (global)
    /// * `len` - Length of the region
    ///
    /// # Returns
    /// `true` if the entire region [start, start+len) is within one contig
    pub fn is_within_one_contig(&self, start: u64, len: usize) -> bool {
        let end = start + len as u64;
        
        if let Some(contig_id) = self.contig_of(start) {
            if let Some(contig_range) = self.contig_range(contig_id) {
                return end <= contig_range.end;
            }
        }
        
        false
    }
    
    /// Get contig information by ID
    pub fn contig_info(&self, id: ContigId) -> Option<&ContigInfo> {
        self.contigs.get(id as usize)
    }
    
    /// Get contig ID by name
    pub fn contig_id(&self, name: &str) -> Option<ContigId> {
        self.name_to_id.get(name).copied()
    }
    
    /// Get total genome size in bases
    pub fn len(&self) -> u64 {
        self.total_size
    }
    
    /// Check if genome is empty
    pub fn is_empty(&self) -> bool {
        self.total_size == 0
    }
    
    /// Get number of contigs
    pub fn num_contigs(&self) -> usize {
        self.contigs.len()
    }
    
    /// Iterator over all contigs
    pub fn contigs(&self) -> impl Iterator<Item = (ContigId, &ContigInfo)> {
        self.contigs
            .iter()
            .enumerate()
            .map(|(id, info)| (id as ContigId, info))
    }
}

/// Clean contig name for consistent indexing
fn clean_contig_name(raw_name: &str) -> String {
    raw_name
        .trim()
        .chars()
        .map(|c| if c.is_whitespace() { '_' } else { c })
        .collect()
}

/// Convert nucleotide character to numeric code
/// Matches the encoding used by the original C++ implementation
fn nucleotide_to_numeric(nucleotide: u8) -> u8 {
    match nucleotide.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 99, // All ambiguous/unknown bases
    }
}

/// Convert numeric code back to nucleotide character
pub fn numeric_to_nucleotide(code: u8) -> char {
    match code {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        _ => 'N',
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;
    
    fn create_test_fasta(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();
        file.flush().unwrap();
        file
    }
    
    #[test]
    fn test_simple_genome_loading() {
        let fasta_content = ">contig1\nATCG\n>contig2\nGCTA\n";
        let file = create_test_fasta(fasta_content);
        
        let genome = Genome::from_fasta(file.path()).unwrap();
        
        assert_eq!(genome.len(), 8);
        assert_eq!(genome.num_contigs(), 2);
        
        // Check sequence encoding
        assert_eq!(genome.slice(0..4), &[0, 3, 1, 2]); // ATCG
        assert_eq!(genome.slice(4..8), &[2, 1, 3, 0]); // GCTA
    }
    
    #[test]
    fn test_contig_boundaries() {
        let fasta_content = ">chr1\nAAAAA\n>chr2\nTTTTT\n";
        let file = create_test_fasta(fasta_content);
        
        let genome = Genome::from_fasta(file.path()).unwrap();
        
        // Test contig identification
        assert_eq!(genome.contig_of(0), Some(0));
        assert_eq!(genome.contig_of(4), Some(0));
        assert_eq!(genome.contig_of(5), Some(1));
        assert_eq!(genome.contig_of(9), Some(1));
        assert_eq!(genome.contig_of(10), None);
        
        // Test contig ranges
        assert_eq!(genome.contig_range(0), Some(0..5));
        assert_eq!(genome.contig_range(1), Some(5..10));
        assert_eq!(genome.contig_range(2), None);
    }
    
    #[test]
    fn test_within_one_contig() {
        let fasta_content = ">chr1\nAAAAA\n>chr2\nTTTTT\n";
        let file = create_test_fasta(fasta_content);
        
        let genome = Genome::from_fasta(file.path()).unwrap();
        
        // Within single contig
        assert!(genome.is_within_one_contig(0, 5));
        assert!(genome.is_within_one_contig(5, 5));
        assert!(genome.is_within_one_contig(1, 3));
        
        // Crosses contig boundary
        assert!(!genome.is_within_one_contig(3, 5));
        assert!(!genome.is_within_one_contig(0, 10));
        
        // Out of bounds
        assert!(!genome.is_within_one_contig(10, 1));
    }
    
    #[test]
    fn test_contig_info() {
        let fasta_content = ">test_contig\nATCGATCG\n";
        let file = create_test_fasta(fasta_content);
        
        let genome = Genome::from_fasta(file.path()).unwrap();
        
        let info = genome.contig_info(0).unwrap();
        assert_eq!(info.name, "test_contig");
        assert_eq!(info.start, 0);
        assert_eq!(info.end, 8);
        assert_eq!(info.length, 8);
        assert_eq!(info.range(), 0..8);
        
        // Test coordinate conversion
        assert_eq!(info.local_to_global(0), Some(0));
        assert_eq!(info.local_to_global(7), Some(7));
        assert_eq!(info.local_to_global(8), None);
        
        assert_eq!(info.global_to_local(0), Some(0));
        assert_eq!(info.global_to_local(7), Some(7));
        assert_eq!(info.global_to_local(8), None);
    }
    
    #[test]
    fn test_name_lookup() {
        let fasta_content = ">contig_A\nATCG\n>contig_B\nGCTA\n";
        let file = create_test_fasta(fasta_content);
        
        let genome = Genome::from_fasta(file.path()).unwrap();
        
        assert_eq!(genome.contig_id("contig_A"), Some(0));
        assert_eq!(genome.contig_id("contig_B"), Some(1));
        assert_eq!(genome.contig_id("nonexistent"), None);
    }
    
    #[test]
    fn test_nucleotide_conversion() {
        assert_eq!(nucleotide_to_numeric(b'A'), 0);
        assert_eq!(nucleotide_to_numeric(b'a'), 0);
        assert_eq!(nucleotide_to_numeric(b'C'), 1);
        assert_eq!(nucleotide_to_numeric(b'G'), 2);
        assert_eq!(nucleotide_to_numeric(b'T'), 3);
        assert_eq!(nucleotide_to_numeric(b'N'), 99);
        assert_eq!(nucleotide_to_numeric(b'X'), 99);
        
        assert_eq!(numeric_to_nucleotide(0), 'A');
        assert_eq!(numeric_to_nucleotide(1), 'C');
        assert_eq!(numeric_to_nucleotide(2), 'G');
        assert_eq!(numeric_to_nucleotide(3), 'T');
        assert_eq!(numeric_to_nucleotide(99), 'N');
    }
    
    #[test]
    fn test_empty_sequences() {
        let fasta_content = ">empty1\n\n>valid\nATCG\n>empty2\n\n";
        let file = create_test_fasta(fasta_content);
        
        let genome = Genome::from_fasta(file.path()).unwrap();
        
        // Should only have one contig (empty ones skipped)
        assert_eq!(genome.num_contigs(), 1);
        assert_eq!(genome.contig_info(0).unwrap().name, "valid");
    }
    
    #[test]
    fn test_error_handling() {
        // Test non-existent file
        let result = Genome::from_fasta("/nonexistent/file.fasta");
        assert!(result.is_err());
        
        // Test duplicate contig names
        let fasta_content = ">dup\nATCG\n>dup\nGCTA\n";
        let file = create_test_fasta(fasta_content);
        let result = Genome::from_fasta(file.path());
        assert!(result.is_err());
        
        // Test completely empty file
        let file = create_test_fasta("");
        let result = Genome::from_fasta(file.path());
        assert!(result.is_err());
    }
}