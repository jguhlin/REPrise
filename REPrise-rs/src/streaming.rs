use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use memmap2::MmapOptions;
use ahash::AHashMap;
use bloom::{BloomFilter, ASMS};
use crossbeam_channel::{Sender, bounded};
use crate::error::Result;
use crate::genome::ContigInfo;

/// Configuration for streaming k-mer processing
#[derive(Debug, Clone)]
pub struct StreamingConfig {
    pub k: usize,
    pub min_frequency: usize,
    pub max_memory_mb: usize,
    pub temp_dir: PathBuf,
    pub bloom_filter_size: usize,
    pub channel_capacity: usize,
}

impl Default for StreamingConfig {
    fn default() -> Self {
        Self {
            k: 15,  // Safe default that won't explode memory
            min_frequency: 2,
            max_memory_mb: 1024,  // 1GB memory limit
            temp_dir: std::env::temp_dir(),
            bloom_filter_size: 10_000_000,  // 10M expected k-mers
            channel_capacity: 1000,
        }
    }
}

/// Represents a k-mer candidate for repeat detection
#[derive(Debug, Clone, PartialEq)]
pub struct KmerCandidate {
    pub kmer: u64,  // Canonical k-mer representation
    pub frequency: usize,
    pub positions: Vec<(u32, u16)>,  // (global_pos, contig_id)
    pub is_reverse: Vec<bool>,
}

/// Disk-based k-mer index using memory-mapped files
pub struct StreamingKmerIndex {
    config: StreamingConfig,
    bloom_filter: BloomFilter,
    temp_files: Vec<PathBuf>,
    kmer_count: usize,
    current_contig_id: u16,
    contig_info: Vec<ContigInfo>,
    
    // Memory-bounded in-memory buffer
    current_batch: AHashMap<u64, Vec<(u32, u16, bool)>>,
    batch_size_mb: f64,
}

impl StreamingKmerIndex {
    pub fn new(config: StreamingConfig) -> Result<Self> {
        let bloom_filter = BloomFilter::with_rate(0.01, config.bloom_filter_size as u32);
        
        Ok(Self {
            config,
            bloom_filter,
            temp_files: Vec::new(),
            kmer_count: 0,
            current_contig_id: 0,
            contig_info: Vec::new(),
            current_batch: AHashMap::new(),
            batch_size_mb: 0.0,
        })
    }
    
    /// Process FASTA file in streaming fashion, one contig at a time
    pub fn build_from_fasta<P: AsRef<Path>>(&mut self, fasta_path: P) -> Result<()> {
        let mut reader = needletail::parse_fastx_file(&fasta_path)?;
        let mut global_position = 0u32;
        
        while let Some(record) = reader.next() {
            let record = record?;
            let contig_name = String::from_utf8_lossy(record.id()).to_string();
            let sequence = record.seq();
            
            println!("Processing contig: {} ({} bp)", contig_name, sequence.len());
            
            // Store contig info
            self.contig_info.push(ContigInfo {
                name: contig_name,
                start: global_position as u64,
                end: (global_position + sequence.len() as u32) as u64,
                length: sequence.len() as u64,
            });
            
            // Process this contig's k-mers
            self.process_contig_kmers(&sequence, global_position)?;
            
            global_position += sequence.len() as u32;
            self.current_contig_id += 1;
            
            // Flush batch if memory limit reached
            if self.batch_size_mb > self.config.max_memory_mb as f64 * 0.8 {
                self.flush_batch_to_disk()?;
            }
        }
        
        // Final flush
        self.flush_batch_to_disk()?;
        
        println!("Processed {} k-mers across {} contigs", self.kmer_count, self.current_contig_id);
        Ok(())
    }
    
    /// Extract k-mers from a single contig
    fn process_contig_kmers(&mut self, sequence: &[u8], start_pos: u32) -> Result<()> {
        if sequence.len() < self.config.k {
            return Ok(());
        }
        
        for i in 0..=sequence.len() - self.config.k {
            let kmer_slice = &sequence[i..i + self.config.k];
            
            // Skip k-mers with N bases
            if kmer_slice.iter().any(|&b| b == b'N' || b == b'n') {
                continue;
            }
            
            let (canonical_kmer, is_reverse) = self.canonicalize_kmer(kmer_slice);
            let global_pos = start_pos + i as u32;
            
            // Add to bloom filter
            self.bloom_filter.insert(&canonical_kmer.to_le_bytes());
            
            // Add to current batch
            self.current_batch
                .entry(canonical_kmer)
                .or_insert_with(Vec::new)
                .push((global_pos, self.current_contig_id, is_reverse));
            
            self.kmer_count += 1;
            
            // Update memory usage estimate
            if self.kmer_count % 1000 == 0 {
                self.update_memory_estimate();
            }
        }
        
        Ok(())
    }
    
    /// Convert k-mer sequence to canonical form (lexicographically smaller of forward/reverse)
    fn canonicalize_kmer(&self, kmer: &[u8]) -> (u64, bool) {
        let forward = self.encode_kmer(kmer);
        let reverse = self.encode_kmer(&self.reverse_complement(kmer));
        
        if forward <= reverse {
            (forward, false)
        } else {
            (reverse, true)
        }
    }
    
    /// Encode DNA sequence as 2-bit packed integer
    fn encode_kmer(&self, kmer: &[u8]) -> u64 {
        let mut encoded = 0u64;
        for &base in kmer {
            encoded = (encoded << 2) | match base.to_ascii_uppercase() {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 0,  // Treat unknown as A
            };
        }
        encoded
    }
    
    /// Get reverse complement of sequence
    fn reverse_complement(&self, seq: &[u8]) -> Vec<u8> {
        seq.iter().rev().map(|&base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => base,
        }).collect()
    }
    
    /// Estimate current memory usage
    fn update_memory_estimate(&mut self) {
        let entries = self.current_batch.len();
        let avg_positions = if entries > 0 {
            self.current_batch.values().map(|v| v.len()).sum::<usize>() / entries
        } else { 0 };
        
        // Rough estimate: HashMap overhead + key + Vec<(u32,u16,bool)>
        self.batch_size_mb = (entries * (8 + 8 + avg_positions * 7)) as f64 / (1024.0 * 1024.0);
    }
    
    /// Flush current batch to disk-based storage
    fn flush_batch_to_disk(&mut self) -> Result<()> {
        if self.current_batch.is_empty() {
            return Ok(());
        }
        
        let temp_file_path = self.config.temp_dir.join(format!("kmers_batch_{}.tmp", self.temp_files.len()));
        let mut file = BufWriter::new(File::create(&temp_file_path)?);
        
        println!("Flushing {} k-mers to disk ({:.2} MB)", self.current_batch.len(), self.batch_size_mb);
        
        // Write batch data in binary format for efficiency
        for (&kmer, positions) in &self.current_batch {
            // Filter by minimum frequency
            if positions.len() < self.config.min_frequency {
                continue;
            }
            
            // Write: kmer(8) + freq(4) + positions
            file.write_all(&kmer.to_le_bytes())?;
            file.write_all(&(positions.len() as u32).to_le_bytes())?;
            
            for &(pos, contig_id, is_rev) in positions {
                file.write_all(&pos.to_le_bytes())?;
                file.write_all(&contig_id.to_le_bytes())?;
                file.write_all(&[if is_rev { 1u8 } else { 0u8 }])?;
            }
        }
        
        file.flush()?;
        drop(file);
        
        self.temp_files.push(temp_file_path);
        self.current_batch.clear();
        self.batch_size_mb = 0.0;
        
        Ok(())
    }
    
    /// Check if k-mer likely exists using Bloom filter
    pub fn may_contain_kmer(&self, kmer: u64) -> bool {
        self.bloom_filter.contains(&kmer.to_le_bytes())
    }
    
    /// Get candidates from disk-based storage for processing pipeline
    pub fn candidate_producer(&self, sender: Sender<KmerCandidate>) -> Result<()> {
        for temp_file in &self.temp_files {
            let file = File::open(temp_file)?;
            let mmap = unsafe { MmapOptions::new().map(&file)? };
            
            let mut cursor = 0;
            while cursor < mmap.len() {
                if cursor + 12 > mmap.len() { break; }
                
                // Read k-mer and frequency
                let kmer = u64::from_le_bytes([
                    mmap[cursor], mmap[cursor+1], mmap[cursor+2], mmap[cursor+3],
                    mmap[cursor+4], mmap[cursor+5], mmap[cursor+6], mmap[cursor+7]
                ]);
                cursor += 8;
                
                let freq = u32::from_le_bytes([
                    mmap[cursor], mmap[cursor+1], mmap[cursor+2], mmap[cursor+3]
                ]) as usize;
                cursor += 4;
                
                // Read positions
                let mut positions = Vec::new();
                let mut is_reverse = Vec::new();
                
                for _ in 0..freq {
                    if cursor + 7 > mmap.len() { break; }
                    
                    let pos = u32::from_le_bytes([
                        mmap[cursor], mmap[cursor+1], mmap[cursor+2], mmap[cursor+3]
                    ]);
                    cursor += 4;
                    
                    let contig_id = u16::from_le_bytes([mmap[cursor], mmap[cursor+1]]);
                    cursor += 2;
                    
                    let is_rev = mmap[cursor] != 0;
                    cursor += 1;
                    
                    positions.push((pos, contig_id));
                    is_reverse.push(is_rev);
                }
                
                if positions.len() >= self.config.min_frequency {
                    let candidate = KmerCandidate {
                        kmer,
                        frequency: freq,
                        positions,
                        is_reverse,
                    };
                    
                    if sender.send(candidate).is_err() {
                        break; // Consumer dropped
                    }
                }
            }
        }
        
        Ok(())
    }
    
    /// Clean up temporary files
    pub fn cleanup(&mut self) -> Result<()> {
        for temp_file in &self.temp_files {
            if temp_file.exists() {
                std::fs::remove_file(temp_file)?;
            }
        }
        self.temp_files.clear();
        Ok(())
    }
}

impl Drop for StreamingKmerIndex {
    fn drop(&mut self) {
        let _ = self.cleanup();
    }
}

/// Producer-consumer pipeline for memory-bounded repeat detection
pub struct StreamingPipeline {
    config: StreamingConfig,
}

impl StreamingPipeline {
    pub fn new(config: StreamingConfig) -> Self {
        Self { config }
    }
    
    /// Run streaming pipeline with bounded memory usage
    pub fn run<P: AsRef<Path>>(&self, fasta_path: P) -> Result<Vec<KmerCandidate>> {
        // Build streaming index
        let mut index = StreamingKmerIndex::new(self.config.clone())?;
        index.build_from_fasta(&fasta_path)?;
        
        // Set up producer-consumer channels
        let (sender, receiver) = bounded::<KmerCandidate>(self.config.channel_capacity);
        
        // Spawn producer thread
        let index_for_producer = index;
        let producer_handle = std::thread::spawn(move || {
            index_for_producer.candidate_producer(sender)
        });
        
        // Consumer: collect candidates for further processing
        let mut candidates = Vec::new();
        while let Ok(candidate) = receiver.recv() {
            candidates.push(candidate);
        }
        
        producer_handle.join().map_err(|_| crate::error::REPriseError::ThreadError("Producer thread failed".into()))?;
        
        println!("Pipeline processed {} candidate k-mers", candidates.len());
        Ok(candidates)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;
    
    #[test]
    fn test_kmer_encoding() {
        let config = StreamingConfig::default();
        let index = StreamingKmerIndex::new(config).unwrap();
        
        let kmer = b"ACGT";
        let encoded = index.encode_kmer(kmer);
        assert_eq!(encoded, 0b00011011); // A=00, C=01, G=10, T=11
    }
    
    #[test]
    fn test_canonical_kmer() {
        let config = StreamingConfig::default();
        let index = StreamingKmerIndex::new(config).unwrap();
        
        let kmer = b"ACGT";
        let rev_comp = b"ACGT"; // ACGT reverse complement is ACGT (palindrome)
        
        let (canonical1, is_rev1) = index.canonicalize_kmer(kmer);
        let (canonical2, is_rev2) = index.canonicalize_kmer(rev_comp);
        
        assert_eq!(canonical1, canonical2);
    }
    
    #[test]
    fn test_streaming_small_fasta() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        writeln!(temp_file, ">test_contig")?;
        writeln!(temp_file, "ACGTACGTACGTACGTACGTACGT")?;
        temp_file.flush()?;
        
        let mut config = StreamingConfig::default();
        config.k = 4;
        config.min_frequency = 1;
        
        let pipeline = StreamingPipeline::new(config);
        let candidates = pipeline.run(temp_file.path())?;
        
        assert!(!candidates.is_empty());
        Ok(())
    }
}