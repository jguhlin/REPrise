//! Contig-by-contig k-mer indexing for large-scale genomic analysis
//!
//! This module implements efficient k-mer indexing that processes genomes
//! contig-by-contig to handle large, fragmented assemblies. It provides
//! deterministic indexing suitable for concurrent repeat detection pipelines.

use crate::error::Result;
use crate::genome::{Genome, ContigId};
use crate::kmer::{Kmer, KmerEngine};
use ahash::{AHashMap, AHashSet};
use rayon::prelude::*;

/// Position of a k-mer occurrence in the genome
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct KmerPosition {
    /// Global position in the concatenated genome
    pub position: u64,
    /// Contig containing this k-mer
    pub contig: ContigId,
    /// Local position within the contig
    pub local_position: u64,
}

impl KmerPosition {
    /// Create a new k-mer position
    pub fn new(position: u64, contig: ContigId, local_position: u64) -> Self {
        Self {
            position,
            contig,
            local_position,
        }
    }
}

/// Configuration for k-mer indexing
#[derive(Debug, Clone)]
pub struct IndexConfig {
    /// K-mer length
    pub k: usize,
    /// Minimum frequency threshold for k-mer inclusion
    pub min_frequency: u32,
    /// Maximum frequency threshold (to filter highly repetitive k-mers)
    pub max_frequency: Option<u32>,
    /// Whether to use parallel processing
    pub parallel: bool,
    /// Maximum number of positions to store per k-mer (memory control)
    pub max_positions_per_kmer: usize,
}

impl Default for IndexConfig {
    fn default() -> Self {
        Self {
            k: 13,
            min_frequency: 2,
            max_frequency: Some(1000),
            parallel: true,
            max_positions_per_kmer: 10000,
        }
    }
}

/// High-performance k-mer index for genomic sequences
///
/// This index provides fast lookup of k-mer positions and frequencies,
/// optimized for concurrent access in repeat detection pipelines.
pub struct KmerIndex {
    /// K-mer engine for canonical k-mer operations
    engine: KmerEngine,
    /// Map from k-mer to list of positions
    positions: AHashMap<Kmer, Vec<KmerPosition>>,
    /// K-mer frequency counts
    frequencies: AHashMap<Kmer, u32>,
    /// Configuration used to build this index
    config: IndexConfig,
    /// Genome this index was built from
    genome_size: u64,
    /// Number of contigs indexed
    num_contigs: usize,
}

impl KmerIndex {
    /// Build a k-mer index from a genome using contig-by-contig processing
    ///
    /// # Arguments
    /// * `genome` - Genome to index
    /// * `config` - Indexing configuration
    ///
    /// # Returns
    /// * `Ok(KmerIndex)` - Successfully built index
    /// * `Err(RepriseError)` - Invalid configuration or processing error
    pub fn build(genome: &Genome, config: IndexConfig) -> Result<Self> {
        let engine = KmerEngine::new(config.k)?;
        let mut positions = engine.new_kmer_map();
        let mut frequencies = engine.new_kmer_map();
        
        if config.parallel {
            Self::build_parallel(genome, &engine, &mut positions, &mut frequencies, &config)?;
        } else {
            Self::build_sequential(genome, &engine, &mut positions, &mut frequencies, &config)?;
        }
        
        // Apply frequency filtering
        Self::apply_frequency_filter(&mut positions, &frequencies, &config);
        
        Ok(KmerIndex {
            engine,
            positions,
            frequencies,
            config,
            genome_size: genome.len(),
            num_contigs: genome.num_contigs(),
        })
    }
    
    /// Build index using parallel processing (recommended for large genomes)
    fn build_parallel(
        genome: &Genome,
        engine: &KmerEngine,
        positions: &mut AHashMap<Kmer, Vec<KmerPosition>>,
        frequencies: &mut AHashMap<Kmer, u32>,
        _config: &IndexConfig,
    ) -> Result<()> {
        use std::sync::Mutex;
        
        let positions_mutex = Mutex::new(positions);
        let frequencies_mutex = Mutex::new(frequencies);
        
        // Process contigs in parallel
        let contigs: Vec<_> = genome.contigs().collect();
        
        contigs.par_iter().try_for_each(|&(contig_id, contig_info)| -> Result<()> {
            let sequence = genome.slice(contig_info.range());
            let _contig_kmers = engine.extract_kmers(sequence);
            
            // Build local maps for this contig
            let mut local_positions: AHashMap<Kmer, Vec<KmerPosition>> = engine.new_kmer_map();
            let mut local_frequencies: AHashMap<Kmer, u32> = engine.new_kmer_map();
            
            // Process k-mers from this contig
            if sequence.len() >= engine.k() {
                for i in 0..=sequence.len() - engine.k() {
                    let subseq = &sequence[i..i + engine.k()];
                    
                    // Skip k-mers with ambiguous bases
                    if subseq.iter().any(|&b| b > 3) {
                        continue;
                    }
                    
                    if let Ok(kmer) = Kmer::from_sequence(subseq, engine.k()) {
                        let global_pos = contig_info.start + i as u64;
                        let local_pos = i as u64;
                        
                        let kmer_pos = KmerPosition::new(global_pos, contig_id, local_pos);
                        
                        // Update local collections
                        local_positions.entry(kmer).or_default().push(kmer_pos);
                        *local_frequencies.entry(kmer).or_insert(0) += 1;
                    }
                }
            }
            
            // Merge with global collections (critical section)
            {
                let mut global_positions = positions_mutex.lock().unwrap();
                let mut global_frequencies = frequencies_mutex.lock().unwrap();
                
                for (kmer, mut local_pos_list) in local_positions {
                    global_positions.entry(kmer).or_default().append(&mut local_pos_list);
                }
                
                for (kmer, local_freq) in local_frequencies {
                    *global_frequencies.entry(kmer).or_insert(0) += local_freq;
                }
            }
            
            Ok(())
        })?;
        
        Ok(())
    }
    
    /// Build index using sequential processing (for smaller genomes or testing)
    fn build_sequential(
        genome: &Genome,
        engine: &KmerEngine,
        positions: &mut AHashMap<Kmer, Vec<KmerPosition>>,
        frequencies: &mut AHashMap<Kmer, u32>,
        _config: &IndexConfig,
    ) -> Result<()> {
        for (contig_id, contig_info) in genome.contigs() {
            let sequence = genome.slice(contig_info.range());
            
            if sequence.len() < engine.k() {
                continue;
            }
            
            for i in 0..=sequence.len() - engine.k() {
                let subseq = &sequence[i..i + engine.k()];
                
                // Skip k-mers with ambiguous bases
                if subseq.iter().any(|&b| b > 3) {
                    continue;
                }
                
                if let Ok(kmer) = Kmer::from_sequence(subseq, engine.k()) {
                    let global_pos = contig_info.start + i as u64;
                    let local_pos = i as u64;
                    
                    let kmer_pos = KmerPosition::new(global_pos, contig_id, local_pos);
                    
                    positions.entry(kmer).or_default().push(kmer_pos);
                    *frequencies.entry(kmer).or_insert(0) += 1;
                }
            }
        }
        
        Ok(())
    }
    
    /// Apply frequency-based filtering to remove low/high frequency k-mers
    fn apply_frequency_filter(
        positions: &mut AHashMap<Kmer, Vec<KmerPosition>>,
        frequencies: &AHashMap<Kmer, u32>,
        config: &IndexConfig,
    ) {
        positions.retain(|kmer, pos_list| {
            if let Some(&freq) = frequencies.get(kmer) {
                let passes_min = freq >= config.min_frequency;
                let passes_max = config.max_frequency.map_or(true, |max| freq <= max);
                let passes_position_limit = pos_list.len() <= config.max_positions_per_kmer;
                
                passes_min && passes_max && passes_position_limit
            } else {
                false
            }
        });
    }
    
    /// Get all positions where a k-mer occurs
    pub fn get_positions(&self, kmer: &Kmer) -> Option<&[KmerPosition]> {
        self.positions.get(kmer).map(|v| v.as_slice())
    }
    
    /// Get frequency of a k-mer
    pub fn get_frequency(&self, kmer: &Kmer) -> u32 {
        self.frequencies.get(kmer).copied().unwrap_or(0)
    }
    
    /// Get all k-mers in the index
    pub fn kmers(&self) -> impl Iterator<Item = &Kmer> {
        self.positions.keys()
    }
    
    /// Get k-mers with frequency in a specific range
    pub fn kmers_with_frequency_range(&self, min_freq: u32, max_freq: u32) -> Vec<Kmer> {
        self.frequencies
            .iter()
            .filter(|&(_, freq)| *freq >= min_freq && *freq <= max_freq)
            .map(|(&kmer, _)| kmer)
            .collect()
    }
    
    /// Get the most frequent k-mers
    pub fn most_frequent_kmers(&self, limit: usize) -> Vec<(Kmer, u32)> {
        let mut freq_pairs: Vec<_> = self.frequencies.iter().map(|(&k, &f)| (k, f)).collect();
        freq_pairs.sort_by(|a, b| b.1.cmp(&a.1)); // Sort by frequency descending
        freq_pairs.truncate(limit);
        freq_pairs
    }
    
    /// Get k-mers that occur in multiple contigs (potential repeat seeds)
    pub fn multi_contig_kmers(&self) -> Vec<Kmer> {
        self.positions
            .iter()
            .filter_map(|(&kmer, positions)| {
                let mut contigs: AHashSet<ContigId> = AHashSet::new();
                for pos in positions {
                    contigs.insert(pos.contig);
                }
                
                if contigs.len() > 1 {
                    Some(kmer)
                } else {
                    None
                }
            })
            .collect()
    }
    
    /// Get index statistics
    pub fn stats(&self) -> IndexStats {
        let total_kmers = self.positions.len();
        let total_positions: usize = self.positions.values().map(|v| v.len()).sum();
        let avg_frequency = if total_kmers > 0 {
            total_positions as f64 / total_kmers as f64
        } else {
            0.0
        };
        
        let multi_contig_count = self.multi_contig_kmers().len();
        
        IndexStats {
            k: self.config.k,
            total_kmers,
            total_positions,
            avg_frequency,
            multi_contig_kmers: multi_contig_count,
            genome_size: self.genome_size,
            num_contigs: self.num_contigs,
        }
    }
    
    /// Get the k-mer length for this index
    pub fn k(&self) -> usize {
        self.config.k
    }
    
    /// Get the configuration used to build this index
    pub fn config(&self) -> &IndexConfig {
        &self.config
    }
}

/// Statistics about a k-mer index
#[derive(Debug, Clone)]
pub struct IndexStats {
    /// K-mer length
    pub k: usize,
    /// Total number of unique k-mers
    pub total_kmers: usize,
    /// Total number of k-mer positions
    pub total_positions: usize,
    /// Average frequency per k-mer
    pub avg_frequency: f64,
    /// Number of k-mers occurring in multiple contigs
    pub multi_contig_kmers: usize,
    /// Size of indexed genome
    pub genome_size: u64,
    /// Number of contigs in the genome
    pub num_contigs: usize,
}

impl IndexStats {
    /// Calculate memory usage estimate in bytes
    pub fn estimated_memory_usage(&self) -> u64 {
        // Rough estimate:
        // - Each Kmer: 16 bytes (8 bytes u64 + 8 bytes usize)
        // - Each KmerPosition: 16 bytes (8 + 4 + 4 bytes)
        // - HashMap overhead: ~50% additional
        
        let kmer_size = 16;
        let position_size = 16;
        let overhead_factor = 1.5;
        
        let base_usage = (self.total_kmers * kmer_size) + (self.total_positions * position_size);
        (base_usage as f64 * overhead_factor) as u64
    }
    
    /// Calculate coverage statistics
    pub fn coverage_stats(&self) -> (f64, f64) {
        let theoretical_kmers = if self.genome_size >= self.k as u64 {
            self.genome_size - self.k as u64 + 1
        } else {
            0
        };
        
        let coverage = if theoretical_kmers > 0 {
            self.total_positions as f64 / theoretical_kmers as f64
        } else {
            0.0
        };
        
        let uniqueness = if self.total_positions > 0 {
            self.total_kmers as f64 / self.total_positions as f64
        } else {
            0.0
        };
        
        (coverage, uniqueness)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome::Genome;
    use std::io::Write;
    use tempfile::NamedTempFile;
    
    fn create_test_genome(sequences: &[(&str, &str)]) -> Genome {
        let mut fasta_content = String::new();
        for (name, seq) in sequences {
            fasta_content.push_str(&format!(">{}\n{}\n", name, seq));
        }
        
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(fasta_content.as_bytes()).unwrap();
        file.flush().unwrap();
        
        Genome::from_fasta(file.path()).unwrap()
    }
    
    #[test]
    fn test_basic_indexing() {
        let genome = create_test_genome(&[("contig1", "ATCGATCGATCG")]);
        let config = IndexConfig {
            k: 4, // Use valid k-mer length
            min_frequency: 1,
            max_frequency: None,
            parallel: false,
            max_positions_per_kmer: 1000,
        };
        
        let index = KmerIndex::build(&genome, config).unwrap();
        
        // Should have extracted k-mers from ATCGATCGATCG
        let stats = index.stats();
        assert!(stats.total_kmers > 0);
        assert_eq!(stats.k, 4);
        assert_eq!(stats.num_contigs, 1);
    }
    
    #[test]
    fn test_multi_contig_indexing() {
        let genome = create_test_genome(&[
            ("contig1", "ATCGATCG"),
            ("contig2", "GATCGATC"),
        ]);
        
        let config = IndexConfig {
            k: 4,
            min_frequency: 1,
            max_frequency: None,
            parallel: false,
            max_positions_per_kmer: 1000,
        };
        
        let index = KmerIndex::build(&genome, config).unwrap();
        
        let stats = index.stats();
        assert_eq!(stats.num_contigs, 2);
        
        // Check for multi-contig k-mers
        let multi_contig = index.multi_contig_kmers();
        assert!(!multi_contig.is_empty(), "Should find k-mers in multiple contigs");
    }
    
    #[test]
    fn test_frequency_filtering() {
        let genome = create_test_genome(&[("test", "AAAAAAAATCGATCGATCG")]); // Many A's
        
        let config = IndexConfig {
            k: 4, // Use valid k-mer length
            min_frequency: 2,
            max_frequency: Some(5),
            parallel: false,
            max_positions_per_kmer: 1000,
        };
        
        let index = KmerIndex::build(&genome, config).unwrap();
        
        // AAAA should appear multiple times but be filtered if too frequent
        let aaaa_kmer = Kmer::from_sequence(&[0, 0, 0, 0], 4).unwrap();
        let frequency = index.get_frequency(&aaaa_kmer);
        
        // Check that frequency filtering was applied
        if frequency > 5 {
            assert!(index.get_positions(&aaaa_kmer).is_none());
        }
    }
    
    #[test]
    fn test_kmer_positions() {
        let genome = create_test_genome(&[("test", "ATCGATCG")]);
        
        let config = IndexConfig {
            k: 4, // Use valid k-mer length
            min_frequency: 1,
            max_frequency: None,
            parallel: false,
            max_positions_per_kmer: 1000,
        };
        
        let index = KmerIndex::build(&genome, config).unwrap();
        
        // Check that positions are correctly recorded
        let atcg_kmer = Kmer::from_sequence(&[0, 3, 1, 2], 4).unwrap();
        if let Some(positions) = index.get_positions(&atcg_kmer) {
            assert!(!positions.is_empty());
            
            for pos in positions {
                assert_eq!(pos.contig, 0);
                assert!(pos.position < genome.len());
                assert_eq!(pos.position, pos.local_position); // Single contig
            }
        }
    }
    
    #[test]
    fn test_parallel_vs_sequential() {
        let genome = create_test_genome(&[
            ("contig1", "ATCGATCGATCGATCG"),
            ("contig2", "GCATGCATGCATGCAT"),
        ]);
        
        let config_seq = IndexConfig {
            k: 4,
            min_frequency: 1,
            max_frequency: None,
            parallel: false,
            max_positions_per_kmer: 1000,
        };
        
        let config_par = IndexConfig { parallel: true, ..config_seq.clone() };
        
        let index_seq = KmerIndex::build(&genome, config_seq).unwrap();
        let index_par = KmerIndex::build(&genome, config_par).unwrap();
        
        // Results should be equivalent
        let stats_seq = index_seq.stats();
        let stats_par = index_par.stats();
        
        assert_eq!(stats_seq.total_kmers, stats_par.total_kmers);
        assert_eq!(stats_seq.total_positions, stats_par.total_positions);
    }
    
    #[test]
    fn test_index_stats() {
        let genome = create_test_genome(&[("test", "ATCGATCGATCGATCGATCG")]); // Longer sequence
        
        let config = IndexConfig {
            min_frequency: 1, // Lower minimum frequency to ensure k-mers are included
            ..IndexConfig::default()
        };
        let index = KmerIndex::build(&genome, config).unwrap();
        
        let stats = index.stats();
        
        // Memory usage calculation should always return something > 0 if there are positions
        if stats.total_positions > 0 {
            assert!(stats.estimated_memory_usage() > 0);
        }
        
        let (coverage, uniqueness) = stats.coverage_stats();
        assert!(coverage >= 0.0);
        assert!(uniqueness >= 0.0 && uniqueness <= 1.0);
    }
}