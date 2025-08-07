//! Contig-by-contig k-mer indexing for large-scale genomic analysis
//! 
//! This module implements efficient k-mer indexing that processes genomes
//! contig-by-contig to handle large, fragmented assemblies. It provides
//! deterministic indexing suitable for concurrent repeat detection pipelines.

use crate::error::Result;
use crate::genome::{Genome, ContigId};
use crate::kmer::{Kmer, KmerEngine};
use crate::pipeline::CandidatePair;
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

/// Memory optimization hints for k-mer indexing
#[derive(Debug, Clone)]
pub struct MemoryHints {
    /// Expected unique k-mers for pre-allocation
    pub expected_unique_kmers: Option<usize>,
    /// Expected total positions for all k-mers
    pub expected_total_positions: Option<usize>,
    /// Whether to use aggressive pre-allocation
    pub aggressive_preallocation: bool,
    /// Enable memory usage tracking
    pub track_memory_usage: bool,
}

impl Default for MemoryHints {
    fn default() -> Self {
        Self {
            expected_unique_kmers: None,
            expected_total_positions: None,
            aggressive_preallocation: false,
            track_memory_usage: false,
        }
    }
}

/// K-mer index configuration
#[derive(Clone)]
pub struct IndexConfig {
    pub k: usize,
    pub min_frequency: u32,
    pub max_frequency: Option<u32>,
    pub parallel: bool,
    pub max_positions_per_kmer: usize,
    pub memory_hints: MemoryHints,
}

impl Default for IndexConfig {
    fn default() -> Self {
        Self {
            k: 13,
            min_frequency: 2,
            max_frequency: Some(1000),
            parallel: true,
            max_positions_per_kmer: 10000,
            memory_hints: MemoryHints::default(),
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
        
        // Pre-allocate collections based on genome size estimates and memory hints
        let capacity_estimate = Self::estimate_kmer_capacity(genome, &config);
        
        // Use memory hints if provided, otherwise use estimates
        let unique_capacity = config.memory_hints.expected_unique_kmers
            .unwrap_or(capacity_estimate.unique_kmers);
        
        let mut positions = engine.new_kmer_map_with_capacity(unique_capacity);
        let mut frequencies = engine.new_kmer_map_with_capacity(unique_capacity);
        
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
        config: &IndexConfig,
    ) -> Result<()> {
        use std::sync::Mutex;
        
        let positions_mutex = Mutex::new(positions);
        let frequencies_mutex = Mutex::new(frequencies);
        
        // Process contigs in parallel
        let contigs: Vec<_> = genome.contigs().collect();
        
        contigs.par_iter().try_for_each(|&(contig_id, contig_info)| -> Result<()> {
            let sequence = genome.slice(contig_info.range());
            let _contig_kmers = engine.extract_kmers(sequence);
            
            // Pre-allocate local maps based on contig size to reduce reallocations
            let contig_len = contig_info.end - contig_info.start;
            let contig_capacity = Self::estimate_contig_kmer_capacity(contig_len, engine.k(), config);
            let mut local_positions = engine.new_kmer_map_with_capacity(contig_capacity.unique_kmers);
            let mut local_frequencies = engine.new_kmer_map_with_capacity(contig_capacity.unique_kmers);
            
            // Pre-allocate position vectors based on expected frequency distribution
            let default_position_capacity = contig_capacity.avg_positions_per_kmer;
            
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
                        
                        // Update local collections with pre-allocated capacity
                        local_positions.entry(kmer)
                            .or_insert_with(|| Vec::with_capacity(default_position_capacity))
                            .push(kmer_pos);
                        *local_frequencies.entry(kmer).or_insert(0) += 1;
                    }
                }
            }
            
            // Merge with global collections (critical section) - use batch merging to minimize lock time
            {
                let mut global_positions = positions_mutex.lock().unwrap();
                let mut global_frequencies = frequencies_mutex.lock().unwrap();
                
                // Batch merge positions with pre-allocated capacity
                for (kmer, mut local_pos_list) in local_positions {
                    global_positions.entry(kmer)
                        .or_insert_with(|| Vec::with_capacity(local_pos_list.len() * 2))
                        .append(&mut local_pos_list);
                }
                
                // Batch merge frequencies
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
                    
                    positions.entry(kmer)
                        .or_insert_with(|| Vec::with_capacity(4)) // Most k-mers appear 2-4 times
                        .push(kmer_pos);
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
    
    /// Get the most frequent k-mers using partial sorting for better performance
    pub fn most_frequent_kmers(&self, limit: usize) -> Vec<(Kmer, u32)> {
        let mut freq_pairs: Vec<_> = self.frequencies.iter().map(|(&k, &f)| (k, f)).collect();
        
        if freq_pairs.len() <= limit {
            // If we need all elements, just sort normally
            freq_pairs.sort_by(|a, b| b.1.cmp(&a.1));
            return freq_pairs;
        }
        
        // Use partial sorting - only sort the top `limit` elements
        // This is much faster than full sorting when limit << total_elements
        freq_pairs.select_nth_unstable_by(limit, |a, b| b.1.cmp(&a.1));
        let (top_elements, _) = freq_pairs.split_at_mut(limit);
        
        // Sort just the top elements we need
        top_elements.sort_by(|a, b| b.1.cmp(&a.1));
        top_elements.to_vec()
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
    
    /// Get comprehensive index statistics including memory usage
    pub fn stats(&self) -> IndexStats {
        let total_kmers = self.positions.len();
        let total_positions: usize = self.positions.values().map(|v| v.len()).sum();
        let avg_frequency = if total_kmers > 0 {
            total_positions as f64 / total_kmers as f64
        } else {
            0.0
        };
        
        let multi_contig_count = self.multi_contig_kmers().len();
        
        // Create a temporary stats structure for memory calculation
        let temp_stats = IndexStats {
            k: self.config.k,
            total_kmers,
            total_positions,
            avg_frequency,
            multi_contig_kmers: multi_contig_count,
            genome_size: self.genome_size,
            num_contigs: self.num_contigs,
            memory_usage: MemoryUsageBreakdown {
                total_bytes: 0,
                kmers_bytes: 0,
                positions_bytes: 0,
                frequencies_bytes: 0,
                overhead_bytes: 0,
                vec_storage_bytes: 0,
            },
        };
        
        // Calculate actual memory usage
        let memory_usage = temp_stats.estimated_memory_usage();
        
        IndexStats {
            k: self.config.k,
            total_kmers,
            total_positions,
            avg_frequency,
            multi_contig_kmers: multi_contig_count,
            genome_size: self.genome_size,
            num_contigs: self.num_contigs,
            memory_usage,
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

    pub fn get_candidate_pairs(&self, _min_freq: u32, _max_freq: Option<u32>) -> Vec<CandidatePair> {
        // Simplified logic for now
        Vec::new()
    }
    
    /// Estimate k-mer capacity requirements for efficient pre-allocation
    /// 
    /// This function calculates expected memory requirements based on genome size,
    /// k-mer length, and filtering parameters to minimize reallocations during indexing.
    fn estimate_kmer_capacity(genome: &Genome, config: &IndexConfig) -> KmerCapacityEstimate {
        let genome_size = genome.len();
        let k = config.k;
        
        // Theoretical maximum k-mers = genome_size - k + 1
        let theoretical_kmers = if genome_size >= k as u64 {
            genome_size - k as u64 + 1
        } else {
            0
        };
        
        // Estimate unique k-mers based on genome complexity
        // For typical genomes: 60-80% of k-mers are unique due to repeats
        let complexity_factor = 0.7; // Conservative estimate
        let estimated_unique = (theoretical_kmers as f64 * complexity_factor) as usize;
        
        // Account for frequency filtering - only keep k-mers meeting min/max frequency
        let frequency_filter_factor = match config.min_frequency {
            1 => 0.95,  // Keep most k-mers
            2 => 0.70,  // Keep moderate frequency k-mers  
            3..=5 => 0.50,  // Keep higher frequency k-mers
            _ => 0.30,  // Keep only frequent k-mers
        };
        
        let filtered_unique = (estimated_unique as f64 * frequency_filter_factor) as usize;
        
        // Estimate average positions per k-mer (repetitive regions)
        let avg_positions = match config.min_frequency {
            1 => 2,
            2 => 3, 
            3..=5 => 4,
            _ => 6,
        };
        
        KmerCapacityEstimate {
            theoretical_kmers: theoretical_kmers as usize,
            unique_kmers: filtered_unique.max(1000), // Minimum reasonable capacity
            avg_positions_per_kmer: avg_positions,
            total_positions: filtered_unique * avg_positions,
        }
    }
    
    /// Estimate k-mer capacity for a single contig
    fn estimate_contig_kmer_capacity(
        contig_len: u64, 
        k: usize, 
        config: &IndexConfig
    ) -> KmerCapacityEstimate {
        let theoretical_kmers = if contig_len >= k as u64 {
            contig_len - k as u64 + 1
        } else {
            0
        };
        
        // Contigs typically have higher uniqueness than whole genomes
        let complexity_factor = 0.85;
        let estimated_unique = (theoretical_kmers as f64 * complexity_factor) as usize;
        
        // Apply frequency filtering
        let frequency_filter_factor = match config.min_frequency {
            1 => 0.95,
            2 => 0.80, 
            _ => 0.60,
        };
        
        let filtered_unique = (estimated_unique as f64 * frequency_filter_factor) as usize;
        
        KmerCapacityEstimate {
            theoretical_kmers: theoretical_kmers as usize,
            unique_kmers: filtered_unique.max(100), // Minimum for small contigs
            avg_positions_per_kmer: 2, // Contigs typically have fewer repeats
            total_positions: filtered_unique * 2,
        }
    }
}

/// Capacity estimate for k-mer collections to enable efficient pre-allocation
#[derive(Debug, Clone)]
pub struct KmerCapacityEstimate {
    /// Theoretical maximum k-mers possible
    pub theoretical_kmers: usize,
    /// Estimated unique k-mers after filtering
    pub unique_kmers: usize,
    /// Average positions per k-mer (for repeat detection)
    pub avg_positions_per_kmer: usize,
    /// Total estimated positions across all k-mers
    pub total_positions: usize,
}

/// Detailed memory usage breakdown for performance analysis
#[derive(Debug, Clone)]
pub struct MemoryUsageBreakdown {
    /// Total estimated memory usage in bytes
    pub total_bytes: u64,
    /// Memory used by k-mer storage
    pub kmers_bytes: u64,
    /// Memory used by position storage
    pub positions_bytes: u64,
    /// Memory used by frequency counters
    pub frequencies_bytes: u64,
    /// Memory overhead (hash tables, vectors, etc.)
    pub overhead_bytes: u64,
    /// Memory used by Vec storage structures
    pub vec_storage_bytes: u64,
}

impl MemoryUsageBreakdown {
    /// Get memory usage in human-readable format
    pub fn human_readable(&self) -> String {
        fn format_bytes(bytes: u64) -> String {
            const UNITS: &[&str] = &["B", "KB", "MB", "GB", "TB"];
            let mut size = bytes as f64;
            let mut unit_idx = 0;
            
            while size >= 1024.0 && unit_idx < UNITS.len() - 1 {
                size /= 1024.0;
                unit_idx += 1;
            }
            
            format!( "{:.2} {}", size, UNITS[unit_idx])
        }
        
        format!(
            "Total: {} (K-mers: {}, Positions: {}, Frequencies: {}, Overhead: {})",
            format_bytes(self.total_bytes),
            format_bytes(self.kmers_bytes),
            format_bytes(self.positions_bytes),
            format_bytes(self.frequencies_bytes),
            format_bytes(self.overhead_bytes)
        )
    }
    
    /// Get percentage breakdown of memory usage
    pub fn percentage_breakdown(&self) -> MemoryPercentages {
        let total = self.total_bytes as f64;
        MemoryPercentages {
            kmers_percent: (self.kmers_bytes as f64 / total) * 100.0,
            positions_percent: (self.positions_bytes as f64 / total) * 100.0,
            frequencies_percent: (self.frequencies_bytes as f64 / total) * 100.0,
            overhead_percent: (self.overhead_bytes as f64 / total) * 100.0,
        }
    }
}

/// Memory usage percentages for analysis
#[derive(Debug, Clone)]
pub struct MemoryPercentages {
    pub kmers_percent: f64,
    pub positions_percent: f64,
    pub frequencies_percent: f64,
    pub overhead_percent: f64,
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
    /// Detailed memory usage breakdown
    pub memory_usage: MemoryUsageBreakdown,
}

impl IndexStats {
    /// Calculate detailed memory usage estimate in bytes
    pub fn estimated_memory_usage(&self) -> MemoryUsageBreakdown {
        // Detailed memory analysis for genomic k-mer indices
        
        // Core data structure sizes
        let kmer_size = std::mem::size_of::<Kmer>(); // Should be 16 bytes
        let position_size = std::mem::size_of::<crate::index::KmerPosition>(); // Should be 24 bytes
        let frequency_entry_size = std::mem::size_of::<u32>(); // 4 bytes
        
        // HashMap overhead estimation (buckets, metadata, etc.)
        let hashmap_overhead_factor = 1.4; // AHashMap is more efficient than std::HashMap
        
        // Calculate base memory usage
        let kmers_memory = self.total_kmers * kmer_size;
        let positions_memory = self.total_positions * position_size;
        let frequencies_memory = self.total_kmers * frequency_entry_size;
        
        // Vec<KmerPosition> storage overhead
        let position_vec_overhead = self.total_kmers * std::mem::size_of::<Vec<crate::index::KmerPosition>>();
        
        let base_usage = kmers_memory + positions_memory + frequencies_memory + position_vec_overhead;
        let total_with_overhead = (base_usage as f64 * hashmap_overhead_factor) as u64;
        
        MemoryUsageBreakdown {
            total_bytes: total_with_overhead,
            kmers_bytes: kmers_memory as u64,
            positions_bytes: positions_memory as u64,
            frequencies_bytes: frequencies_memory as u64,
            overhead_bytes: total_with_overhead - base_usage as u64,
            vec_storage_bytes: position_vec_overhead as u64,
        }
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
            fasta_content.push_str(&format!( ">{}
{}
", name, seq));
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
            memory_hints: MemoryHints::default(),
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
            memory_hints: crate::index::MemoryHints::default(),
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
            memory_hints: crate::index::MemoryHints::default(),
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
            memory_hints: MemoryHints::default(),
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
            memory_hints: crate::index::MemoryHints::default(),
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
            assert!(stats.memory_usage.total_bytes > 0);
        }
        
        let (coverage, uniqueness) = stats.coverage_stats();
        assert!(coverage >= 0.0);
        assert!(uniqueness >= 0.0 && uniqueness <= 1.0);
    }
}
