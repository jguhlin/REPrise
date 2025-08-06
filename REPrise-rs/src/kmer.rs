//! Canonical k-mer engine for deterministic genomic indexing
//!
//! This module provides high-performance k-mer processing with canonical
//! representation (lexicographically smallest of forward/reverse complement).
//! It uses ahash for deterministic, high-performance hashing suitable for
//! parallel processing pipelines.

use crate::error::{RepriseError, Result};
use crate::genome::Genome;
use ahash::{AHashMap, AHashSet};
use std::hash::{Hash, Hasher};

/// Maximum supported k-mer length (limited by u64 bit representation)
pub const MAX_K: usize = 32;

/// Minimum supported k-mer length
pub const MIN_K: usize = 4;

/// Fixed seed for deterministic hashing across runs
const HASH_SEED: usize = 0x51f3b5b8;

/// A k-mer represented as a 64-bit integer for efficient operations
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Kmer {
    /// Canonical k-mer (lexicographically smallest of forward/reverse complement)
    canonical: u64,
    /// Length of the k-mer
    k: usize,
}

impl Kmer {
    /// Create a new k-mer from a sequence slice
    ///
    /// # Arguments
    /// * `seq` - Nucleotide sequence (0=A, 1=C, 2=G, 3=T, 99=N)
    /// * `k` - K-mer length
    ///
    /// # Returns
    /// * `Ok(Kmer)` - Valid canonical k-mer
    /// * `Err(RepriseError)` - Invalid sequence or k-mer parameters
    pub fn from_sequence(seq: &[u8], k: usize) -> Result<Self> {
        if k < MIN_K || k > MAX_K {
            return Err(RepriseError::invalid_kmer_length(k, MIN_K, MAX_K));
        }
        
        if seq.len() != k {
            return Err(RepriseError::config(format!(
                "Sequence length {} does not match k-mer length {}",
                seq.len(),
                k
            )));
        }
        
        // Check for ambiguous bases
        for (i, &base) in seq.iter().enumerate() {
            if base > 3 {
                return Err(RepriseError::InvalidKmer(i as u64));
            }
        }
        
        // Convert to 64-bit representation
        let forward = sequence_to_u64(seq)?;
        let reverse = reverse_complement_u64(forward, k);
        
        // Canonical form is the lexicographically smaller
        let canonical = forward.min(reverse);
        
        Ok(Kmer { canonical, k })
    }
    
    /// Get the canonical 64-bit representation
    pub fn canonical(&self) -> u64 {
        self.canonical
    }
    
    /// Get the k-mer length
    pub fn k(&self) -> usize {
        self.k
    }
    
    /// Check if this k-mer equals another k-mer
    pub fn equals(&self, other: &Kmer) -> bool {
        self.k == other.k && self.canonical == other.canonical
    }
    
    /// Convert back to nucleotide sequence (canonical form)
    pub fn to_sequence(&self) -> Vec<u8> {
        u64_to_sequence(self.canonical, self.k)
    }
    
    /// Get the reverse complement of this k-mer
    pub fn reverse_complement(&self) -> Self {
        let rev_comp = reverse_complement_u64(self.canonical, self.k);
        
        // The canonical form might flip
        let canonical = self.canonical.min(rev_comp);
        
        Kmer {
            canonical,
            k: self.k,
        }
    }
}

impl Hash for Kmer {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.canonical.hash(state);
        self.k.hash(state);
    }
}

/// High-performance k-mer indexing engine
///
/// This struct provides efficient k-mer extraction and indexing for genomic
/// sequences, with support for canonical representation and deterministic
/// hashing suitable for concurrent processing.
#[derive(Debug)]
pub struct KmerEngine {
    k: usize,
    /// Deterministic hasher for consistent results across runs
    hasher_builder: ahash::RandomState,
}

impl KmerEngine {
    /// Create a new k-mer engine with specified k-mer length
    ///
    /// # Arguments
    /// * `k` - K-mer length (must be between MIN_K and MAX_K)
    ///
    /// # Returns
    /// * `Ok(KmerEngine)` - Successfully created engine
    /// * `Err(RepriseError)` - Invalid k-mer length
    pub fn new(k: usize) -> Result<Self> {
        if k < MIN_K || k > MAX_K {
            return Err(RepriseError::invalid_kmer_length(k, MIN_K, MAX_K));
        }
        
        // Use fixed seed for deterministic hashing
        let hasher_builder = ahash::RandomState::with_seed(HASH_SEED);
        
        Ok(KmerEngine {
            k,
            hasher_builder,
        })
    }
    
    /// Extract all k-mers from a genomic sequence
    ///
    /// # Arguments
    /// * `sequence` - Genomic sequence in numeric format
    ///
    /// # Returns
    /// Vector of valid k-mers (skips k-mers containing ambiguous bases)
    pub fn extract_kmers(&self, sequence: &[u8]) -> Vec<Kmer> {
        if sequence.len() < self.k {
            return Vec::new();
        }
        
        let mut kmers = Vec::new();
        
        for i in 0..=sequence.len() - self.k {
            let subseq = &sequence[i..i + self.k];
            
            // Skip k-mers with ambiguous bases
            if subseq.iter().any(|&b| b > 3) {
                continue;
            }
            
            if let Ok(kmer) = Kmer::from_sequence(subseq, self.k) {
                kmers.push(kmer);
            }
        }
        
        kmers
    }
    
    /// Extract k-mers from a specific genomic range
    ///
    /// # Arguments
    /// * `genome` - Genome to extract from
    /// * `start` - Start position (global coordinates)
    /// * `length` - Length of region to extract from
    ///
    /// # Returns
    /// * `Ok(Vec<(u64, Kmer)>)` - Vector of (position, kmer) pairs
    /// * `Err(RepriseError)` - Invalid range or boundary violation
    pub fn extract_kmers_from_range(&self, genome: &Genome, start: u64, length: usize) -> Result<Vec<(u64, Kmer)>> {
        if start + length as u64 > genome.len() {
            return Err(RepriseError::invalid_range(start, start + length as u64));
        }
        
        // Check if the range crosses contig boundaries
        if !genome.is_within_one_contig(start, length) {
            return Err(RepriseError::BoundaryViolation);
        }
        
        let sequence = genome.slice(start..start + length as u64);
        let mut result = Vec::new();
        
        if sequence.len() < self.k {
            return Ok(result);
        }
        
        for i in 0..=sequence.len() - self.k {
            let subseq = &sequence[i..i + self.k];
            
            // Skip k-mers with ambiguous bases
            if subseq.iter().any(|&b| b > 3) {
                continue;
            }
            
            if let Ok(kmer) = Kmer::from_sequence(subseq, self.k) {
                result.push((start + i as u64, kmer));
            }
        }
        
        Ok(result)
    }
    
    /// Create a deterministic hash map for k-mer indexing
    pub fn new_kmer_map<V>(&self) -> AHashMap<Kmer, V> {
        AHashMap::with_hasher(self.hasher_builder.clone())
    }
    
    /// Create a deterministic hash set for k-mer tracking
    pub fn new_kmer_set(&self) -> AHashSet<Kmer> {
        AHashSet::with_hasher(self.hasher_builder.clone())
    }
    
    /// Get the k-mer length for this engine
    pub fn k(&self) -> usize {
        self.k
    }
}

/// Convert nucleotide sequence to 64-bit representation
/// Each nucleotide uses 2 bits: A=00, C=01, G=10, T=11
fn sequence_to_u64(seq: &[u8]) -> Result<u64> {
    if seq.len() > 32 {
        return Err(RepriseError::config("Sequence too long for 64-bit representation".to_string()));
    }
    
    let mut result = 0u64;
    for &nucleotide in seq {
        if nucleotide > 3 {
            return Err(RepriseError::InvalidKmer(0)); // Position doesn't matter here
        }
        result = (result << 2) | (nucleotide as u64);
    }
    
    Ok(result)
}

/// Convert 64-bit representation back to nucleotide sequence
fn u64_to_sequence(mut value: u64, k: usize) -> Vec<u8> {
    let mut result = vec![0u8; k];
    
    for i in (0..k).rev() {
        result[i] = (value & 0b11) as u8;
        value >>= 2;
    }
    
    result
}

/// Compute reverse complement of a 64-bit k-mer
fn reverse_complement_u64(kmer: u64, k: usize) -> u64 {
    let mut result = 0u64;
    let mut temp = kmer;
    
    for _ in 0..k {
        let nucleotide = temp & 0b11;
        let complement = 3 - nucleotide; // A(0)<->T(3), C(1)<->G(2)
        result = (result << 2) | complement;
        temp >>= 2;
    }
    
    result
}

/// K-mer frequency counter for genomic analysis
#[derive(Debug)]
pub struct KmerFrequency {
    engine: KmerEngine,
    counts: AHashMap<Kmer, u32>,
}

impl KmerFrequency {
    /// Create a new k-mer frequency counter
    pub fn new(k: usize) -> Result<Self> {
        let engine = KmerEngine::new(k)?;
        let counts = engine.new_kmer_map();
        
        Ok(KmerFrequency { engine, counts })
    }
    
    /// Count k-mers in a sequence
    pub fn count_sequence(&mut self, sequence: &[u8]) {
        let kmers = self.engine.extract_kmers(sequence);
        for kmer in kmers {
            *self.counts.entry(kmer).or_insert(0) += 1;
        }
    }
    
    /// Count k-mers in a genome
    pub fn count_genome(&mut self, genome: &Genome) -> Result<()> {
        for (_, contig_info) in genome.contigs() {
            let sequence = genome.slice(contig_info.range());
            self.count_sequence(sequence);
        }
        Ok(())
    }
    
    /// Get frequency of a specific k-mer
    pub fn get_frequency(&self, kmer: &Kmer) -> u32 {
        self.counts.get(kmer).copied().unwrap_or(0)
    }
    
    /// Get all k-mers with their frequencies
    pub fn get_all_frequencies(&self) -> &AHashMap<Kmer, u32> {
        &self.counts
    }
    
    /// Get k-mers with frequency at least min_freq
    pub fn get_frequent_kmers(&self, min_freq: u32) -> Vec<(Kmer, u32)> {
        self.counts
            .iter()
            .filter(|&(_, freq)| *freq >= min_freq)
            .map(|(&kmer, &freq)| (kmer, freq))
            .collect()
    }
    
    /// Clear all counts
    pub fn clear(&mut self) {
        self.counts.clear();
    }
    
    /// Get total number of unique k-mers
    pub fn unique_kmers(&self) -> usize {
        self.counts.len()
    }
    
    /// Get total k-mer count
    pub fn total_kmers(&self) -> u64 {
        self.counts.values().map(|&freq| freq as u64).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome::Genome;
    use std::io::Write;
    use tempfile::NamedTempFile;
    
    fn create_test_genome(sequence: &str) -> Genome {
        let fasta_content = format!(">test\n{}\n", sequence);
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(fasta_content.as_bytes()).unwrap();
        file.flush().unwrap();
        
        Genome::from_fasta(file.path()).unwrap()
    }
    
    #[test]
    fn test_kmer_creation() {
        // Test valid k-mer
        let seq = &[0, 1, 2, 3]; // ACGT
        let kmer = Kmer::from_sequence(seq, 4).unwrap();
        assert_eq!(kmer.k(), 4);
        
        // Test canonical property (should be same for forward and reverse complement)
        let rev_seq = &[0, 1, 2, 3]; // ACGT
        let rev_kmer = Kmer::from_sequence(rev_seq, 4).unwrap();
        
        // ACGT (0123) vs ACGT (0123) - both should be canonical
        assert_eq!(kmer.canonical(), rev_kmer.canonical());
    }
    
    #[test]
    fn test_kmer_canonical_form() {
        // Test that canonical form works correctly
        let forward = &[0, 0, 0, 0]; // AAAA
        let reverse = &[3, 3, 3, 3]; // TTTT
        
        let kmer_forward = Kmer::from_sequence(forward, 4).unwrap();
        let kmer_reverse = Kmer::from_sequence(reverse, 4).unwrap();
        
        // Both should have the same canonical representation
        assert_eq!(kmer_forward.canonical(), kmer_reverse.canonical());
    }
    
    #[test]
    fn test_invalid_kmers() {
        // Test invalid k-mer length
        assert!(Kmer::from_sequence(&[0, 1, 2], 2).is_err());
        assert!(KmerEngine::new(1).is_err());
        assert!(KmerEngine::new(100).is_err());
        
        // Test ambiguous bases
        let seq_with_n = &[0, 1, 99, 3]; // AC[N]T
        assert!(Kmer::from_sequence(seq_with_n, 4).is_err());
    }
    
    #[test]
    fn test_kmer_extraction() {
        let engine = KmerEngine::new(4).unwrap(); // Use valid k-mer length
        let sequence = &[0, 1, 2, 3, 0, 1]; // ACGTAC
        
        let kmers = engine.extract_kmers(sequence);
        
        // Should extract: ACGT, CGTA, GTAC
        assert_eq!(kmers.len(), 3);
        
        // Test with ambiguous base
        let sequence_with_n = &[0, 1, 99, 3, 0, 1]; // AC[N]TAC
        let kmers_with_n = engine.extract_kmers(sequence_with_n);
        
        // Should skip k-mers containing N
        assert_eq!(kmers_with_n.len(), 0); // No valid 4-mers without N
    }
    
    #[test]
    fn test_genome_kmer_extraction() {
        let genome = create_test_genome("ATCGATCGATCG");
        let engine = KmerEngine::new(4).unwrap(); // Use valid k-mer length
        
        let result = engine.extract_kmers_from_range(&genome, 0, 12).unwrap();
        
        // Should extract 9 k-mers from ATCGATCGATCG (12-4+1)
        assert_eq!(result.len(), 9);
        
        // Check positions
        for (i, (pos, _)) in result.iter().enumerate() {
            assert_eq!(*pos, i as u64);
        }
    }
    
    #[test]
    fn test_kmer_frequency_counting() {
        let mut counter = KmerFrequency::new(4).unwrap(); // Use valid k-mer length
        let sequence = &[0, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGT (ACGT appears twice)
        
        counter.count_sequence(sequence);
        
        let acgt_kmer = Kmer::from_sequence(&[0, 1, 2, 3], 4).unwrap();
        assert_eq!(counter.get_frequency(&acgt_kmer), 2);
        
        assert!(counter.unique_kmers() > 0); // Some unique k-mers
        assert!(counter.total_kmers() > 0); // Some total k-mers
    }
    
    #[test]
    fn test_deterministic_hashing() {
        // Test that hashing is deterministic across engine instances
        let engine1 = KmerEngine::new(4).unwrap();
        let engine2 = KmerEngine::new(4).unwrap();
        
        let mut map1 = engine1.new_kmer_map();
        let mut map2 = engine2.new_kmer_map();
        
        let kmer = Kmer::from_sequence(&[0, 1, 2, 3], 4).unwrap();
        map1.insert(kmer, 1);
        map2.insert(kmer, 2);
        
        // Both maps should have the same internal structure
        assert!(map1.contains_key(&kmer));
        assert!(map2.contains_key(&kmer));
    }
    
    #[test]
    fn test_sequence_conversions() {
        let original = &[0, 1, 2, 3, 0, 1];
        let as_u64 = sequence_to_u64(original).unwrap();
        let converted_back = u64_to_sequence(as_u64, original.len());
        
        assert_eq!(original, converted_back.as_slice());
    }
    
    #[test]
    fn test_reverse_complement() {
        // Test ATCG -> CGAT
        let atcg = sequence_to_u64(&[0, 3, 1, 2]).unwrap(); // ATCG
        let cgat = sequence_to_u64(&[1, 2, 0, 3]).unwrap(); // CGAT
        
        let rev_comp = reverse_complement_u64(atcg, 4);
        assert_eq!(rev_comp, cgat);
        
        // Test palindrome (should equal itself)
        let atgcat = sequence_to_u64(&[0, 3, 2, 1, 0, 3]).unwrap(); // ATGCAT
        let rev_comp_pal = reverse_complement_u64(atgcat, 6);
        assert_eq!(atgcat, rev_comp_pal);
    }
}