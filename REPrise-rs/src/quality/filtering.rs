//! Quality filtering for k-mer candidates
//! 
//! Implements entropy filtering, frequency stratification, and basic quality checks
//! to eliminate low-quality candidates before expensive extension alignment.

use crate::quality::QualityError;
use std::collections::HashMap;

/// Quality filter for k-mer pre-screening
pub struct QualityFilter {
    pub max_entropy: f64,         // C++ MAXENTROPY = -0.7
    pub min_length: usize,        // C++ MINLENGTH = 50
    pub min_freq: usize,          // C++ MINFREQ = 3
    pub gc_content_range: (f64, f64), // Biological GC content limits
    pub max_homopolymer: usize,   // Maximum homopolymer run length
    pub max_simple_repeat: f64,   // Maximum simple repeat content
}

impl Default for QualityFilter {
    fn default() -> Self {
        Self {
            max_entropy: -0.7,           // Match C++ threshold
            min_length: 50,              // Match C++ minimum
            min_freq: 3,                 // Match C++ minimum
            gc_content_range: (0.2, 0.8), // 20-80% GC content
            max_homopolymer: 8,          // Max homopolymer run
            max_simple_repeat: 0.8,      // Max 80% simple repeats
        }
    }
}

impl QualityFilter {
    pub fn new(max_entropy: f64, min_length: usize, min_freq: usize) -> Self {
        Self {
            max_entropy,
            min_length,
            min_freq,
            ..Default::default()
        }
    }

    /// Pre-filter k-mer before expensive extension alignment
    pub fn prefilter_kmer(&self, kmer: &[u8], frequency: usize) -> bool {
        // Basic frequency check
        if frequency < self.min_freq {
            return false;
        }

        // Entropy filtering (matches C++ compute_entropy)
        if self.compute_entropy(kmer) <= self.max_entropy {
            return false;
        }

        // GC content check
        if !self.check_gc_content(kmer) {
            return false;
        }

        // Homopolymer run check
        if self.has_excessive_homopolymer(kmer) {
            return false;
        }

        // Simple repeat check
        if self.has_excessive_simple_repeats(kmer) {
            return false;
        }

        true
    }

    /// Compute entropy of k-mer sequence (matches C++ implementation)
    pub fn compute_entropy(&self, sequence: &[u8]) -> f64 {
        if sequence.is_empty() {
            return 0.0;
        }

        let mut base_counts = [0; 4]; // A, C, G, T
        let mut valid_bases = 0;

        // Count bases
        for &base in sequence {
            if base < 4 {
                base_counts[base as usize] += 1;
                valid_bases += 1;
            }
        }

        if valid_bases == 0 {
            return 0.0;
        }

        // Calculate Shannon entropy
        let mut entropy = 0.0;
        for &count in &base_counts {
            if count > 0 {
                let probability = count as f64 / valid_bases as f64;
                entropy -= probability * probability.log2();
            }
        }

        entropy
    }

    /// Check GC content is within biological range
    fn check_gc_content(&self, sequence: &[u8]) -> bool {
        let mut gc_count = 0;
        let mut total_count = 0;

        for &base in sequence {
            if base < 4 {
                if base == 1 || base == 2 { // C or G
                    gc_count += 1;
                }
                total_count += 1;
            }
        }

        if total_count == 0 {
            return false;
        }

        let gc_content = gc_count as f64 / total_count as f64;
        gc_content >= self.gc_content_range.0 && gc_content <= self.gc_content_range.1
    }

    /// Check for excessive homopolymer runs
    fn has_excessive_homopolymer(&self, sequence: &[u8]) -> bool {
        if sequence.len() < 2 {
            return false;
        }

        let mut current_base = sequence[0];
        let mut current_run = 1;

        for &base in &sequence[1..] {
            if base == current_base {
                current_run += 1;
                if current_run > self.max_homopolymer {
                    return true;
                }
            } else {
                current_base = base;
                current_run = 1;
            }
        }

        false
    }

    /// Check for excessive simple repeat content
    fn has_excessive_simple_repeats(&self, sequence: &[u8]) -> bool {
        if sequence.len() < 4 {
            return false;
        }

        // Check for dinucleotide repeats (most common simple repeats)
        for period in 2..=4 {
            let repeat_content = self.calculate_repeat_content(sequence, period);
            if repeat_content > self.max_simple_repeat {
                return true;
            }
        }

        false
    }

    /// Calculate content of repeats with given period
    fn calculate_repeat_content(&self, sequence: &[u8], period: usize) -> f64 {
        if sequence.len() < period * 2 {
            return 0.0;
        }

        let mut matching_bases = 0;
        let mut total_comparable = 0;

        for i in period..sequence.len() {
            let current_base = sequence[i];
            let period_base = sequence[i - period];
            
            if current_base < 4 && period_base < 4 {
                if current_base == period_base {
                    matching_bases += 1;
                }
                total_comparable += 1;
            }
        }

        if total_comparable == 0 {
            return 0.0;
        }

        matching_bases as f64 / total_comparable as f64
    }

    /// Advanced filtering for candidate stratification
    pub fn stratify_candidates(&self, candidates: &[(Vec<u8>, usize)]) -> Vec<Vec<(Vec<u8>, usize)>> {
        let mut strata = vec![Vec::new(); 5]; // 5 quality strata

        for (kmer, freq) in candidates {
            if !self.prefilter_kmer(kmer, *freq) {
                continue; // Skip low-quality k-mers
            }

            // Classify into strata based on frequency and quality
            let quality_score = self.calculate_quality_score(kmer, *freq);
            let stratum = match quality_score {
                q if q >= 0.8 => 0, // Highest quality
                q if q >= 0.6 => 1, // High quality
                q if q >= 0.4 => 2, // Medium quality
                q if q >= 0.2 => 3, // Low quality
                _ => 4,             // Lowest quality
            };

            strata[stratum].push((kmer.clone(), *freq));
        }

        // Sort each stratum by frequency (descending)
        for stratum in &mut strata {
            stratum.sort_by(|a, b| b.1.cmp(&a.1));
        }

        strata
    }

    /// Calculate composite quality score for k-mer
    fn calculate_quality_score(&self, kmer: &[u8], frequency: usize) -> f64 {
        if kmer.is_empty() {
            return 0.0;
        }

        // Entropy component (normalized to 0-1)
        let entropy = self.compute_entropy(kmer);
        let entropy_score = (entropy / 2.0).min(1.0).max(0.0); // Max entropy ~2 for 4-base alphabet

        // Frequency component (log-scaled)
        let freq_score = if frequency <= 1 {
            0.0
        } else {
            (frequency as f64).ln() / 10.0
        }.min(1.0);

        // GC content component
        let gc_content = self.calculate_gc_content(kmer);
        let gc_score = 1.0 - ((gc_content - 0.5).abs() * 2.0); // Penalty for extreme GC

        // Length component
        let length_score = (kmer.len() as f64 / 13.0).min(1.0); // Normalized to typical k-mer length

        // Composite score (weighted average)
        0.3 * entropy_score + 0.3 * freq_score + 0.2 * gc_score + 0.2 * length_score
    }

    /// Calculate GC content as fraction
    fn calculate_gc_content(&self, sequence: &[u8]) -> f64 {
        let mut gc_count = 0;
        let mut total_count = 0;

        for &base in sequence {
            if base < 4 {
                if base == 1 || base == 2 { // C or G
                    gc_count += 1;
                }
                total_count += 1;
            }
        }

        if total_count == 0 {
            return 0.5; // Neutral GC content for empty/invalid sequence
        }

        gc_count as f64 / total_count as f64
    }

    /// Get filtering statistics
    pub fn get_filter_stats(&self, candidates: &[(Vec<u8>, usize)]) -> FilterStatistics {
        let total_candidates = candidates.len();
        let mut passed_entropy = 0;
        let mut passed_gc = 0;
        let mut passed_homopolymer = 0;
        let mut passed_simple_repeat = 0;
        let mut passed_frequency = 0;
        let mut passed_all = 0;

        for (kmer, freq) in candidates {
            if *freq >= self.min_freq {
                passed_frequency += 1;
            }
            
            if self.compute_entropy(kmer) > self.max_entropy {
                passed_entropy += 1;
            }
            
            if self.check_gc_content(kmer) {
                passed_gc += 1;
            }
            
            if !self.has_excessive_homopolymer(kmer) {
                passed_homopolymer += 1;
            }
            
            if !self.has_excessive_simple_repeats(kmer) {
                passed_simple_repeat += 1;
            }
            
            if self.prefilter_kmer(kmer, *freq) {
                passed_all += 1;
            }
        }

        FilterStatistics {
            total_candidates,
            passed_frequency,
            passed_entropy,
            passed_gc,
            passed_homopolymer,
            passed_simple_repeat,
            passed_all,
            filter_efficiency: passed_all as f64 / total_candidates as f64,
        }
    }
}

/// Filter statistics for analysis
#[derive(Debug, Clone)]
pub struct FilterStatistics {
    pub total_candidates: usize,
    pub passed_frequency: usize,
    pub passed_entropy: usize,
    pub passed_gc: usize,
    pub passed_homopolymer: usize,
    pub passed_simple_repeat: usize,
    pub passed_all: usize,
    pub filter_efficiency: f64,
}

/// Frequency-based k-mer processor for managing computational load
pub struct FrequencyProcessor {
    max_candidates_per_stratum: usize,
    max_total_candidates: usize,
}

impl Default for FrequencyProcessor {
    fn default() -> Self {
        Self {
            max_candidates_per_stratum: 1000,  // Limit per quality level
            max_total_candidates: 5000,        // Total processing limit
        }
    }
}

impl FrequencyProcessor {
    pub fn new(max_per_stratum: usize, max_total: usize) -> Self {
        Self {
            max_candidates_per_stratum: max_per_stratum,
            max_total_candidates: max_total,
        }
    }

    /// Process candidates in frequency-stratified manner
    pub fn process_stratified_candidates(
        &self,
        strata: &[Vec<(Vec<u8>, usize)>],
    ) -> Vec<(Vec<u8>, usize)> {
        let mut processed_candidates = Vec::new();
        let mut total_processed = 0;

        // Process high-quality strata first
        for stratum in strata {
            let candidates_to_take = std::cmp::min(
                self.max_candidates_per_stratum,
                self.max_total_candidates - total_processed,
            );

            if candidates_to_take == 0 {
                break;
            }

            let stratum_candidates = stratum
                .iter()
                .take(candidates_to_take)
                .cloned()
                .collect::<Vec<_>>();

            total_processed += stratum_candidates.len();
            processed_candidates.extend(stratum_candidates);

            if total_processed >= self.max_total_candidates {
                break;
            }
        }

        processed_candidates
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_entropy_calculation() {
        let filter = QualityFilter::default();
        
        // High entropy sequence (all different bases)
        let high_entropy = vec![0, 1, 2, 3]; // ACGT
        let entropy = filter.compute_entropy(&high_entropy);
        assert!(entropy > 1.5); // Should be near 2.0
        
        // Low entropy sequence (all same base)
        let low_entropy = vec![0, 0, 0, 0]; // AAAA
        let entropy = filter.compute_entropy(&low_entropy);
        assert!(entropy < 0.1); // Should be near 0.0
    }

    #[test]
    fn test_gc_content_filtering() {
        let filter = QualityFilter::default();
        
        // Normal GC content (50%)
        let normal_gc = vec![0, 1, 2, 3]; // ACGT = 50% GC
        assert!(filter.check_gc_content(&normal_gc));
        
        // Extreme GC content (100% GC)
        let high_gc = vec![1, 2, 1, 2]; // CGCG = 100% GC
        assert!(!filter.check_gc_content(&high_gc));
        
        // No GC content (0% GC)
        let no_gc = vec![0, 3, 0, 3]; // ATAT = 0% GC
        assert!(!filter.check_gc_content(&no_gc));
    }

    #[test]
    fn test_homopolymer_detection() {
        let filter = QualityFilter::default();
        
        // Normal sequence
        let normal = vec![0, 1, 2, 3, 0, 1]; // ACGTAC
        assert!(!filter.has_excessive_homopolymer(&normal));
        
        // Long homopolymer run
        let homopolymer = vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; // 10 A's
        assert!(filter.has_excessive_homopolymer(&homopolymer));
    }

    #[test]
    fn test_prefiltering() {
        let filter = QualityFilter::default();
        
        // High-quality k-mer
        let good_kmer = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 1]; // Mixed sequence
        assert!(filter.prefilter_kmer(&good_kmer, 10));
        
        // Low-frequency k-mer
        assert!(!filter.prefilter_kmer(&good_kmer, 1));
        
        // Low-entropy k-mer
        let bad_kmer = vec![0; 13]; // All A's
        assert!(!filter.prefilter_kmer(&bad_kmer, 10));
    }
}