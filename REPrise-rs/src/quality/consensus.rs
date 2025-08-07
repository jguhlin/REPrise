//! Multi-way consensus sequence building
//! 
//! Builds consensus sequences from multiple alignment results using frequency-based
//! base calling and quality assessment.

use crate::quality::{QualityError, ExtensionResult, AlignmentResult};
use std::collections::HashMap;

/// Consensus sequence builder using multi-way alignment
pub struct ConsensusBuilder {
    pub min_consensus_freq: f64,  // Minimum frequency for consensus base calling
    pub min_quality_score: f64,   // Minimum per-position quality score
    pub max_gap_fraction: f64,    // Maximum fraction of gaps allowed
    pub min_coverage: usize,      // Minimum number of sequences covering position
}

impl Default for ConsensusBuilder {
    fn default() -> Self {
        Self {
            min_consensus_freq: 0.5,   // 50% agreement for consensus
            min_quality_score: 0.3,    // Minimum position quality
            max_gap_fraction: 0.3,     // Max 30% gaps
            min_coverage: 3,           // At least 3 sequences
        }
    }
}

/// Consensus sequence with quality metrics
#[derive(Debug, Clone)]
pub struct ConsensusSequence {
    pub sequence: Vec<u8>,
    pub left_extension: Vec<u8>,
    pub seed_sequence: Vec<u8>,
    pub right_extension: Vec<u8>,
    pub position_qualities: Vec<f64>,
    pub overall_quality: f64,
    pub coverage: Vec<usize>,
    pub length: usize,
    pub gap_positions: Vec<bool>,
}

/// Position-specific consensus information
#[derive(Debug)]
struct PositionConsensus {
    base_counts: [usize; 4], // A, C, G, T counts
    gap_count: usize,
    total_coverage: usize,
    consensus_base: u8,
    quality_score: f64,
}

impl ConsensusBuilder {
    pub fn new(min_consensus_freq: f64, min_quality_score: f64) -> Self {
        Self {
            min_consensus_freq,
            min_quality_score,
            max_gap_fraction: 0.3,
            min_coverage: 3,
        }
    }

    /// Build consensus sequence from extension result
    pub fn build_consensus(
        &self,
        extension_result: &ExtensionResult,
    ) -> Result<ConsensusSequence, QualityError> {
        if extension_result.alignments.len() < self.min_coverage {
            return Err(QualityError::InsufficientData);
        }

        // Extract sequences and build consensus
        let aligned_sequences = self.extract_aligned_sequences(&extension_result.alignments)?;
        let consensus_positions = self.build_position_consensus(&aligned_sequences)?;
        
        // Split into left extension, seed, right extension
        let (left_ext, seed_seq, right_ext) = self.split_consensus_regions(
            &consensus_positions,
            &extension_result
        )?;

        // Calculate overall quality metrics
        let overall_quality = self.calculate_overall_quality(&consensus_positions);
        let position_qualities: Vec<f64> = consensus_positions
            .iter()
            .map(|pos| pos.quality_score)
            .collect();
        
        let coverage: Vec<usize> = consensus_positions
            .iter()
            .map(|pos| pos.total_coverage)
            .collect();

        let gap_positions: Vec<bool> = consensus_positions
            .iter()
            .map(|pos| pos.gap_count > pos.total_coverage / 2)
            .collect();

        // Build final sequence
        let sequence: Vec<u8> = consensus_positions
            .iter()
            .map(|pos| pos.consensus_base)
            .collect();

        Ok(ConsensusSequence {
            sequence: sequence.clone(),
            left_extension: left_ext,
            seed_sequence: seed_seq,
            right_extension: right_ext,
            position_qualities,
            overall_quality,
            coverage,
            length: sequence.len(),
            gap_positions,
        })
    }

    /// Extract aligned sequences from alignment results
    fn extract_aligned_sequences(
        &self,
        alignments: &[AlignmentResult],
    ) -> Result<Vec<Vec<u8>>, QualityError> {
        if alignments.is_empty() {
            return Err(QualityError::InsufficientData);
        }

        // Find maximum sequence length for alignment
        let max_length = alignments
            .iter()
            .map(|a| a.aligned_sequence.len())
            .max()
            .unwrap_or(0);

        if max_length == 0 {
            return Err(QualityError::ConsensusBuilding("No sequence data".to_string()));
        }

        // Pad sequences to same length and handle reverse complements
        let mut aligned_sequences = Vec::new();
        
        for alignment in alignments {
            let mut seq = alignment.aligned_sequence.clone();
            
            // Apply reverse complement if needed
            if alignment.is_reverse {
                seq = self.reverse_complement(&seq);
            }
            
            // Pad to maximum length
            while seq.len() < max_length {
                seq.push(4); // Use 4 as gap character
            }
            
            aligned_sequences.push(seq);
        }

        Ok(aligned_sequences)
    }

    /// Build position-wise consensus from aligned sequences
    fn build_position_consensus(
        &self,
        aligned_sequences: &[Vec<u8>],
    ) -> Result<Vec<PositionConsensus>, QualityError> {
        if aligned_sequences.is_empty() {
            return Err(QualityError::InsufficientData);
        }

        let seq_length = aligned_sequences[0].len();
        let mut consensus_positions = Vec::new();

        for pos in 0..seq_length {
            let mut base_counts = [0; 4]; // A, C, G, T
            let mut gap_count = 0;
            let mut total_coverage = 0;

            // Count bases at this position across all sequences
            for sequence in aligned_sequences {
                if pos < sequence.len() {
                    let base = sequence[pos];
                    if base < 4 {
                        base_counts[base as usize] += 1;
                        total_coverage += 1;
                    } else {
                        gap_count += 1;
                    }
                }
            }

            // Skip positions with insufficient coverage
            if total_coverage < self.min_coverage {
                continue;
            }

            // Find consensus base (most frequent)
            let (consensus_base, max_count) = base_counts
                .iter()
                .enumerate()
                .max_by_key(|&(_, count)| count)
                .map(|(base, &count)| (base as u8, count))
                .unwrap_or((0, 0));

            // Calculate position quality
            let consensus_frequency = max_count as f64 / total_coverage as f64;
            let gap_fraction = gap_count as f64 / (total_coverage + gap_count) as f64;
            
            // Quality based on consensus frequency and gap content
            let quality_score = consensus_frequency * (1.0 - gap_fraction);
            
            // Only include high-quality positions
            if consensus_frequency >= self.min_consensus_freq 
                && gap_fraction <= self.max_gap_fraction
                && quality_score >= self.min_quality_score {
                
                consensus_positions.push(PositionConsensus {
                    base_counts,
                    gap_count,
                    total_coverage,
                    consensus_base,
                    quality_score,
                });
            }
        }

        if consensus_positions.is_empty() {
            return Err(QualityError::ConsensusBuilding("No high-quality positions".to_string()));
        }

        Ok(consensus_positions)
    }

    /// Split consensus into left extension, seed, and right extension regions
    fn split_consensus_regions(
        &self,
        consensus_positions: &[PositionConsensus],
        extension_result: &ExtensionResult,
    ) -> Result<(Vec<u8>, Vec<u8>, Vec<u8>), QualityError> {
        let total_length = consensus_positions.len();
        let left_length = extension_result.left_extension.len();
        let right_length = extension_result.right_extension.len();
        
        // Estimate seed length (total - extensions, or use original seed if available)
        let seed_length = if total_length > left_length + right_length {
            total_length - left_length - right_length
        } else {
            total_length / 3 // Fallback: assume seed is roughly 1/3
        };

        let mut left_ext = Vec::new();
        let mut seed_seq = Vec::new();
        let mut right_ext = Vec::new();

        for (i, pos_consensus) in consensus_positions.iter().enumerate() {
            if i < left_length {
                left_ext.push(pos_consensus.consensus_base);
            } else if i < left_length + seed_length {
                seed_seq.push(pos_consensus.consensus_base);
            } else {
                right_ext.push(pos_consensus.consensus_base);
            }
        }

        Ok((left_ext, seed_seq, right_ext))
    }

    /// Calculate overall consensus quality
    fn calculate_overall_quality(&self, consensus_positions: &[PositionConsensus]) -> f64 {
        if consensus_positions.is_empty() {
            return 0.0;
        }

        let avg_position_quality: f64 = consensus_positions
            .iter()
            .map(|pos| pos.quality_score)
            .sum::<f64>() / consensus_positions.len() as f64;

        let avg_coverage: f64 = consensus_positions
            .iter()
            .map(|pos| pos.total_coverage)
            .sum::<usize>() as f64 / consensus_positions.len() as f64;

        let length_factor = (consensus_positions.len() as f64).sqrt() / 10.0;
        let coverage_factor = (avg_coverage / self.min_coverage as f64).min(2.0);

        // Composite quality score
        (avg_position_quality * coverage_factor * length_factor).min(1.0)
    }

    /// Generate reverse complement of sequence
    fn reverse_complement(&self, sequence: &[u8]) -> Vec<u8> {
        sequence
            .iter()
            .rev()
            .map(|&base| match base {
                0 => 3, // A -> T
                1 => 2, // C -> G  
                2 => 1, // G -> C
                3 => 0, // T -> A
                _ => base, // Gap or unknown
            })
            .collect()
    }

    /// Validate consensus sequence meets quality thresholds
    pub fn validate_consensus(&self, consensus: &ConsensusSequence) -> bool {
        // Check minimum length
        if consensus.length < 50 {
            return false;
        }

        // Check overall quality
        if consensus.overall_quality < self.min_quality_score {
            return false;
        }

        // Check gap content
        let gap_fraction = consensus.gap_positions.iter().filter(|&&is_gap| is_gap).count() as f64 
            / consensus.gap_positions.len() as f64;
        
        if gap_fraction > self.max_gap_fraction {
            return false;
        }

        // Check minimum coverage
        let avg_coverage = consensus.coverage.iter().sum::<usize>() as f64 / consensus.coverage.len() as f64;
        if avg_coverage < self.min_coverage as f64 {
            return false;
        }

        true
    }

    /// Get detailed consensus statistics
    pub fn get_consensus_stats(&self, consensus: &ConsensusSequence) -> ConsensusStatistics {
        let gc_content = consensus.sequence
            .iter()
            .filter(|&&base| base == 1 || base == 2) // C or G
            .count() as f64 / consensus.sequence.len() as f64;

        let avg_quality = consensus.position_qualities.iter().sum::<f64>() 
            / consensus.position_qualities.len() as f64;
        
        let avg_coverage = consensus.coverage.iter().sum::<usize>() as f64 
            / consensus.coverage.len() as f64;

        let gap_fraction = consensus.gap_positions.iter().filter(|&&is_gap| is_gap).count() as f64 
            / consensus.gap_positions.len() as f64;

        ConsensusStatistics {
            length: consensus.length,
            gc_content,
            avg_quality,
            avg_coverage,
            gap_fraction,
            overall_quality: consensus.overall_quality,
        }
    }
}

/// Consensus sequence statistics
#[derive(Debug, Clone)]
pub struct ConsensusStatistics {
    pub length: usize,
    pub gc_content: f64,
    pub avg_quality: f64,
    pub avg_coverage: f64,
    pub gap_fraction: f64,
    pub overall_quality: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quality::extension::*;

    #[test]
    fn test_consensus_builder_creation() {
        let builder = ConsensusBuilder::default();
        assert_eq!(builder.min_consensus_freq, 0.5);
        assert_eq!(builder.min_coverage, 3);
    }

    #[test]
    fn test_reverse_complement() {
        let builder = ConsensusBuilder::default();
        let seq = vec![0, 1, 2, 3]; // ACGT
        let rev_comp = builder.reverse_complement(&seq);
        assert_eq!(rev_comp, vec![3, 2, 1, 0]); // TGCA
    }

    #[test]
    fn test_consensus_validation() {
        let builder = ConsensusBuilder::default();
        let consensus = ConsensusSequence {
            sequence: vec![0; 60], // 60 bp sequence
            left_extension: vec![],
            seed_sequence: vec![0; 60],
            right_extension: vec![],
            position_qualities: vec![0.8; 60],
            overall_quality: 0.7,
            coverage: vec![5; 60],
            length: 60,
            gap_positions: vec![false; 60],
        };
        
        assert!(builder.validate_consensus(&consensus));
    }
}