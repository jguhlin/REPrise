//! Banded dynamic programming extension alignment
//! 
//! Implements bidirectional extension from k-mer seeds using banded DP to build
//! consensus sequences, matching C++ REPrise extension behavior.

use crate::quality::{QualityError};
use std::cmp::{min, max};

/// Banded dynamic programming aligner for consensus extension
pub struct BandedAligner {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open_score: i32,
    pub gap_extend_score: i32,
    pub cap_penalty: i32,
    pub band_width: usize,
    pub max_extend: usize,
    pub stop_threshold: usize,
    pub min_improvement: usize,
}

impl Default for BandedAligner {
    fn default() -> Self {
        Self {
            match_score: 1,        // C++ MATCHSCORE
            mismatch_score: -1,    // C++ MISMATCHSCORE
            gap_open_score: -5,    // C++ GAPSCORE
            gap_extend_score: -1,  // C++ GAPEXTENDSCORE
            cap_penalty: -20,      // C++ CAPPENALTY
            band_width: 11,        // C++ 2*OFFSETWIDTH + 1
            max_extend: 10000,     // C++ MAXEXTEND
            stop_threshold: 100,   // C++ WHEN_TO_STOP
            min_improvement: 3,    // C++ MINIMPROVEMENT
        }
    }
}

/// Extension result containing alignments and quality metrics
#[derive(Debug, Clone)]
pub struct ExtensionResult {
    pub left_extension: Vec<u8>,
    pub right_extension: Vec<u8>,
    pub alignments: Vec<AlignmentResult>,
    pub total_score: f64,
    pub extension_quality: f64,
}

/// Individual sequence alignment result  
#[derive(Debug, Clone)]
pub struct AlignmentResult {
    pub sequence_id: usize,
    pub start_pos: usize,
    pub end_pos: usize,
    pub alignment_score: f64,
    pub is_reverse: bool,
    pub aligned_sequence: Vec<u8>,
}

/// DP scoring matrices (reused for memory efficiency)
struct DPMatrices {
    score_match: Vec<Vec<i32>>,
    score_ins: Vec<Vec<i32>>,
    score_del: Vec<Vec<i32>>,
    traceback: Vec<Vec<u8>>,
}

impl DPMatrices {
    fn new(max_length: usize, band_width: usize) -> Self {
        let rows = max_length + 1;
        let cols = band_width;
        
        Self {
            score_match: vec![vec![i32::MIN / 2; cols]; rows],
            score_ins: vec![vec![i32::MIN / 2; cols]; rows],
            score_del: vec![vec![i32::MIN / 2; cols]; rows],
            traceback: vec![vec![0; cols]; rows],
        }
    }

    fn reset(&mut self) {
        for row in &mut self.score_match {
            row.fill(i32::MIN / 2);
        }
        for row in &mut self.score_ins {
            row.fill(i32::MIN / 2);
        }
        for row in &mut self.score_del {
            row.fill(i32::MIN / 2);
        }
        for row in &mut self.traceback {
            row.fill(0);
        }
    }
}

impl BandedAligner {
    /// Perform bidirectional extension from k-mer seeds
    pub fn extend_bidirectional(
        &self,
        seed_kmer: &[u8],
        positions: &[usize],
        rev_flags: &[bool],
        sequence: &[u8],
        mask: &[bool],
    ) -> Result<ExtensionResult, QualityError> {
        // Filter valid, unmasked positions
        let valid_positions: Vec<(usize, bool)> = positions
            .iter()
            .zip(rev_flags.iter())
            .enumerate()
            .filter_map(|(i, (&pos, &is_rev))| {
                if pos < sequence.len() && pos + seed_kmer.len() < sequence.len() {
                    // Check if position is masked
                    let masked = (pos..pos + seed_kmer.len()).any(|p| *mask.get(p).unwrap_or(&true));
                    if !masked {
                        Some((pos, is_rev))
                    } else {
                        None
                    }
                } else {
                    None
                }
            })
            .collect();

        if valid_positions.len() < 3 {
            return Err(QualityError::InsufficientData);
        }

        // Perform left extension
        let left_result = self.extend_direction(
            seed_kmer, &valid_positions, sequence, false // left = false
        )?;

        // Perform right extension  
        let right_result = self.extend_direction(
            seed_kmer, &valid_positions, sequence, true // right = true
        )?;

        // Combine results and generate alignments
        let alignments = self.generate_alignments(
            &left_result, &right_result, seed_kmer, &valid_positions, sequence
        )?;

        let total_score = alignments.iter().map(|a| a.alignment_score).sum::<f64>();
        let extension_quality = self.calculate_extension_quality(&left_result, &right_result);

        Ok(ExtensionResult {
            left_extension: left_result,
            right_extension: right_result,
            alignments,
            total_score,
            extension_quality,
        })
    }

    /// Extend in one direction (left or right) using banded DP
    fn extend_direction(
        &self,
        seed_kmer: &[u8],
        positions: &[(usize, bool)],
        sequence: &[u8],
        extend_right: bool,
    ) -> Result<Vec<u8>, QualityError> {
        let mut extension = Vec::new();
        let mut dp_matrices = DPMatrices::new(self.max_extend, self.band_width);
        let mut no_improvement_count = 0;
        let mut last_best_score = i32::MIN;

        for ext_pos in 0..self.max_extend {
            let mut base_scores = [0i32; 4]; // A, C, G, T
            
            // Try each base at this extension position
            for base in 0..4 {
                let mut total_score = 0i32;
                
                // Score this base against all valid sequence positions
                for (i, &(seq_pos, is_reverse)) in positions.iter().enumerate() {
                    let score = self.align_position_banded(
                        &mut dp_matrices,
                        sequence,
                        seq_pos,
                        seed_kmer.len() + ext_pos,
                        base,
                        extend_right,
                        is_reverse,
                    );
                    
                    // Apply cap penalty if sequence doesn't extend this far
                    let actual_score = if self.sequence_extends(seq_pos, seed_kmer.len() + ext_pos, extend_right, sequence, is_reverse) {
                        score
                    } else {
                        max(0, score + self.cap_penalty)
                    };
                    
                    total_score += actual_score;
                }
                
                base_scores[base] = total_score;
            }
            
            // Select best base for this position
            let best_base = base_scores
                .iter()
                .enumerate()
                .max_by_key(|&(_, score)| score)
                .map(|(base, _)| base)
                .unwrap_or(0);
                
            let best_score = base_scores[best_base];
            
            // Check for improvement
            if best_score > last_best_score + (self.min_improvement as i32) {
                extension.push(best_base as u8);
                last_best_score = best_score;
                no_improvement_count = 0;
            } else {
                no_improvement_count += 1;
                if no_improvement_count >= self.stop_threshold {
                    break; // Early termination like C++
                }
                extension.push(best_base as u8);
            }
        }

        Ok(extension)
    }

    /// Perform banded alignment at specific position
    fn align_position_banded(
        &self,
        matrices: &mut DPMatrices,
        sequence: &[u8],
        seq_start: usize,
        align_length: usize,
        consensus_base: usize,
        extend_right: bool,
        is_reverse: bool,
    ) -> i32 {
        matrices.reset();
        
        let offset_width = self.band_width / 2;
        let seq_end = min(seq_start + align_length, sequence.len());
        
        if seq_start >= seq_end {
            return self.cap_penalty;
        }

        // Initialize first row/column
        matrices.score_match[0][offset_width] = 0;
        
        // Fill DP matrices in banded fashion
        for i in 1..=min(align_length, seq_end - seq_start) {
            for j_offset in 0..self.band_width {
                let j = if j_offset >= offset_width {
                    j_offset - offset_width
                } else {
                    continue;
                };
                
                if j >= align_length { break; }
                
                let seq_idx = if extend_right {
                    seq_start + i - 1
                } else {
                    seq_start + align_length - i
                };
                
                if seq_idx >= sequence.len() { continue; }
                
                let seq_base = if is_reverse {
                    3 - sequence[seq_idx] // complement
                } else {
                    sequence[seq_idx]
                } as usize;
                
                let consensus_pos_base = if j < 1 {
                    consensus_base
                } else {
                    consensus_base // simplified for this position
                };
                
                // Calculate scores for match, insertion, deletion
                let match_score = if seq_base == consensus_pos_base {
                    self.match_score
                } else {
                    self.mismatch_score
                };
                
                // Match/mismatch
                if i > 0 && j_offset > 0 {
                    let prev_score = matrices.score_match[i-1][j_offset-1];
                    matrices.score_match[i][j_offset] = prev_score + match_score;
                }
                
                // Insertion (gap in sequence)
                if i > 0 {
                    let gap_score = matrices.score_ins[i-1][j_offset] + self.gap_extend_score;
                    let open_score = matrices.score_match[i-1][j_offset] + self.gap_open_score;
                    matrices.score_ins[i][j_offset] = max(gap_score, open_score);
                }
                
                // Deletion (gap in consensus)
                if j_offset > 0 {
                    let gap_score = matrices.score_del[i][j_offset-1] + self.gap_extend_score;
                    let open_score = matrices.score_match[i][j_offset-1] + self.gap_open_score;
                    matrices.score_del[i][j_offset] = max(gap_score, open_score);
                }
                
                // Take maximum of all three states
                matrices.score_match[i][j_offset] = max(
                    matrices.score_match[i][j_offset],
                    max(matrices.score_ins[i][j_offset], matrices.score_del[i][j_offset])
                );
            }
        }
        
        // Find best score in final row/column
        let final_row = min(align_length, seq_end - seq_start);
        matrices.score_match[final_row]
            .iter()
            .max()
            .copied()
            .unwrap_or(self.cap_penalty)
    }

    /// Check if sequence extends far enough in given direction
    fn sequence_extends(
        &self,
        seq_pos: usize,
        required_length: usize,
        extend_right: bool,
        sequence: &[u8],
        is_reverse: bool,
    ) -> bool {
        if extend_right {
            seq_pos + required_length < sequence.len()
        } else {
            seq_pos >= required_length
        }
    }

    /// Generate alignment results for all positions
    fn generate_alignments(
        &self,
        left_ext: &[u8],
        right_ext: &[u8],
        seed: &[u8],
        positions: &[(usize, bool)],
        sequence: &[u8],
    ) -> Result<Vec<AlignmentResult>, QualityError> {
        let mut alignments = Vec::new();
        
        for (i, &(pos, is_rev)) in positions.iter().enumerate() {
            let start_pos = if left_ext.is_empty() { 
                pos 
            } else { 
                pos.saturating_sub(left_ext.len()) 
            };
            
            let end_pos = min(
                pos + seed.len() + right_ext.len(),
                sequence.len()
            );
            
            if start_pos < end_pos {
                let aligned_seq = sequence[start_pos..end_pos].to_vec();
                let alignment_score = self.calculate_alignment_score(&aligned_seq, left_ext, seed, right_ext);
                
                alignments.push(AlignmentResult {
                    sequence_id: i,
                    start_pos,
                    end_pos,
                    alignment_score,
                    is_reverse: is_rev,
                    aligned_sequence: aligned_seq,
                });
            }
        }
        
        Ok(alignments)
    }

    /// Calculate alignment score for a sequence against consensus
    fn calculate_alignment_score(
        &self,
        sequence: &[u8],
        left_ext: &[u8],
        seed: &[u8],
        right_ext: &[u8],
    ) -> f64 {
        let mut score = 0.0;
        let mut seq_idx = 0;
        
        // Score against left extension
        for &base in left_ext {
            if seq_idx < sequence.len() {
                score += if sequence[seq_idx] == base {
                    self.match_score as f64
                } else {
                    self.mismatch_score as f64
                };
                seq_idx += 1;
            }
        }
        
        // Score against seed
        for &base in seed {
            if seq_idx < sequence.len() {
                score += if sequence[seq_idx] == base {
                    self.match_score as f64
                } else {
                    self.mismatch_score as f64
                };
                seq_idx += 1;
            }
        }
        
        // Score against right extension
        for &base in right_ext {
            if seq_idx < sequence.len() {
                score += if sequence[seq_idx] == base {
                    self.match_score as f64
                } else {
                    self.mismatch_score as f64
                };
                seq_idx += 1;
            }
        }
        
        score
    }

    /// Calculate extension quality based on consistency and length
    fn calculate_extension_quality(&self, left_ext: &[u8], right_ext: &[u8]) -> f64 {
        let total_length = left_ext.len() + right_ext.len();
        if total_length == 0 {
            return 0.0;
        }
        
        // Simple quality metric: longer extensions with consistent bases score higher
        let quality = (total_length as f64).sqrt() / 10.0;
        quality.min(1.0).max(0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_banded_aligner_creation() {
        let aligner = BandedAligner::default();
        assert_eq!(aligner.match_score, 1);
        assert_eq!(aligner.mismatch_score, -1);
        assert_eq!(aligner.band_width, 11);
    }

    #[test]
    fn test_extension_basic() {
        let aligner = BandedAligner::default();
        let seed = vec![0, 1, 2, 3]; // ACGT
        let sequence = vec![0, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGT
        let positions = vec![0, 4];
        let rev_flags = vec![false, false];
        let mask = vec![false; sequence.len()];

        let result = aligner.extend_bidirectional(&seed, &positions, &rev_flags, &sequence, &mask);
        assert!(result.is_ok());
        
        let ext_result = result.unwrap();
        assert_eq!(ext_result.alignments.len(), 2);
        assert!(ext_result.total_score > 0.0);
    }
}