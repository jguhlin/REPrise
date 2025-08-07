//! Progressive masking system with boundary refinement
//! 
//! Implements sophisticated masking strategy to prevent overlapping family calls
//! while maintaining optimal coverage, matching C++ masking_align behavior.

use crate::quality::{QualityError, ConsensusSequence, ExtensionResult, AlignmentResult, QualityMetrics};
use std::collections::HashMap;
use std::cmp::{min, max};

/// Progressive masker with boundary refinement and overlap resolution
pub struct ProgressiveMasker {
    pub tandem_dist: usize,        // C++ TANDEMDIST = 500
    pub mask_buffer: usize,        // Buffer around masked regions
    pub min_overlap_penalty: f64,  // Penalty for overlapping families
    pub quality_priority: bool,    // Use quality for overlap resolution
}

impl Default for ProgressiveMasker {
    fn default() -> Self {
        Self {
            tandem_dist: 500,           // Match C++ TANDEMDIST
            mask_buffer: 10,            // Small buffer to prevent edge effects
            min_overlap_penalty: 0.5,   // Reduce quality for overlaps
            quality_priority: true,     // Higher quality families take precedence
        }
    }
}

/// Masked element with refined boundaries
#[derive(Debug, Clone)]
pub struct MaskedElement {
    pub start_pos: usize,
    pub end_pos: usize,
    pub alignment_score: f64,
    pub is_reverse: bool,
    pub sequence_id: usize,
    pub refined_boundaries: bool,
}

/// Overlap resolution result
#[derive(Debug)]
struct OverlapResolution {
    pub keep_current: bool,
    pub conflict_regions: Vec<(usize, usize)>,
    pub quality_adjustment: f64,
}

impl ProgressiveMasker {
    pub fn new(tandem_dist: usize) -> Self {
        Self {
            tandem_dist,
            ..Default::default()
        }
    }

    /// Mask family elements with boundary refinement and overlap resolution
    pub fn mask_by_family(
        &mut self,
        consensus: &ConsensusSequence,
        extension_result: &ExtensionResult,
        mask: &mut [bool],
    ) -> Result<Vec<MaskedElement>, QualityError> {
        let mut masked_elements = Vec::new();

        // Step 1: Refine element boundaries using consensus alignment
        let refined_elements = self.refine_element_boundaries(
            consensus, &extension_result.alignments, mask
        )?;

        // Step 2: Resolve overlaps with existing masked regions
        let resolved_elements = self.resolve_mask_overlaps(&refined_elements, mask)?;

        // Step 3: Apply progressive masking with quality priority
        for element in resolved_elements {
            if self.should_mask_element(&element, mask) {
                self.apply_element_mask(&element, mask)?;
                masked_elements.push(element);
            }
        }

        // Step 4: Validate masking consistency
        self.validate_masking_consistency(&masked_elements, mask)?;

        Ok(masked_elements)
    }

    /// Refine element boundaries using individual consensus alignment (like C++ masking_align)
    fn refine_element_boundaries(
        &self,
        consensus: &ConsensusSequence,
        alignments: &[AlignmentResult],
        current_mask: &[bool],
    ) -> Result<Vec<MaskedElement>, QualityError> {
        let mut refined_elements = Vec::new();

        for (i, alignment) in alignments.iter().enumerate() {
            // Skip if element region is already heavily masked
            if self.is_region_heavily_masked(alignment.start_pos, alignment.end_pos, current_mask) {
                continue;
            }

            // Perform individual alignment to refine boundaries
            let refined_boundaries = self.perform_individual_alignment(
                consensus, alignment, current_mask
            )?;

            refined_elements.push(MaskedElement {
                start_pos: refined_boundaries.0,
                end_pos: refined_boundaries.1,
                alignment_score: alignment.alignment_score,
                is_reverse: alignment.is_reverse,
                sequence_id: alignment.sequence_id,
                refined_boundaries: true,
            });
        }

        Ok(refined_elements)
    }

    /// Perform individual element alignment for boundary refinement (C++ masking_align equivalent)
    fn perform_individual_alignment(
        &self,
        consensus: &ConsensusSequence,
        alignment: &AlignmentResult,
        current_mask: &[bool],
    ) -> Result<(usize, usize), QualityError> {
        let consensus_seq = &consensus.sequence;
        let element_seq = &alignment.aligned_sequence;

        // Use banded DP for precise boundary detection
        let alignment_result = self.align_with_consensus(
            consensus_seq, element_seq, alignment.start_pos, alignment.is_reverse
        )?;

        // Adjust boundaries based on alignment quality
        let (start_pos, end_pos) = self.optimize_boundaries(
            &alignment_result, alignment.start_pos, alignment.end_pos, current_mask
        );

        // Ensure boundaries are within valid range
        let start_pos = max(0, start_pos);
        let end_pos = min(current_mask.len(), end_pos);

        if start_pos >= end_pos {
            return Err(QualityError::MaskingError("Invalid boundaries after refinement".to_string()));
        }

        Ok((start_pos, end_pos))
    }

    /// Align element sequence with consensus using simplified banded DP
    fn align_with_consensus(
        &self,
        consensus: &[u8],
        element: &[u8],
        start_hint: usize,
        is_reverse: bool,
    ) -> Result<AlignmentMatrix, QualityError> {
        let consensus_len = consensus.len();
        let element_len = element.len();
        
        if consensus_len == 0 || element_len == 0 {
            return Err(QualityError::MaskingError("Empty sequences for alignment".to_string()));
        }

        // Simple DP matrix for boundary refinement
        let mut dp = vec![vec![0i32; element_len + 1]; consensus_len + 1];
        
        // Initialize gap penalties
        for i in 1..=consensus_len {
            dp[i][0] = -(i as i32) * 2; // Gap penalty
        }
        for j in 1..=element_len {
            dp[0][j] = -(j as i32) * 2; // Gap penalty
        }

        // Fill DP matrix
        for i in 1..=consensus_len {
            for j in 1..=element_len {
                let consensus_base = if is_reverse {
                    3 - consensus[consensus_len - i] // Reverse complement
                } else {
                    consensus[i - 1]
                };
                
                let element_base = element[j - 1];
                
                let match_score = if consensus_base == element_base { 2 } else { -1 };
                
                let diagonal = dp[i-1][j-1] + match_score;
                let up = dp[i-1][j] - 2; // Gap in element
                let left = dp[i][j-1] - 2; // Gap in consensus
                
                dp[i][j] = max(max(diagonal, up), left);
            }
        }

        // Find best alignment boundaries
        let (best_i, best_j) = self.find_best_alignment_end(&dp);
        
        Ok(AlignmentMatrix {
            matrix: dp,
            best_end: (best_i, best_j),
            consensus_len,
            element_len,
        })
    }

    /// Find optimal alignment endpoint
    fn find_best_alignment_end(&self, dp: &[Vec<i32>]) -> (usize, usize) {
        let rows = dp.len();
        let cols = dp[0].len();
        let mut best_score = i32::MIN;
        let mut best_pos = (0, 0);

        // Check last row and column for best local alignment end
        for j in 0..cols {
            if dp[rows-1][j] > best_score {
                best_score = dp[rows-1][j];
                best_pos = (rows-1, j);
            }
        }
        
        for i in 0..rows {
            if dp[i][cols-1] > best_score {
                best_score = dp[i][cols-1];
                best_pos = (i, cols-1);
            }
        }

        best_pos
    }

    /// Optimize element boundaries based on alignment quality
    fn optimize_boundaries(
        &self,
        alignment: &AlignmentMatrix,
        original_start: usize,
        original_end: usize,
        current_mask: &[bool],
    ) -> (usize, usize) {
        // Traceback from best alignment end to find optimal boundaries
        let (end_i, end_j) = alignment.best_end;
        let (start_i, start_j) = self.traceback_to_start(alignment, end_i, end_j);
        
        // Convert matrix coordinates to genome coordinates
        let genome_start = original_start + start_i;
        let genome_end = original_start + end_i;
        
        // Adjust for masked regions (avoid extending into masked areas)
        let adjusted_start = self.adjust_boundary_for_mask(genome_start, current_mask, true);
        let adjusted_end = self.adjust_boundary_for_mask(genome_end, current_mask, false);
        
        (adjusted_start, adjusted_end)
    }

    /// Traceback through alignment matrix to find start position
    fn traceback_to_start(&self, alignment: &AlignmentMatrix, end_i: usize, end_j: usize) -> (usize, usize) {
        let mut i = end_i;
        let mut j = end_j;
        
        // Simple traceback to find reasonable start (could be more sophisticated)
        while i > 0 && j > 0 && alignment.matrix[i][j] > 0 {
            // Find direction of best predecessor
            let diagonal = if i > 0 && j > 0 { alignment.matrix[i-1][j-1] } else { i32::MIN };
            let up = if i > 0 { alignment.matrix[i-1][j] } else { i32::MIN };
            let left = if j > 0 { alignment.matrix[i][j-1] } else { i32::MIN };
            
            if diagonal >= up && diagonal >= left && i > 0 && j > 0 {
                i -= 1;
                j -= 1;
            } else if up >= left && i > 0 {
                i -= 1;
            } else if j > 0 {
                j -= 1;
            } else {
                break;
            }
        }
        
        (i, j)
    }

    /// Adjust boundary to avoid masked regions
    fn adjust_boundary_for_mask(&self, position: usize, mask: &[bool], is_start: bool) -> usize {
        if position >= mask.len() {
            return if is_start { mask.len().saturating_sub(1) } else { mask.len() };
        }

        let mut adjusted_pos = position;
        
        if is_start {
            // For start boundary, move right to avoid masked region
            while adjusted_pos < mask.len() && mask[adjusted_pos] {
                adjusted_pos += 1;
            }
        } else {
            // For end boundary, move left to avoid masked region
            while adjusted_pos > 0 && mask[adjusted_pos.saturating_sub(1)] {
                adjusted_pos -= 1;
            }
        }
        
        adjusted_pos
    }

    /// Resolve overlaps between new elements and existing mask
    fn resolve_mask_overlaps(
        &self,
        elements: &[MaskedElement],
        current_mask: &[bool],
    ) -> Result<Vec<MaskedElement>, QualityError> {
        let mut resolved_elements = Vec::new();

        for element in elements {
            let overlap_resolution = self.assess_mask_overlap(element, current_mask);
            
            if overlap_resolution.keep_current {
                // Apply quality adjustment for partial overlaps
                let mut adjusted_element = element.clone();
                adjusted_element.alignment_score *= (1.0 - overlap_resolution.quality_adjustment);
                resolved_elements.push(adjusted_element);
            }
            // If keep_current is false, element is rejected due to excessive overlap
        }

        Ok(resolved_elements)
    }

    /// Assess overlap between element and existing mask
    fn assess_mask_overlap(&self, element: &MaskedElement, mask: &[bool]) -> OverlapResolution {
        let total_length = element.end_pos - element.start_pos;
        let mut overlap_length = 0;
        let mut conflict_regions = Vec::new();
        let mut current_conflict_start = None;

        // Count overlapping positions
        for pos in element.start_pos..element.end_pos {
            if pos < mask.len() && mask[pos] {
                overlap_length += 1;
                
                // Track conflict regions
                if current_conflict_start.is_none() {
                    current_conflict_start = Some(pos);
                }
            } else {
                // End of conflict region
                if let Some(conflict_start) = current_conflict_start {
                    conflict_regions.push((conflict_start, pos));
                    current_conflict_start = None;
                }
            }
        }
        
        // Close final conflict region if needed
        if let Some(conflict_start) = current_conflict_start {
            conflict_regions.push((conflict_start, element.end_pos));
        }

        let overlap_fraction = overlap_length as f64 / total_length as f64;
        
        // Decision logic: keep element if overlap is not excessive
        let keep_current = overlap_fraction < 0.5; // Allow up to 50% overlap
        
        // Quality penalty based on overlap amount
        let quality_adjustment = overlap_fraction * self.min_overlap_penalty;

        OverlapResolution {
            keep_current,
            conflict_regions,
            quality_adjustment,
        }
    }

    /// Check if element should be masked (tandem repeat filtering)
    fn should_mask_element(&self, element: &MaskedElement, mask: &[bool]) -> bool {
        // Check for tandem repeat proximity (C++ tandem distance filtering)
        if !self.check_tandem_distance(element, mask) {
            return false;
        }

        // Minimum length requirement
        let element_length = element.end_pos - element.start_pos;
        if element_length < 30 { // Minimum meaningful element size
            return false;
        }

        // Quality threshold
        if element.alignment_score < 0.0 { // Negative scores indicate poor alignment
            return false;
        }

        true
    }

    /// Check tandem repeat distance constraint
    fn check_tandem_distance(&self, element: &MaskedElement, mask: &[bool]) -> bool {
        let search_start = element.start_pos.saturating_sub(self.tandem_dist);
        let search_end = min(element.end_pos + self.tandem_dist, mask.len());

        // Look for nearby masked regions (potential tandem repeats)
        for pos in search_start..element.start_pos {
            if pos < mask.len() && mask[pos] {
                return false; // Too close to existing repeat
            }
        }
        
        for pos in element.end_pos..search_end {
            if pos < mask.len() && mask[pos] {
                return false; // Too close to existing repeat
            }
        }

        true
    }

    /// Apply element mask to the global mask array
    fn apply_element_mask(&self, element: &MaskedElement, mask: &mut [bool]) -> Result<(), QualityError> {
        if element.start_pos >= mask.len() || element.end_pos > mask.len() {
            return Err(QualityError::MaskingError("Element boundaries exceed mask length".to_string()));
        }

        // Apply mask with buffer
        let mask_start = element.start_pos.saturating_sub(self.mask_buffer);
        let mask_end = min(element.end_pos + self.mask_buffer, mask.len());

        for pos in mask_start..mask_end {
            mask[pos] = true;
        }

        Ok(())
    }

    /// Check if region is heavily masked (>80%)
    fn is_region_heavily_masked(&self, start: usize, end: usize, mask: &[bool]) -> bool {
        if start >= mask.len() || end > mask.len() || start >= end {
            return true;
        }

        let masked_count = (start..end).filter(|&pos| mask[pos]).count();
        let total_length = end - start;
        
        masked_count as f64 / total_length as f64 > 0.8
    }

    /// Validate masking consistency after family processing
    fn validate_masking_consistency(&self, elements: &[MaskedElement], mask: &[bool]) -> Result<(), QualityError> {
        // Check that all elements are properly masked
        for element in elements {
            let mut properly_masked = true;
            
            for pos in element.start_pos..element.end_pos {
                if pos < mask.len() && !mask[pos] {
                    properly_masked = false;
                    break;
                }
            }
            
            if !properly_masked {
                return Err(QualityError::MaskingError(
                    format!("Element {}:{} not properly masked", element.start_pos, element.end_pos)
                ));
            }
        }

        Ok(())
    }

    /// Get masking statistics for analysis
    pub fn get_masking_stats(&self, elements: &[MaskedElement], mask: &[bool]) -> MaskingStatistics {
        let total_elements = elements.len();
        let total_masked_bases = mask.iter().filter(|&&is_masked| is_masked).count();
        let mask_density = total_masked_bases as f64 / mask.len() as f64;
        
        let avg_element_length = if total_elements > 0 {
            elements.iter().map(|e| e.end_pos - e.start_pos).sum::<usize>() as f64 / total_elements as f64
        } else {
            0.0
        };

        let avg_alignment_score = if total_elements > 0 {
            elements.iter().map(|e| e.alignment_score).sum::<f64>() / total_elements as f64
        } else {
            0.0
        };

        MaskingStatistics {
            total_elements,
            total_masked_bases,
            mask_density,
            avg_element_length,
            avg_alignment_score,
        }
    }
}

/// Alignment matrix for boundary refinement
struct AlignmentMatrix {
    matrix: Vec<Vec<i32>>,
    best_end: (usize, usize),
    consensus_len: usize,
    element_len: usize,
}

/// Masking statistics for analysis
#[derive(Debug, Clone)]
pub struct MaskingStatistics {
    pub total_elements: usize,
    pub total_masked_bases: usize,
    pub mask_density: f64,
    pub avg_element_length: f64,
    pub avg_alignment_score: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_progressive_masker_creation() {
        let masker = ProgressiveMasker::default();
        assert_eq!(masker.tandem_dist, 500);
        assert_eq!(masker.mask_buffer, 10);
    }

    #[test]
    fn test_tandem_distance_check() {
        let masker = ProgressiveMasker::new(100);
        let mut mask = vec![false; 1000];
        
        // Create element
        let element = MaskedElement {
            start_pos: 500,
            end_pos: 550,
            alignment_score: 10.0,
            is_reverse: false,
            sequence_id: 0,
            refined_boundaries: true,
        };
        
        // Should pass with clean mask
        assert!(masker.check_tandem_distance(&element, &mask));
        
        // Add nearby masked region
        mask[450] = true;
        assert!(!masker.check_tandem_distance(&element, &mask));
    }

    #[test]
    fn test_overlap_assessment() {
        let masker = ProgressiveMasker::default();
        let mut mask = vec![false; 1000];
        
        // Create element
        let element = MaskedElement {
            start_pos: 100,
            end_pos: 200,
            alignment_score: 10.0,
            is_reverse: false,
            sequence_id: 0,
            refined_boundaries: true,
        };
        
        // No overlap - should keep
        let resolution = masker.assess_mask_overlap(&element, &mask);
        assert!(resolution.keep_current);
        
        // Heavy overlap - should reject
        for i in 100..180 {
            mask[i] = true;
        }
        let resolution = masker.assess_mask_overlap(&element, &mask);
        assert!(!resolution.keep_current);
    }

    #[test]
    fn test_boundary_adjustment() {
        let masker = ProgressiveMasker::default();
        let mut mask = vec![false; 1000];
        
        // Mask some positions
        mask[100] = true;
        mask[101] = true;
        
        // Adjust start boundary - should move right
        let adjusted_start = masker.adjust_boundary_for_mask(100, &mask, true);
        assert!(adjusted_start > 101);
        
        // Adjust end boundary - should move left  
        let adjusted_end = masker.adjust_boundary_for_mask(102, &mask, false);
        assert!(adjusted_end <= 100);
    }
}