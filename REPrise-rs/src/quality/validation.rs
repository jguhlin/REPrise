//! Family validation and quality scoring
//! 
//! Implements multi-criteria validation for repeat families and quality trend analysis
//! to determine when processing should terminate (matching C++ behavior).

use crate::quality::{QualityError, ConsensusSequence, ExtensionResult};
use std::collections::VecDeque;

/// Family validator with multi-criteria quality assessment
pub struct FamilyValidator {
    pub quality_threshold: f64,    // Minimum composite quality score
    pub min_length: usize,         // Minimum consensus length (C++ MINLENGTH)
    pub min_coverage: usize,       // Minimum element coverage
    pub max_gap_fraction: f64,     // Maximum gap content
    pub gc_content_range: (f64, f64), // Acceptable GC content range
}

impl Default for FamilyValidator {
    fn default() -> Self {
        Self {
            quality_threshold: 0.6,     // Require 60% quality
            min_length: 50,             // Match C++ MINLENGTH
            min_coverage: 3,            // Minimum 3 elements
            max_gap_fraction: 0.3,      // Max 30% gaps
            gc_content_range: (0.15, 0.85), // Broader than filter (biological range)
        }
    }
}

impl FamilyValidator {
    pub fn new(quality_threshold: f64) -> Self {
        Self {
            quality_threshold,
            ..Default::default()
        }
    }

    /// Assess comprehensive quality of a consensus sequence
    pub fn assess_quality(
        &self,
        consensus: &ConsensusSequence,
        extension_result: &ExtensionResult,
    ) -> Result<QualityMetrics, QualityError> {
        // Length quality component
        let length_score = self.calculate_length_score(consensus.length);
        
        // Frequency/coverage quality component
        let avg_coverage = consensus.coverage.iter().sum::<usize>() as f64 / consensus.coverage.len() as f64;
        let frequency_score = self.calculate_frequency_score(avg_coverage);
        
        // Consensus quality component (position-wise agreement)
        let consensus_score = consensus.overall_quality;
        
        // Coverage consistency component
        let coverage_score = self.calculate_coverage_consistency(&consensus.coverage);
        
        // GC content component  
        let gc_score = self.calculate_gc_score(consensus);
        
        // Gap content component
        let gap_score = self.calculate_gap_score(consensus);
        
        // Extension quality component
        let extension_score = extension_result.extension_quality;
        
        // Alignment score component
        let alignment_score = if extension_result.alignments.is_empty() {
            0.0
        } else {
            let total_score: f64 = extension_result.alignments.iter().map(|a| a.alignment_score).sum();
            let avg_score = total_score / extension_result.alignments.len() as f64;
            (avg_score / 100.0).min(1.0).max(0.0) // Normalize alignment scores
        };
        
        // Composite quality score (weighted combination)
        let composite_score = 0.25 * length_score +
                             0.20 * frequency_score +
                             0.20 * consensus_score +
                             0.10 * coverage_score +
                             0.10 * gc_score +
                             0.05 * gap_score +
                             0.05 * extension_score +
                             0.05 * alignment_score;
        
        Ok(QualityMetrics {
            length_score,
            frequency_score,
            consensus_quality: consensus_score,
            coverage_score,
            gc_score,
            gap_score,
            extension_score,
            alignment_score,
            composite_score,
        })
    }

    /// Check if family meets acceptance threshold
    pub fn meets_threshold(&self, quality: &QualityMetrics) -> bool {
        quality.composite_score >= self.quality_threshold
    }

    /// Validate consensus meets basic requirements
    pub fn validate_basic_requirements(&self, consensus: &ConsensusSequence) -> bool {
        // Length requirement
        if consensus.length < self.min_length {
            return false;
        }

        // Coverage requirement
        let avg_coverage = consensus.coverage.iter().sum::<usize>() as f64 / consensus.coverage.len() as f64;
        if avg_coverage < self.min_coverage as f64 {
            return false;
        }

        // Gap content requirement
        let gap_fraction = consensus.gap_positions.iter().filter(|&&is_gap| is_gap).count() as f64 
            / consensus.gap_positions.len() as f64;
        if gap_fraction > self.max_gap_fraction {
            return false;
        }

        // GC content requirement
        let gc_content = self.calculate_gc_content_fraction(consensus);
        if gc_content < self.gc_content_range.0 || gc_content > self.gc_content_range.1 {
            return false;
        }

        true
    }

    /// Calculate length quality score (longer is generally better, with diminishing returns)
    fn calculate_length_score(&self, length: usize) -> f64 {
        if length < self.min_length {
            return 0.0;
        }
        
        // Sigmoid-like function: rapid increase up to ~200bp, then diminishing returns
        let normalized_length = (length as f64 - self.min_length as f64) / 150.0;
        (1.0 - (-normalized_length).exp()).min(1.0)
    }

    /// Calculate frequency/coverage quality score
    fn calculate_frequency_score(&self, avg_coverage: f64) -> f64 {
        if avg_coverage < self.min_coverage as f64 {
            return 0.0;
        }
        
        // Log scale for frequency: high frequency indicates strong repeat signal
        let log_coverage = (avg_coverage.max(1.0)).ln();
        (log_coverage / 5.0).min(1.0) // Normalize to reasonable range
    }

    /// Calculate coverage consistency across positions
    fn calculate_coverage_consistency(&self, coverage: &[usize]) -> f64 {
        if coverage.is_empty() {
            return 0.0;
        }

        let mean: f64 = coverage.iter().sum::<usize>() as f64 / coverage.len() as f64;
        
        if mean == 0.0 {
            return 0.0;
        }

        // Calculate coefficient of variation (std dev / mean)
        let variance: f64 = coverage.iter()
            .map(|&x| {
                let diff = x as f64 - mean;
                diff * diff
            })
            .sum::<f64>() / coverage.len() as f64;
        
        let std_dev = variance.sqrt();
        let cv = std_dev / mean;
        
        // Convert CV to score (lower variation = higher score)
        (1.0 - cv.min(1.0)).max(0.0)
    }

    /// Calculate GC content quality score
    fn calculate_gc_score(&self, consensus: &ConsensusSequence) -> f64 {
        let gc_content = self.calculate_gc_content_fraction(consensus);
        
        // Penalty for extreme GC content (optimal around 40-60%)
        let optimal_gc = 0.5;
        let deviation = (gc_content - optimal_gc).abs();
        let penalty = deviation * 2.0; // Linear penalty
        
        (1.0 - penalty).max(0.0)
    }

    /// Calculate gap content quality score
    fn calculate_gap_score(&self, consensus: &ConsensusSequence) -> f64 {
        if consensus.gap_positions.is_empty() {
            return 1.0;
        }

        let gap_fraction = consensus.gap_positions.iter().filter(|&&is_gap| is_gap).count() as f64 
            / consensus.gap_positions.len() as f64;
        
        // Linear penalty for gap content
        (1.0 - (gap_fraction / self.max_gap_fraction)).max(0.0)
    }

    /// Calculate GC content as fraction
    fn calculate_gc_content_fraction(&self, consensus: &ConsensusSequence) -> f64 {
        let mut gc_count = 0;
        let mut total_count = 0;

        for &base in &consensus.sequence {
            if base < 4 {
                if base == 1 || base == 2 { // C or G
                    gc_count += 1;
                }
                total_count += 1;
            }
        }

        if total_count == 0 {
            return 0.5; // Neutral GC content
        }

        gc_count as f64 / total_count as f64
    }

    /// Get detailed validation report
    pub fn get_validation_report(&self, quality: &QualityMetrics, consensus: &ConsensusSequence) -> ValidationReport {
        let basic_requirements_met = self.validate_basic_requirements(consensus);
        let meets_threshold = self.meets_threshold(quality);
        
        let issues = self.identify_quality_issues(quality, consensus);
        
        ValidationReport {
            overall_pass: basic_requirements_met && meets_threshold,
            basic_requirements_met,
            meets_quality_threshold: meets_threshold,
            composite_score: quality.composite_score,
            quality_threshold: self.quality_threshold,
            identified_issues: issues,
            recommendations: self.generate_recommendations(quality, consensus),
        }
    }

    /// Identify specific quality issues
    fn identify_quality_issues(&self, quality: &QualityMetrics, consensus: &ConsensusSequence) -> Vec<String> {
        let mut issues = Vec::new();

        if consensus.length < self.min_length {
            issues.push(format!("Consensus too short: {} < {}", consensus.length, self.min_length));
        }

        if quality.frequency_score < 0.3 {
            issues.push("Low element frequency/coverage".to_string());
        }

        if quality.consensus_quality < 0.5 {
            issues.push("Poor position-wise consensus agreement".to_string());
        }

        if quality.coverage_score < 0.5 {
            issues.push("Inconsistent coverage across positions".to_string());
        }

        let gc_content = self.calculate_gc_content_fraction(consensus);
        if gc_content < 0.2 || gc_content > 0.8 {
            issues.push(format!("Extreme GC content: {:.2}", gc_content));
        }

        let gap_fraction = consensus.gap_positions.iter().filter(|&&is_gap| is_gap).count() as f64 
            / consensus.gap_positions.len() as f64;
        if gap_fraction > self.max_gap_fraction {
            issues.push(format!("Excessive gap content: {:.2}", gap_fraction));
        }

        issues
    }

    /// Generate improvement recommendations
    fn generate_recommendations(&self, quality: &QualityMetrics, consensus: &ConsensusSequence) -> Vec<String> {
        let mut recommendations = Vec::new();

        if quality.length_score < 0.5 {
            recommendations.push("Consider adjusting extension parameters for longer consensus".to_string());
        }

        if quality.frequency_score < 0.5 {
            recommendations.push("Increase minimum frequency threshold to improve signal".to_string());
        }

        if quality.consensus_quality < 0.5 {
            recommendations.push("Improve consensus building parameters for better agreement".to_string());
        }

        if quality.coverage_score < 0.5 {
            recommendations.push("Check for uneven sequence coverage or alignment issues".to_string());
        }

        recommendations
    }
}

/// Comprehensive quality metrics for repeat families
#[derive(Debug, Clone)]
pub struct QualityMetrics {
    pub length_score: f64,
    pub frequency_score: f64,
    pub consensus_quality: f64,
    pub coverage_score: f64,
    pub gc_score: f64,
    pub gap_score: f64,
    pub extension_score: f64,
    pub alignment_score: f64,
    pub composite_score: f64,
}

/// Validation report with detailed analysis
#[derive(Debug)]
pub struct ValidationReport {
    pub overall_pass: bool,
    pub basic_requirements_met: bool,
    pub meets_quality_threshold: bool,
    pub composite_score: f64,
    pub quality_threshold: f64,
    pub identified_issues: Vec<String>,
    pub recommendations: Vec<String>,
}

/// Quality trend tracker for early termination
pub struct QualityTracker {
    recent_scores: VecDeque<f64>,
    window_size: usize,
    stagnation_threshold: usize,
    min_improvement_rate: f64,
    families_processed: usize,
}

impl QualityTracker {
    pub fn new() -> Self {
        Self {
            recent_scores: VecDeque::new(),
            window_size: 10,              // Track last 10 families
            stagnation_threshold: 10,     // Stop after 10 families without improvement
            min_improvement_rate: 0.01,   // Minimum improvement rate
            families_processed: 0,
        }
    }

    /// Update with new family quality score
    pub fn update(&mut self, quality_score: f64) {
        self.recent_scores.push_back(quality_score);
        self.families_processed += 1;

        // Keep window size limited
        while self.recent_scores.len() > self.window_size {
            self.recent_scores.pop_front();
        }
    }

    /// Check if processing should terminate due to quality plateau
    pub fn should_terminate(&self) -> bool {
        // Need minimum families before considering termination
        if self.families_processed < 20 {
            return false;
        }

        // Need full window for trend analysis
        if self.recent_scores.len() < self.window_size {
            return false;
        }

        // Check for improvement trend
        let recent_avg = self.recent_scores.iter().skip(self.window_size / 2).sum::<f64>() 
            / (self.window_size / 2) as f64;
        
        let early_avg = self.recent_scores.iter().take(self.window_size / 2).sum::<f64>() 
            / (self.window_size / 2) as f64;

        let improvement_rate = (recent_avg - early_avg) / early_avg.max(0.01);

        // Terminate if no significant improvement
        improvement_rate < self.min_improvement_rate
    }

    /// Get current processing statistics
    pub fn get_statistics(&self) -> QualityStatistics {
        let current_avg_quality = if self.recent_scores.is_empty() {
            0.0
        } else {
            self.recent_scores.iter().sum::<f64>() / self.recent_scores.len() as f64
        };

        let quality_trend = if self.recent_scores.len() >= 2 {
            let first_half = self.recent_scores.iter().take(self.recent_scores.len() / 2).sum::<f64>() 
                / (self.recent_scores.len() / 2) as f64;
            let second_half = self.recent_scores.iter().skip(self.recent_scores.len() / 2).sum::<f64>() 
                / (self.recent_scores.len() - self.recent_scores.len() / 2) as f64;
            second_half - first_half
        } else {
            0.0
        };

        QualityStatistics {
            families_processed: self.families_processed,
            current_avg_quality,
            quality_trend,
            should_terminate: self.should_terminate(),
        }
    }

    pub fn families_processed(&self) -> usize {
        self.families_processed
    }
}

/// Quality processing statistics
#[derive(Debug, Clone)]
pub struct QualityStatistics {
    pub families_processed: usize,
    pub current_avg_quality: f64,
    pub quality_trend: f64,
    pub should_terminate: bool,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_family_validator_creation() {
        let validator = FamilyValidator::default();
        assert_eq!(validator.quality_threshold, 0.6);
        assert_eq!(validator.min_length, 50);
    }

    #[test]
    fn test_length_score_calculation() {
        let validator = FamilyValidator::default();
        
        // Too short
        assert_eq!(validator.calculate_length_score(30), 0.0);
        
        // Minimum length
        assert!(validator.calculate_length_score(50) > 0.0);
        
        // Good length
        assert!(validator.calculate_length_score(100) > validator.calculate_length_score(50));
    }

    #[test]
    fn test_quality_tracker() {
        let mut tracker = QualityTracker::new();
        
        // Add some quality scores
        for i in 0..15 {
            tracker.update(0.5 + (i as f64) * 0.01); // Improving quality
        }
        
        assert!(!tracker.should_terminate()); // Should continue with improving trend
        
        // Add stagnant scores
        for _ in 0..15 {
            tracker.update(0.6); // Flat quality
        }
        
        assert!(tracker.should_terminate()); // Should terminate due to stagnation
    }

    #[test]
    fn test_basic_validation() {
        let validator = FamilyValidator::default();
        
        // Create a valid consensus
        let consensus = ConsensusSequence {
            sequence: vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3; 10], // 200 bp mixed
            left_extension: vec![],
            seed_sequence: vec![0; 60],
            right_extension: vec![],
            position_qualities: vec![0.8; 200],
            overall_quality: 0.8,
            coverage: vec![5; 200],
            length: 200,
            gap_positions: vec![false; 200],
        };
        
        assert!(validator.validate_basic_requirements(&consensus));
    }
}