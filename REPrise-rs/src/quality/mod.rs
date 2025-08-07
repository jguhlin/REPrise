//! Quality filtering pipeline for REPrise repeat family detection
//! 
//! This module implements the sophisticated quality filtering system that transforms
//! raw k-mer matches into biologically meaningful repeat families, matching C++ REPrise behavior.

pub mod extension;
pub mod consensus;
pub mod filtering;
pub mod validation;
pub mod masking;

pub use extension::*;
pub use consensus::*;
pub use filtering::*;
pub use validation::*;
pub use masking::*;

use std::collections::BinaryHeap;

/// Main quality pipeline coordinating all filtering stages
pub struct QualityPipeline {
    pub filter: QualityFilter,
    pub aligner: BandedAligner,
    pub consensus_builder: ConsensusBuilder,
    pub validator: FamilyValidator,
    pub masker: ProgressiveMasker,
    pub quality_tracker: QualityTracker,
}

impl QualityPipeline {
    pub fn new() -> Self {
        Self {
            filter: QualityFilter::default(),
            aligner: BandedAligner::default(),
            consensus_builder: ConsensusBuilder::default(),
            validator: FamilyValidator::default(),
            masker: ProgressiveMasker::default(),
            quality_tracker: QualityTracker::new(),
        }
    }

    pub fn with_params(
        max_entropy: f64,
        min_length: usize,
        min_freq: usize,
        quality_threshold: f64
    ) -> Self {
        Self {
            filter: QualityFilter::new(max_entropy, min_length, min_freq),
            aligner: BandedAligner::default(),
            consensus_builder: ConsensusBuilder::default(),
            validator: FamilyValidator::new(quality_threshold),
            masker: ProgressiveMasker::default(),
            quality_tracker: QualityTracker::new(),
        }
    }

    /// Process a candidate k-mer through the complete quality pipeline
    pub fn process_candidate_kmer(
        &mut self,
        kmer: &[u8],
        positions: &[usize],
        rev_flags: &[bool],
        sequence: &[u8],
        mask: &mut [bool],
        quality_params: &QualityParams,
    ) -> Result<Option<ValidatedFamily>, QualityError> {
        // Stage 1: Pre-filtering (entropy, frequency, basic checks)
        if !self.filter.prefilter_kmer(kmer, positions.len()) {
            return Ok(None);
        }

        // Stage 2: Extension alignment (banded DP)
        let extension_result = self.aligner.extend_bidirectional(
            kmer, positions, rev_flags, sequence, mask
        )?;

        // Stage 3: Consensus building
        let consensus = self.consensus_builder.build_consensus(&extension_result)?;

        // Stage 4: Family validation
        let quality = self.validator.assess_quality(&consensus, &extension_result)?;
        if !self.validator.meets_threshold(&quality) {
            return Ok(None);
        }

        // Stage 5: Quality tracking and early termination
        self.quality_tracker.update(quality.composite_score);
        
        // Stage 6: Progressive masking
        let masked_elements = self.masker.mask_by_family(&consensus, &extension_result, mask)?;

        Ok(Some(ValidatedFamily {
            consensus,
            quality,
            elements: masked_elements,
            family_id: self.quality_tracker.families_processed(),
        }))
    }

    /// Check if processing should terminate based on quality trends
    pub fn should_terminate(&self) -> bool {
        self.quality_tracker.should_terminate()
    }

    /// Get quality statistics for current processing session
    pub fn get_statistics(&self) -> QualityStatistics {
        self.quality_tracker.get_statistics()
    }
}

/// Parameters for quality assessment
#[derive(Debug, Clone)]
pub struct QualityParams {
    pub max_entropy: f64,
    pub min_length: usize,
    pub min_freq: usize,
    pub quality_threshold: f64,
    pub tandem_dist: usize,
    pub max_families: usize,
}

impl Default for QualityParams {
    fn default() -> Self {
        Self {
            max_entropy: -0.7,     // Match C++ MAXENTROPY
            min_length: 50,        // Match C++ MINLENGTH
            min_freq: 3,           // Match C++ MINFREQ
            quality_threshold: 0.6, // Composite quality threshold
            tandem_dist: 500,      // Match C++ TANDEMDIST
            max_families: 100,     // Conservative maximum
        }
    }
}

/// Validated repeat family output
#[derive(Debug, Clone)]
pub struct ValidatedFamily {
    pub consensus: ConsensusSequence,
    pub quality: QualityMetrics,
    pub elements: Vec<MaskedElement>,
    pub family_id: usize,
}

/// Error types for quality processing
#[derive(Debug)]
pub enum QualityError {
    ExtensionFailed(String),
    ConsensusBuilding(String),
    ValidationError(String),
    MaskingError(String),
    InsufficientData,
}

impl std::fmt::Display for QualityError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            QualityError::ExtensionFailed(msg) => write!(f, "Extension failed: {}", msg),
            QualityError::ConsensusBuilding(msg) => write!(f, "Consensus building failed: {}", msg),
            QualityError::ValidationError(msg) => write!(f, "Validation error: {}", msg),
            QualityError::MaskingError(msg) => write!(f, "Masking error: {}", msg),
            QualityError::InsufficientData => write!(f, "Insufficient data for quality assessment"),
        }
    }
}

impl std::error::Error for QualityError {}