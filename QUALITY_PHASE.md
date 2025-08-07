# QUALITY_PHASE: Comprehensive Plan for C++ Parity in Rust REPrise

## Executive Summary

The Rust REPrise implementation successfully identifies repeat patterns with 4.3x better performance than C++, but currently finds 1000+ unfiltered families instead of C++'s 45 high-quality families. This document outlines a comprehensive plan to implement the sophisticated quality filtering system that transforms raw k-mer matches into biologically meaningful repeat families.

## Current State Analysis

### What Works ✅
- **Core Algorithm**: K-mer detection, suffix arrays, frequency analysis (99%+ correlation with C++)
- **Performance**: 4.3x faster than C++ (6.2s vs 26.3s)
- **Seed Detection**: Identifies same high-frequency patterns as C++
- **Basic Infrastructure**: FASTA loading, chromosome parsing, basic masking

### Critical Gap ❌
- **Quality Filtering**: Accepts all k-mers ≥3 frequency without validation
- **Extension Alignment**: Missing banded dynamic programming for consensus building
- **Family Validation**: No biological significance testing
- **Progressive Masking**: Incomplete masking strategy allows overlap

## C++ Quality System Analysis

### Multi-Stage Filtering Pipeline
1. **Entropy Filtering**: `compute_entropy(query) > -0.7` (eliminates low-complexity)
2. **Frequency Threshold**: `MINFREQ = 3` (minimum occurrences)  
3. **Tandem Removal**: `TANDEMDIST = 500` (prevents tandem repeat overlap)
4. **Masking Awareness**: Excludes already-processed regions
5. **Extension Quality**: Banded DP with improvement requirements
6. **Length Validation**: `MINLENGTH = 50bp` consensus minimum
7. **Progressive Masking**: Prevents family overlap

### Termination Logic
- Processes families until quality stops improving
- Early termination after stagnant score improvement
- Results in ~45 families due to diminishing returns, not arbitrary limits

### Key Functions Requiring Rust Implementation
- **`extend()`**: Banded dynamic programming extension
- **`masking_align()`**: Individual element boundary refinement  
- **`build_repeat_families()`**: Main coordination logic
- **Quality scoring**: Multi-criteria family validation

## Implementation Architecture

### Phase 1: Core Extension Alignment (Weeks 1-2)

#### 1.1 Banded Dynamic Programming
```rust
// New module: src/quality/extension.rs
pub struct BandedAligner {
    match_score: i32,        // Default: 1
    mismatch_score: i32,     // Default: -1  
    gap_open_score: i32,     // Default: -5
    gap_extend_score: i32,   // Default: -1
    cap_penalty: i32,        // Default: -20
    band_width: usize,       // Default: 11 (2*OFFSETWIDTH + 1)
}

impl BandedAligner {
    pub fn extend_bidirectional(
        &self, 
        seed: &[u8], 
        positions: &[usize], 
        sequence: &[u8]
    ) -> ConsensusResult {
        // Implement left + right extension with quality termination
    }
}
```

#### 1.2 Consensus Building
```rust
// New module: src/quality/consensus.rs  
pub struct ConsensusBuilder {
    min_improvement: usize,  // Default: 3
    max_extend: usize,       // Default: 10000
    stop_threshold: usize,   // Default: 100
}

impl ConsensusBuilder {
    pub fn build_consensus(
        &self,
        alignments: &[AlignmentResult]
    ) -> Option<ConsensusSequence> {
        // Multi-way alignment with base selection by max score
    }
}
```

### Phase 2: Quality Filtering Pipeline (Weeks 2-3)

#### 2.1 Entropy and Complexity Filtering
```rust
// New module: src/quality/filtering.rs
pub struct QualityFilter {
    max_entropy: f64,        // Default: -0.7
    min_length: usize,       // Default: 50  
    min_freq: usize,         // Default: 3
    gc_content_range: (f64, f64), // Default: (0.2, 0.8)
}

impl QualityFilter {
    pub fn prefilter_kmer(&self, kmer: &[u8]) -> bool {
        // Implement entropy calculation matching C++
        let entropy = compute_entropy(kmer);
        entropy > self.max_entropy
    }
    
    pub fn validate_family(&self, family: &RepeatFamily) -> bool {
        // Multi-criteria validation: length, frequency, GC content
    }
}
```

#### 2.2 Progressive Quality Assessment
```rust
// New module: src/quality/validation.rs
pub struct FamilyValidator {
    quality_threshold: f64,   // Minimum family quality score
    improvement_window: usize, // Track recent quality trends
}

pub struct QualityMetrics {
    pub length_score: f64,
    pub frequency_score: f64,  
    pub consensus_quality: f64,
    pub coverage_score: f64,
    pub composite_score: f64,
}
```

### Phase 3: Progressive Masking System (Week 3-4)

#### 3.1 Advanced Masking Coordination
```rust
// New module: src/quality/masking.rs
pub struct ProgressiveMasker {
    tandem_dist: usize,      // Default: 500
    mask_buffer: usize,      // Buffer around masked regions
}

impl ProgressiveMasker {
    pub fn mask_by_family(&mut self, 
        family: &ValidatedFamily, 
        mask: &mut [bool]
    ) {
        // Implement sophisticated boundary detection
        // Handle overlapping families with quality priority
    }
    
    pub fn refine_boundaries(&self,
        element: &RepeatElement,
        consensus: &[u8],
        sequence: &[u8]
    ) -> (usize, usize) {
        // Implement masking_align equivalent
    }
}
```

### Phase 4: Integration and Optimization (Week 4-5)

#### 4.1 Main Pipeline Integration
```rust
// Enhanced src/bin/reprise_cpp_port.rs
pub struct QualityPipeline {
    filter: QualityFilter,
    aligner: BandedAligner,  
    consensus_builder: ConsensusBuilder,
    validator: FamilyValidator,
    masker: ProgressiveMasker,
}

impl QualityPipeline {
    pub fn process_families(&mut self, 
        kmers: &mut BinaryHeap<KmerEntry>,
        sequence: &[u8],
        mask: &mut [bool]
    ) -> Vec<ValidatedFamily> {
        let mut families = Vec::new();
        let mut quality_tracker = QualityTracker::new();
        
        while let Some(candidate) = kmers.pop() {
            // 1. Entropy and frequency filtering
            if !self.filter.prefilter_kmer(&candidate.kmer) {
                continue;
            }
            
            // 2. Extension alignment
            let extension_result = self.aligner.extend_bidirectional(
                &candidate.kmer, &candidate.positions, sequence
            )?;
            
            // 3. Consensus building  
            let consensus = self.consensus_builder.build_consensus(
                &extension_result.alignments
            )?;
            
            // 4. Family validation
            let quality = self.validator.assess_quality(&consensus);
            if !self.validator.meets_threshold(&quality) {
                continue;
            }
            
            // 5. Quality trend analysis
            quality_tracker.update(quality.composite_score);
            if quality_tracker.should_terminate() {
                break; // Early termination like C++
            }
            
            // 6. Progressive masking
            self.masker.mask_by_family(&consensus, mask);
            
            families.push(ValidatedFamily {
                consensus,
                quality,
                family_id: families.len(),
            });
        }
        
        families
    }
}
```

## Expected Outcomes

### Family Count Convergence
- **Current**: 1000+ unfiltered k-mers
- **Phase 1**: ~500 families (with basic extension)
- **Phase 2**: ~150 families (with quality filtering)  
- **Phase 3**: ~75 families (with progressive masking)
- **Phase 4**: ~45 families (with C++-equivalent termination)

### Quality Improvements
- **Entropy Filtering**: Eliminates low-complexity false positives
- **Extension Validation**: Ensures biological significance
- **Progressive Masking**: Prevents overlapping family calls
- **Early Termination**: Stops when quality plateaus

### Performance Characteristics
- **Maintained Speed**: Quality filtering happens after fast k-mer detection
- **Memory Efficient**: Reuse alignment matrices, early candidate elimination
- **Scalable**: Architecture supports larger genomes with same quality

## Implementation Timeline

### Week 1-2: Core Extension Alignment
- [ ] Implement `BandedAligner` with C++-equivalent scoring
- [ ] Add bidirectional extension with quality termination
- [ ] Implement `ConsensusBuilder` for multi-way alignment
- [ ] Unit tests for alignment accuracy vs C++

### Week 2-3: Quality Filtering Pipeline  
- [ ] Implement entropy calculation matching C++ `compute_entropy`
- [ ] Add multi-criteria family validation
- [ ] Implement quality scoring and trend analysis
- [ ] Integration tests for filtering accuracy

### Week 3-4: Progressive Masking System
- [ ] Implement `ProgressiveMasker` with boundary refinement
- [ ] Add overlap resolution with quality-based priority
- [ ] Implement `masking_align` equivalent for precise boundaries
- [ ] Validate masking prevents family overlap

### Week 4-5: Integration and Validation
- [ ] Integrate all components into main pipeline
- [ ] Add early termination logic matching C++ behavior
- [ ] Performance optimization and memory efficiency
- [ ] Comprehensive testing against C++ equivalence

## Success Criteria

### Functional Parity
- [ ] **Family Count**: ~45 families (±10% of C++ output)
- [ ] **Quality Correlation**: Top 10 families match C++ frequencies within 5%
- [ ] **Length Distribution**: Family lengths match C++ distribution
- [ ] **No Overlaps**: Families are non-overlapping like C++

### Performance Maintenance  
- [ ] **Speed**: Maintain <10s runtime (vs C++ 26s)
- [ ] **Memory**: Keep memory usage reasonable (<1GB)
- [ ] **Scalability**: Architecture supports larger genomes

### Quality Validation
- [ ] **Biological Significance**: Families represent real repetitive elements
- [ ] **Consensus Quality**: Generated sequences have high quality scores
- [ ] **Reproducibility**: Consistent results across runs

## Risk Mitigation

### Technical Risks
- **Complexity**: Phase implementation allows incremental validation
- **Performance**: Early filtering maintains speed advantages  
- **Memory**: Matrix reuse and bounded processing prevent bloat

### Integration Risks  
- **Breaking Changes**: Maintain backward compatibility with existing interfaces
- **Testing**: Comprehensive unit tests for each component
- **Validation**: Continuous comparison against C++ reference

## Long-term Impact

This quality phase transforms the Rust REPrise from a fast k-mer detector into a complete, production-ready repeat detection system that:

- **Matches C++ Quality**: Same biological accuracy and family count
- **Maintains Performance Advantage**: 4x+ speed improvement retained
- **Enables Production Use**: High-quality families suitable for downstream analysis
- **Supports Scalability**: Architecture handles larger genomes efficiently

The completed system will provide the genomics community with a high-performance, high-quality alternative to the C++ implementation while maintaining full functional parity.

---

*Document Version: 1.0*  
*Date: August 7, 2025*  
*Objective: Transform Rust REPrise from 1000+ unfiltered families to C++-equivalent 45 high-quality families*