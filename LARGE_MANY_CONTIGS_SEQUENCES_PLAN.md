# REPrise Scaling Implementation Guide (Definitive Production Blueprint v2.2)

**Version:** 2.2  
**Status:** Final, Ready for Implementation  
**Contact:** Project Lead

## 1. Executive Summary

This document is the **definitive production blueprint** for re-architecting the REPrise repeat-finding tool. It is a complete and actionable guide for the implementation engineer, incorporating multiple rounds of rigorous review to ensure safety, scalability, and correctness.

The primary objective is to overcome the severe limitations of the original C++ implementation, enabling REPrise to process massive and highly fragmented genomes (e.g., 12 Gbp, >250,000 contigs) efficiently. This will be achieved via a **dual-strategy architecture**:

*   **Path A (Fidelity): Suffix Array (SA) Based.** A high-fidelity port of the original logic. Resource-intensive but serves as a validation gold standard.
*   **Path B (Speed & Scale): SA-Free Heuristics.** The **recommended default path**, prioritizing speed and accessibility on commodity hardware.

This blueprint is built on hardened architectural principles and now includes concrete implementation skeletons for the highest-risk components, along with specific recommendations for high-performance libraries to accelerate development and maximize performance.

## 2. Recommended Crates & Libraries

To ensure a high-quality, performant, and maintainable implementation, the use of the following best-in-class crates is strongly recommended.

| Crate | Purpose | Why it's Recommended |
| :--- | :--- | :--- |
| **`needletail`** | FASTA/Q Parsing | Extremely fast, memory-efficient, and robust parser for biological sequences. |
| **`block-aligner`** | Gapped Alignment | SIMD-accelerated (SSE/AVX2) Smith-Waterman alignment. Orders of magnitude faster than a scalar implementation. |
| **`ahash`** | Hashing | A very fast, DoS-resistant hashing algorithm. Use with a fixed seed for deterministic `HashMap`s. |
| **`rayon`** | Data Parallelism | The de-facto standard for easy and powerful data parallelism in Rust. |
| **`crossbeam`** | Concurrency Primitives | Provides high-performance, lock-free channels and other utilities essential for our producer-consumer pipeline. |
| **`loom`** | Concurrency Testing | A tool for exhaustively testing concurrent code to find subtle race conditions. **Essential for validating the `Bitmask`**. |
| **`clap`** | CLI Parsing | The standard for building ergonomic, feature-rich, and self-documenting command-line interfaces. |
| **`indicatif`** | Progress Bars | Provides clean, configurable progress bars for long-running operations. |
| **`thiserror`** | Error Handling | Boilerplate-free derivation of the `std::error::Error` trait for our custom error enum. |
| **`mimalloc`** | (Optional) Allocator | A drop-in, high-performance memory allocator that can provide a global speed boost in highly concurrent applications. |

## 3. High-Risk Implementation Skeletons (Priority 1)

The following components must be implemented and stress-tested first, as the entire project's safety and performance depend on them.

### 3.1. The `Bitmask` and `ClaimGuard`

This is the cornerstone of safe parallelism. The `Bitmask` provides a memory-efficient atomic bitmap, and the `ClaimGuard` ensures panic-safe, automatic release of claims.

```rust
// In src/mask.rs
use std::sync::atomic::{AtomicU64, Ordering};
use std::ops::Range;

// A memory-efficient, thread-safe bitmask for tracking claimed genomic regions.
pub struct Bitmask {
    data: Vec<AtomicU64>,
    len_bases: u64,
}

impl Bitmask {
    const BITS_PER_WORD: u64 = 64;

    pub fn new(len_bases: u64) -> Self {
        let num_words = (len_bases + Self::BITS_PER_WORD - 1) / Self::BITS_PER_WORD;
        let data = (0..num_words).map(|_| AtomicU64::new(0)).collect();
        Self { data, len_bases }
    }

    /// Attempts to atomically claim all bits in a range.
    /// This is an all-or-nothing operation with internal rollback.
    pub fn claim_range(&self, range: &Range<u64>) -> bool {
        if range.start >= range.end || range.end > self.len_bases { return false; }

        let mut claimed_words = Vec::new();
        for word_idx in (range.start / Self::BITS_PER_WORD)..=((range.end - 1) / Self::BITS_PER_WORD) {
            let word_start = word_idx * Self::BITS_PER_WORD;
            let word_end = word_start + Self::BITS_PER_WORD;

            let claim_start = range.start.max(word_start);
            let claim_end = range.end.min(word_end);

            let mut bitmask_for_word = 0u64;
            for i in claim_start..claim_end {
                bitmask_for_word |= 1u64 << (i % Self::BITS_PER_WORD);
            }

            let atomic_word = &self.data[word_idx as usize];
            let mut current_val = atomic_word.load(Ordering::Acquire);
            loop {
                // Check for conflict inside the loop to avoid race conditions
                if (current_val & bitmask_for_word) != 0 {
                    // Conflict detected. Roll back all previous claims.
                    for (idx, mask) in &claimed_words {
                        self.data[*idx].fetch_and(!*mask, Ordering::Release);
                    }
                    return false;
                }
                
                match atomic_word.compare_exchange_weak(
                    current_val, current_val | bitmask_for_word,
                    Ordering::AcqRel, Ordering::Acquire,
                ) {
                    Ok(_) => {
                        claimed_words.push((word_idx as usize, bitmask_for_word));
                        break;
                    }
                    Err(new_current_val) => current_val = new_current_val,
                }
            }
        }
        true
    }

    /// Releases all bits in a range. Idempotent.
    pub fn release_range(&self, range: &Range<u64>) {
        if range.start >= range.end || range.end > self.len_bases { return; }
        
        for word_idx in (range.start / Self::BITS_PER_WORD)..=((range.end - 1) / Self::BITS_PER_WORD) {
            let word_start = word_idx * Self::BITS_PER_WORD;
            let claim_start = range.start.max(word_start);
            let claim_end = range.end.min(word_start + Self::BITS_PER_WORD);
            
            let mut bitmask = 0u64;
            for i in claim_start..claim_end {
                bitmask |= 1u64 << (i % Self::BITS_PER_WORD);
            }
            
            self.data[word_idx as usize].fetch_and(!bitmask, Ordering::Release);
        }
    }

    /// A cheap, non-atomic check for pre-filtering.
    pub fn is_range_free(&self, range: &Range<u64>) -> bool {
        if range.start >= range.end || range.end > self.len_bases { return false; }
        
        for word_idx in (range.start / Self::BITS_PER_WORD)..=((range.end - 1) / Self::BITS_PER_WORD) {
            let word_start = word_idx * Self::BITS_PER_WORD;
            let claim_start = range.start.max(word_start);
            let claim_end = range.end.min(word_start + Self::BITS_PER_WORD);
            
            let mut bitmask = 0u64;
            for i in claim_start..claim_end {
                bitmask |= 1u64 << (i % Self::BITS_PER_WORD);
            }
            
            if (self.data[word_idx as usize].load(Ordering::Relaxed) & bitmask) != 0 {
                return false;
            }
        }
        true
    }
}

/// A panic-safe RAII guard for atomically claiming a region.
#[must_use = "Claim must be held by a variable to ensure its lifetime"]
pub struct ClaimGuard<'a> {
    mask: &'a Bitmask,
    claimed_range: Option<Range<u64>>,
}

impl<'a> Drop for ClaimGuard<'a> {
    fn drop(&mut self) {
        if let Some(range) = self.claimed_range.take() {
            self.mask.release_range(&range);
        }
    }
}

impl<'a> ClaimGuard<'a> {
    pub fn new(mask: &'a Bitmask, range: Range<u64>) -> Option<Self> {
        if mask.claim_range(&range) {
            Some(Self { mask, claimed_range: Some(range) })
        } else {
            None
        }
    }
}
```

### 3.2. Streaming Candidate Pair Pipeline

This producer-consumer model is essential to prevent OOM errors on large genomes.

```rust
// In src/pipeline.rs
use crossbeam::channel;
use rayon::scope;
use std::sync::Arc;

const CHANNEL_CAPACITY: usize = 1_024_000; // Bounded channel to control memory.

pub fn run_pipeline(num_threads: usize, genome: Arc<Genome>, mask: Arc<Bitmask>) {
    let (tx, rx) = channel::bounded::<(u64, u64)>(CHANNEL_CAPACITY);

    scope(|s| {
        // Producer thread
        s.spawn(move |_| {
            // `generate_pairs` reads from the k-mer index or SA.
            // It pre-filters pairs using `mask.is_range_free` before sending.
            generate_pairs(tx, &genome, &mask);
        });

        // Worker threads
        for _ in 0..num_threads {
            let rx_clone = rx.clone();
            let genome_clone = Arc::clone(&genome);
            let mask_clone = Arc::clone(&mask);
            s.spawn(move |_| {
                worker(rx_clone, &genome_clone, &mask_clone);
            });
        }
    });
}
```

## 4. Core Architecture & Implementation Details

### 4.1. The `Genome` API

The binary **requires a 64-bit target**. A static assertion will enforce this.

```rust
// In src/main.rs or lib.rs
#[cfg(target_pointer_width = "32")]
compile_error!("REPrise requires a 64-bit target to handle large genomes.");

// In src/genome.rs
use std::ops::Range;

pub type ContigId = u32;

pub struct Genome { /* ... */ }

impl Genome {
    /// Builds a Genome from a FASTA file using the `needletail` crate.
    pub fn from_fasta(path: &Path) -> Result<Self, anyhow::Error>;
    pub fn slice(&self, range: Range<u64>) -> &[u8];
    pub fn contig_of(&self, pos: u64) -> Option<ContigId>;
    pub fn contig_range(&self, id: ContigId) -> Option<Range<u64>>;
    pub fn is_within_one_contig(&self, start: u64, len: usize) -> bool;
}
```

### 4.2. K-mer Engine & Alignment

*   **K-mer Indexing:** Use a **contig-by-contig** pass. The indexer must use the canonical k-mer logic (handling reverse-complements and skipping ambiguous bases).
*   **Alignment:** The core `bounded_gapped_alignment` function should be implemented using the **`block-aligner`** crate for SIMD-accelerated performance.

### 4.3. Error Handling

A custom error enum will provide clear, specific error types.

```rust
// In src/error.rs
use thiserror::Error;

#[derive(Error, Debug)]
pub enum RepriseError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Invalid k-mer containing non-ACGT base at position {0}")]
    InvalidKmer(u64),
    #[error("Alignment attempted across contig boundary")]
    BoundaryViolation,
    #[error("External tool '{0}' failed: {1}")]
    ExternalTool(String, String),
    // ... other error variants
}
```

## 5. Implementation Roadmap & Timeline (Next 2 Weeks)

| Day | Milestone |
| :-- | :--- |
| **1-2** | Implement `Bitmask` (claim/release) and its unit tests. |
| **3-4** | Implement `ClaimGuard` and a panic-safety test using `std::panic::catch_unwind`. |
| **5-6** | Write the `loom` stress test for the mask and integrate it into CI. |
| **7** | Build the bounded channel producer-consumer stub; verify memory usage stays low on a 1 Gbp mock genome. |
| **8-9** | Implement the canonical k-mer engine (using `ahash` with a fixed seed) and a simple indexer for k ≤ 32. |
| **10** | Wire the k-mer indexer into the producer to generate real candidate pairs. |
| **11-12**| Define the custom `RepriseError` enum and refactor early functions to return `Result`. |
| **13-14**| Merge documentation files, bump version, and run the first full integration test on a 100 Mbp dataset. |

## 6. Testing, Documentation, & Project Hygiene

*   **Testing:** A `loom` stress test for the `Bitmask` is non-negotiable. A `cargo bench` suite with small, medium, and large datasets is required. A nightly CI job on a ≥500 Mbp synthetic genome will catch scaling regressions.
*   **Documentation:**
    *   **`README.md`:** Must include a "Hardware Requirements" table and a "Recommended Crates" section.
    *   **`OUTPUT_FORMATS.md`:** Must detail all output file columns and coordinate systems (0-based, half-open).
    *   **`CONTRIBUTING.md`:** Must define code style and PR process.
    *   **`LICENSES.md`:** Must list licenses of all third-party dependencies and clarify that external tools are not bundled.

## 7. Command-Line Interface (CLI) Design

```rust
use clap::{Parser, ValueEnum};

#[derive(ValueEnum, Clone, Debug, Copy)]
pub enum Strategy {
    /// Suffix Array path: High fidelity, high resource usage (requires >2x genome size on SSD).
    #[value(name = "sa")]
    SuffixArray,
    /// Heuristic path: Faster, lower resource usage. Recommended for all new analyses.
    #[value(name = "heuristic")]
    Heuristic,
}

#[derive(Parser, Debug)]
#[command(author, version, about = "REPrise: A tool for de-novo repeat identification in large genomes.")]
pub struct Args {
    /// Input FASTA file path.
    #[arg(short, long, value_name = "FILE")]
    pub input: std::path::PathBuf,

    /// Output file prefix. The program will write <PREFIX>.bed, <PREFIX>.masked.fa, etc.
    #[arg(short, long, value_name = "PREFIX")]
    pub output: String,

    /// The seed-pairing strategy to use.
    #[arg(long, value_enum, default_value_t = Strategy::Heuristic)]
    pub strategy: Strategy,
    
    /// Number of parallel threads to use. [default: 0, meaning all available cores]
    #[arg(long, short = 'p', value_name = "INT", default_value_t = 0)]
    pub threads: usize,

    // ... other REPrise parameters like -k, -minfreq, etc. will be added here ...
}
```

