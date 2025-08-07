# GEMINI.md

This file provides guidance to Gemini when working with code in this repository.

## Project Overview

REPrise is a de novo interspersed repeat detection tool using inexact seeding. This repository contains:

1.  **Original C++ implementation** (`REPrise.cpp`, `REPrise.hpp`) - The reference implementation. *Make no changes to this code.*
2.  **Rust port** (`REPrise-rs/`) - A port of the C++ implementation that is being evolved for large-scale genomics.
3.  **Evaluation and testing tools** - JSON-based equivalence testing between C++ and Rust versions.

## Next Generation Architecture (REPrise v2.0)

To support extremely large (>10 Gbp) and highly fragmented (>100,000 contigs) genomes, REPrise is moving beyond 1:1 C++ parity towards a new, highly scalable architecture. The focus is on performance, memory efficiency, and robust handling of complex genomic data.

### Core Principles

*   **Memory Bounded:** All operations will be designed to work within a configurable memory footprint, avoiding out-of-memory errors on large datasets.
*   **Streaming First:** Genomes will be processed as streams of contigs or chunks, never requiring the full sequence to be loaded into memory.
*   **Concurrency by Design:** The pipeline will be built on a producer-consumer model to maximize parallelism and throughput.
*   **Contig-Awareness:** All algorithms will respect contig boundaries to ensure biological accuracy.

### Key Architectural Changes

1.  **Dual-Strategy Repeat Finding:**
    *   **`heuristic` (Default):** A new, SA-free approach using k-mer hashing and indexing for maximum scalability. This will be the recommended strategy for large or fragmented genomes.
    *   **`suffix-array` (Legacy):** The existing high-fidelity port of the C++ implementation will be maintained for smaller genomes and validation purposes.

2.  **Producer-Consumer Pipeline:**
    *   A `crossbeam` channel will be used to stream candidate repeat pairs from producers (k-mer indexers) to consumers (alignment workers).
    *   The channel will be bounded, providing backpressure to control memory usage.
    *   Work will be batched to reduce channel overhead.

3.  **Atomic Region Masking:**
    *   A custom `Bitmask` with atomic operations (`AtomicU64`) will be used to track claimed genomic regions.
    *   `ClaimGuard` (RAII) will ensure that claimed regions are automatically and safely released, even in the event of a panic. This is critical for safe concurrency.

4.  **Contig-Aware Genome Representation:**
    *   The `Genome` API will manage global vs. local coordinates, allowing algorithms to operate on a unified position space while respecting contig boundaries.
    *   The `needletail` crate will be used for high-performance, streaming FASTA parsing.

### Recommended Crates for v2.0

*   **`needletail`**: For fast, memory-efficient FASTA/Q parsing.
*   **`ahash`**: For a fast, DoS-resistant hashing algorithm for k-mer indexing.
*   **`rayon`**: For data parallelism in tasks like initial k-mer counting.
*   **`crossbeam`**: For high-performance, lock-free channels for the main pipeline.
*   **`clap`**: For a robust and user-friendly CLI.
*   **`thiserror`**: For clean, structured error handling.

## Development Workflow

### v2.0 Development
1.  Implement the `heuristic` strategy using the producer-consumer pipeline.
2.  Develop robust k-mer indexing and candidate generation for the producer.
3.  Implement the alignment and scoring logic in the consumer.
4.  Ensure the `Bitmask` and `ClaimGuard` are used for safe concurrent processing.
5.  Add comprehensive benchmarks for large and fragmented genomes.

### Legacy Maintenance
*   The C++ implementation and the `suffix-array` Rust strategy are considered feature-complete.
*   Changes to the legacy code should only be for bug fixes or to maintain compatibility with the testing harness.

## Testing Strategy

*   **Unit tests**: For individual components of the new architecture.
*   **Integration tests**: For the full `heuristic` pipeline.
*   **Equivalence testing**: The existing JSON-based tests will be used to ensure the `suffix-array` strategy remains in parity with the C++ version.
*   **Large-scale benchmarks**: New tests will be added to validate the performance and scalability of the `heuristic` strategy on large datasets.