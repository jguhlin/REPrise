# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

REPrise is a de novo interspersed repeat detection tool using inexact seeding. This repository contains:

1. **Original C++ implementation** (`REPrise.cpp`, `REPrise.hpp`) - The reference implementation. *Make no changes to this code.*
2. **Rust port** (`REPrise-rs/`) - A 1:1 parity port of the C++ implementation with additional CLI interface
3. **Evaluation and testing tools** - JSON-based equivalence testing between C++ and Rust versions

## Build Commands

### C++ Version
```bash
# Build the main REPrise executable
make

# Build the JSON equivalence test utility  
make eq-cpp

# Clean build artifacts
make clean
```

### Rust Version
```bash
# Run the CLI with default test file
cargo run --manifest-path REPrise-rs/Cargo.toml -- --fasta test/tst.fa

# Run with output to file
cargo run --manifest-path REPrise-rs/Cargo.toml -- --fasta test/tst.fa --out results.tsv

# Run tests
cargo test --manifest-path REPrise-rs/Cargo.toml

# Run integration equivalence tests
cargo test --manifest-path REPrise-rs/Cargo.toml --test equiv -- --nocapture
```

## Testing and Equivalence Validation

The project uses a three-step equivalence testing approach:

1. **Build C++ JSON emitter**: `make eq-cpp`
2. **Run Rust integration harness**: `cargo test --test equiv -- --nocapture`  
3. **Compare outputs**: `pixi run python eval/compare_equiv.py`

Expected result: "EQUIVALENCE PASS" when schemas and outputs match.

## Architecture

### Key Concepts
- **Inexact Seeding**: Uses k-mers with allowed edit distance for initial repeat detection
- **Extension Alignment**: Extends seed matches using banded dynamic programming
- **Masking Strategy**: Prevents overlapping repeat calls through progressive masking

### Core Data Structures
- **Suffix Array**: Used for efficient k-mer occurrence lookup
- **Cache Table**: Pre-computed k-mer variants for inexact matching
- **Mask Array**: Boolean array tracking masked genome regions

### Critical Function Categories

1. **Index Building**:
   - `store_cache`: Builds k-mer cache for inexact matching
   - `build_sortedkmers`: Creates frequency-sorted k-mer priority queue
   - `build_sequence`: Loads FASTA and creates padded numeric sequence

2. **Repeat Detection**:
   - `findkmer`: Locates k-mer occurrences using cache and suffix array
   - `find_bestseed`: Selects highest frequency unmasked seed
   - `masking_align`/`mask_extention_score`: Performs banded DP extension

3. **Masking and Filtering**:
   - `removetandem`: Filters closely spaced occurrences
   - `removemasked`: Excludes already-masked regions
   - `maskbyseed`/`maskbyrepeat`: Updates mask array

## Development Workflow

### Function Porting (C++ → Rust)
1. Check specs in `function_spec.md`
2. Implement in `REPrise-rs/src/alg/repeat.rs`
3. Add unit tests in same file under `cfg(test)`
4. Update status in `conversion_status.md` and `docs/rep_conversion_status.md`
5. Test equivalence with full pipeline

### Key Files to Monitor
- **Status tracking**: `conversion_status.md`, `docs/rep_conversion_status.md`
- **Function specs**: `function_spec.md`
- **Core algorithms**: `REPrise-rs/src/alg/repeat.rs`
- **Integration tests**: `REPrise-rs/tests/equiv.rs`
- **CLI interface**: `REPrise-rs/src/main.rs`

## Testing Strategy

- **Unit tests**: Co-located with implementation in Rust modules
- **Integration tests**: JSON-based equivalence between C++ and Rust
- **CLI tests**: Use `test/tst.fa` as canonical test input
- **Deterministic output**: All operations maintain consistent ordering for reproducible results

The project prioritizes 1:1 behavioral parity between C++ reference and Rust port, with comprehensive testing to ensure equivalence at both function and system levels.
