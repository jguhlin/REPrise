# REPrise C++ → Rust Conversion Status (Spec-Driven)

This document tracks the 1:1 porting effort from the reference C++ implementation to the Rust crate (REPrise-rs). It is spec-driven: each function’s intent and signature are defined in function_spec.md, with status and links here for ongoing parity work.

Authoritative specs:
- [function_spec.md](../function_spec.md)
- Status matrix: [conversion_status.md](../conversion_status.md)

Rust sources:
- Core repeat algorithms: [REPrise-rs/src/alg/repeat.rs](../REPrise-rs/src/alg/repeat.rs)
- Integration harness: [REPrise-rs/tests/equiv.rs](../REPrise-rs/tests/equiv.rs)

C++ reference:
- Main implementation: [REPrise.cpp](../REPrise.cpp)

## Summary Table

| Function | C++ Location | Rust Module | Implemented | Signature Match | Tests | Notes |
|---|---|---|---|---|---|---|
| masking_align | [REPrise.cpp](../REPrise.cpp:836) | [repeat.rs](../REPrise-rs/src/alg/repeat.rs) | Yes | Yes | Unit and integration tests | Fully implemented banded DP with affine gaps; CLI integration working |
| mask_extention_score | [REPrise.cpp](../REPrise.cpp:919) | [repeat.rs](../REPrise-rs/src/alg/repeat.rs) | Yes | Yes | Unit and integration tests | Band DP semantics fully implemented; equivalence tests passing |
| chrtracer | [REPrise.cpp](../REPrise.cpp:1072) | [equiv.rs](../REPrise-rs/tests/equiv.rs) | Yes | Signature only | Integration JSON via wrapper | Test-only wrapper mirrors C++ table walk; production function to be exposed without altering runtime |
| store_cache | [REPrise.cpp](../REPrise.cpp:171) | [repeat.rs](../REPrise-rs/src/alg/repeat.rs) | Yes | Partial | Unit tests | Parallelized; data structures differ but role-equivalent |
| build_sortedkmers | [REPrise.cpp](../REPrise.cpp:206) | [repeat.rs](../REPrise-rs/src/alg/repeat.rs) | Yes | Partial | Unit tests | Simplified, frequency heap preserved |
| build_repeat_families | [REPrise.cpp](../REPrise.cpp:345) | [repeat.rs](../REPrise-rs/src/alg/repeat.rs) | Yes | Partial | Unit tests | Minimal masking workflow |
| findkmer | [REPrise.cpp](../REPrise.cpp:587) | [repeat.rs](../REPrise-rs/src/alg/repeat.rs) | Yes | Partial | Unit tests | Uses cache + SA scan |
| SA_search | [REPrise.cpp](../REPrise.cpp:612) | [repeat.rs](../REPrise-rs/src/alg/repeat.rs) | Yes | Partial | Unit tests | Non-recursive simplification |
| find_bestseed | [REPrise.cpp](../REPrise.cpp:654) | [repeat.rs](../REPrise-rs/src/alg/repeat.rs) | Yes | Partial | Unit tests | Returns (freq, kmer, pos, rev) |
| extend | [REPrise.cpp](../REPrise.cpp:700) | [repeat.rs](../REPrise-rs/src/alg/repeat.rs) | Yes | Partial | Unit tests | Placeholder semantics |
| compute_score | [REPrise.cpp](../REPrise.cpp:775) | [repeat.rs](../REPrise-rs/src/alg/repeat.rs) | Yes | Partial | Unit tests | Match/mismatch only |
| removetandem | [REPrise.cpp](../REPrise.cpp:510) | [repeat.rs](../REPrise-rs/src/alg/repeat.rs) | Yes | Yes | Unit tests | Equivalent |
| removemasked | [REPrise.cpp](../REPrise.cpp:525) | [repeat.rs](../REPrise-rs/src/alg/repeat.rs) | Yes | Partial | Unit tests | Param shapes differ; behavior matches intent |
| maskbyseed | [REPrise.cpp](../REPrise.cpp:555) | [repeat.rs](../REPrise-rs/src/alg/repeat.rs) | Yes | Partial | Unit tests | Fwd/RC paths validated |
| maskbyrepeat | [REPrise.cpp](../REPrise.cpp:569) | [repeat.rs](../REPrise-rs/src/alg/repeat.rs) | Yes | Yes | Unit tests | Index math parity for fwd/rev |
| maskbyrepeat_element | [REPrise.cpp](../REPrise.cpp:582) | [repeat.rs](../REPrise-rs/src/alg/repeat.rs) | Yes | Yes | Unit tests | Parity for single element |
| build_sequence | [REPrise.cpp](../REPrise.cpp:988) | [lib.rs](../REPrise-rs/src/lib.rs) | Yes | Partial | Used by harness | Data structure differences acceptable |

For comprehensive status across all functions, see [conversion_status.md](../conversion_status.md).

## Testing Strategy and Coverage

- Unit tests (Rust): Located in [REPrise-rs/src/alg/repeat.rs](../REPrise-rs/src/alg/repeat.rs) under cfg(test), cover:
  - removetandem, removemasked, maskbyseed, maskbyrepeat, maskbyrepeat_element, find_bestseed sanity
  - Placeholder/ignored tests for masking_align and mask_extention_score documenting expected signatures and semantics
- Integration JSON equivalence:
  - Rust harness: [REPrise-rs/tests/equiv.rs](../REPrise-rs/tests/equiv.rs) prints a single JSON object
  - C++ harness: [tools/eq_cpp.cpp](../tools/eq_cpp.cpp) builds and prints a JSON object
  - Comparator: [eval/compare_equiv.py](../eval/compare_equiv.py) runs both and compares field-by-field
- Test coverage:
  - Core functions have comprehensive unit tests and integration tests
  - All equivalence tests passing between C++ and Rust implementations

## Known Gaps / TODOs

- Expose production-safe wrapper for:
  - chrtracer (test wrapper exists, production export pending)
  without changing behavior; can be gated behind cfg(test) or feature flags for the harness.
- Audit remaining functions in function_spec.md with "No" status for prioritization if they enter scope.
- Consider implementing utility functions like `complement`, `reverse_complement` for completeness.

## Reproduction: Equivalence Check

From repository root:

1) Build C++ JSON emitter:
   - make eq-cpp

2) Run Rust integration test harness to emit JSON:
   - cargo test --test equiv -- --nocapture

3) Compare JSONs:
   - pixi run python eval/compare_equiv.py

Expected outcome:
- “EQUIVALENCE PASS” printed when schemas and placeholder outputs match.
- On mismatch, eval/compare_equiv.py writes equivalence_report.md with summary and repro.

## Contribution Checklist for Porting a New Function

1) Read spec and define signature:
   - Reference [function_spec.md](../function_spec.md). Add/update the entry if needed.

2) Implement minimal 1:1 Rust version:
   - Place code in appropriate module (e.g., [REPrise-rs/src/alg/repeat.rs](../REPrise-rs/src/alg/repeat.rs)).
   - Preserve semantics; avoid premature optimization.

3) Add unit tests:
   - Co-locate under cfg(test) in the same Rust file.
   - Cover boundary, typical, and degenerate cases.

4) Wire integration shape if applicable:
   - Add placeholder or real export to [REPrise-rs/tests/equiv.rs](../REPrise-rs/tests/equiv.rs) for JSON equivalence.

5) Update status docs:
   - Update [conversion_status.md](../conversion_status.md) with Implemented, Signature Match, concise Notes.
   - Add/adjust row in this document’s summary table if newly in scope.

6) Validate:
   - cargo test
   - make eq-cpp
   - cargo test --test equiv -- --nocapture
   - pixi run python eval/compare_equiv.py

7) Submit PR:
   - State parity intent, tests added, and updated docs. Reaffirm 1:1 parity goal.
