# REPrise Function Conversion Status

| Function | Implemented | Signature Match | Notes |
|----------|-------------|-----------------|-------|
| store_cache | Yes | Partial | Rust implementation uses different approach with parallel processing and different data structures |
| build_sortedkmers | Yes | Yes | Ported tie-break and counting semantics. Deterministic order: freq desc, then lexicographic k-mer asc. Rust: [REPrise-rs/src/alg/repeat.rs](REPrise-rs/src/alg/repeat.rs:164) C++: [REPrise.cpp](REPrise.cpp:206) |
| build_repeat_families | Yes | Partial | Rust implementation is simplified |
| findkmer | Yes | Yes | Parity on forward enumeration using cache begin/end; invalid bases filtered; deterministic ascending genomic positions. Rust: [REPrise-rs/src/alg/repeat.rs](REPrise-rs/src/alg/repeat.rs:145) C++: [REPrise.cpp](REPrise.cpp:587) |
| SA_search | Yes | Partial | Rust implementation is simplified version of recursive C++ version |
| find_bestseed | Yes | Partial | Returns frequency, seed, and computed pos/rev locally; behavior matches C++ selection logic |
| extend | Yes | Partial | Rust implementation is a placeholder returning zero extension |
| compute_score | Yes | Partial | Rust implementation is simplified match/mismatch scoring |
| removetandem | Yes | Yes | Functionally equivalent |
| removemasked | Yes | Yes | C++ parity for forward and RC masking windows. Rust: [REPrise-rs/src/alg/repeat.rs](REPrise-rs/src/alg/repeat.rs:97) C++: [REPrise.cpp](REPrise.cpp:525) |
| maskbyseed | Yes | Partial | Rust implementation has slightly different parameters |
| maskbyrepeat | Yes | Yes | Matches C++ masking index math for fwd/rev |
| maskbyrepeat_element | Yes | Yes | Matches C++ element masking behavior |
| masking_align | Yes | Yes | Fully implemented banded DP with affine gaps; integrated into CLI and passing equivalence tests. Rust: [REPrise-rs/src/alg/repeat.rs](REPrise-rs/src/alg/repeat.rs:406) C++: [REPrise.cpp](REPrise.cpp:836) |
| mask_extention_score | Yes | Yes | Fully implemented with band semantics and affine gaps; passing equivalence tests. Rust: [REPrise-rs/src/alg/repeat.rs](REPrise-rs/src/alg/repeat.rs:293) C++: [REPrise.cpp](REPrise.cpp:919) |
| build_sequence | Yes | Partial | Rust implementation follows similar logic but uses different data structures |
| chrtracer | Yes | Signature only | Provided via test-only wrapper; production export pending. C++ [REPrise.cpp](REPrise.cpp:1072); Rust harness wrapper in [REPrise-rs/tests/equiv.rs](REPrise-rs/tests/equiv.rs). |
| allocate_space | No | N/A | Not implemented in Rust |
| freespace | No | N/A | Not implemented in Rust |
| num_to_char | No | N/A | Not implemented in Rust |
| char_to_num | Yes | Partial | Rust implementation is private function with similar functionality |
| complement | No | N/A | Not implemented in Rust |
| reverse_complement | No | N/A | Not implemented in Rust |
| compute_entropy | No | N/A | Not implemented in Rust |
| default_k | No | N/A | Not implemented in Rust |
| display_time | No | N/A | Not implemented in Rust |
| print_usage | No | N/A | Not implemented in Rust |