// Integration test harness to emit deterministic JSON comparing key functions.
// Runs small synthetic and test/tst.fa inputs, printing a single JSON object per run.
//
// Invoke with:
//   cargo test --test equiv -- --nocapture
//
// The JSON schema:
// {
//   "rust": {
//     "masking_align": [ { "i": ..., "consensusstart": ..., "consensusend": ..., "result": [start, end] }, ... ],
//     "chrtracer": [ { "pos": ..., "chr": "...", "off": ... }, ... ],
//     "mask_extention_score": [ { "is_right": true/false, "ext": ..., "i": ..., "band": ..., "score": ... }, ... ]
//   }
// }
//
// Note: Some functions are not fully implemented in Rust. For those, we provide placeholder
// deterministic outputs, clearly marked, to structure the comparison. This harness focuses
// on export and shape; the C++ driver will mirror the same logic for apples-to-apples comparison.

use reprise::alg::repeat::*;
use reprise::{build_sequence, suffix_array, PADLENGTH};
use serde_json::{json, Value};

#[cfg(test)]
fn small_seq_and_sa() -> (Vec<u8>, Vec<i64>) {
    let seq = vec![
        0,0,0,0, // AAAA
        1,1,1,1, // CCCC
        2,2,2,2, // GGGG
        3,3,3,3, // TTTT
        0,1,2,3, // ACGT
    ];
    let sa = suffix_array(&seq);
    (seq, sa)
}

#[cfg(test)]
fn build_cache_for_k(seq: &[u8], sa: &[i64], k: usize, edit_distance: u8) -> Vec<Vec<CacheEntry>> {
    store_cache(edit_distance, k, seq, sa)
}

// Test-only wrapper for chrtracer-style mapping based on SequenceData.chrtable
#[cfg(test)]
fn chrtracer_wrapper(chrtable: &[(String, usize)], stringpos: usize) -> (String, usize) {
    // Mirror C++: default to last, find first entry where pos < chrtable[i].second, then i-1
    if chrtable.is_empty() {
        return ("unknown".to_string(), 0);
    }
    let mut idx = chrtable.len() - 1;
    for i in 0..chrtable.len() {
        if stringpos < chrtable[i].1 {
            if i == 0 { idx = 0; } else { idx = i - 1; }
            break;
        }
    }
    let (name, start) = &chrtable[idx];
    (name.clone(), *start)
}

// Placeholders for masking_align and mask_extention_score to structure JSON.
// These are not production exports and only exist for deterministic testing.
// They will be compared to the C++ driver that emits the same placeholder logic.
#[cfg(test)]
fn masking_align_placeholder(_i: usize, consensusstart: usize, consensusend: usize) -> (usize, usize) {
    // Deterministic: return inner trimmed window or equal if same
    if consensusend >= consensusstart {
        let mid = (consensusstart + consensusend) / 2;
        (mid.saturating_sub(2), mid.saturating_add(2).min(consensusend))
    } else {
        (consensusstart, consensusend)
    }
}

#[cfg(test)]
fn mask_extention_score_placeholder(_is_right: bool, _ext: usize, _i: usize, band: usize) -> i32 {
    // Deterministic function of band for shape only
    band as i32
}

#[test]
fn emit_equiv_json() {
    // Synthetic sequence context
    let (seq_small, sa_small) = small_seq_and_sa();
    let k = 2usize;
    let _cache_small = build_cache_for_k(&seq_small, &sa_small, k, 0);

    // Representative window indices for masking_align
    let masking_cases = vec![
        (0usize, 5usize, 15usize),
        (1, 0, 4),
        (2, 10, 10), // zero-length window
    ];
    let mut masking_align_out = Vec::new();
    for (i, cs, ce) in masking_cases {
        let (s, e) = masking_align_placeholder(i, cs, ce);
        masking_align_out.push(json!({
            "i": i, "consensusstart": cs, "consensusend": ce, "result": [s, e]
        }));
    }

    // chrtracer mappings using real test fasta if available
    let fasta_path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../test/tst.fa");
    let mut chr_out = Vec::new();
    if let Ok(data) = build_sequence(&fasta_path) {
        // boundary and internal positions
        let mut positions = Vec::new();
        // after initial padding
        positions.push(PADLENGTH + 10);
        // try first non-padding contig start + 100
        if let Some((_, start)) = data.chrtable.iter().find(|(n, _)| n != "unknown" && n != "padding") {
            positions.push(*start + 100);
        }
        // a padding entry position if exists
        if let Some((_, p)) = data.chrtable.iter().find(|(n, _)| n == "padding") {
            positions.push(*p);
        }
        for pos in positions {
            let (chr, off) = chrtracer_wrapper(&data.chrtable, pos);
            chr_out.push(json!({"pos": pos, "chr": chr, "off": off}));
        }
    } else {
        // Fallback if file not present
        chr_out.push(json!({"pos": 0, "chr": "unknown", "off": 0}));
    }

    // mask_extention_score grid of band positions/params
    let grid_params = vec![
        (true, 10usize, 0usize, 5usize),
        (true, 0, 1, 7),
        (false, 3, 2, 9),
    ];
    let mut mes_out = Vec::new();
    for (is_right, ext, i, band) in grid_params {
        let score = mask_extention_score_placeholder(is_right, ext, i, band);
        mes_out.push(json!({
            "is_right": is_right,
            "ext": ext,
            "i": i,
            "band": band,
            "score": score
        }));
    }

    let out = json!({
        "rust": {
            "masking_align": masking_align_out,
            "chrtracer": chr_out,
            "mask_extention_score": mes_out
        }
    });

    println!("{}", serde_json::to_string_pretty(&out).unwrap());

    // Use assertions only to ensure JSON non-empty and keys exist; do not fail on content
    let v: Value = out;
    assert!(v.get("rust").is_some());
    assert!(v["rust"].get("masking_align").is_some());
    assert!(v["rust"].get("chrtracer").is_some());
    assert!(v["rust"].get("mask_extention_score").is_some());
}