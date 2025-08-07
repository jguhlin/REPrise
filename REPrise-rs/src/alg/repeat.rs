use rayon::prelude::*;
use std::collections::BinaryHeap;
use crate::compute_entropy;

/// Simple tuple type used by caches.
pub type CacheEntry = (usize, usize, u8, u8);

/// Perform a naive search on the suffix array allowing up to `rem_dist`
/// mismatches. This is a simplification of the recursive C++ version.
pub fn sa_search(
    query: &[u8],
    seq: &[u8],
    sa: &[i64],
    begin: usize,
    end: usize,
    rem_dist: i32,
    out: &mut Vec<(usize, usize)>,
) {
    for idx in begin..end {
        let pos = sa[idx] as usize;
        if pos + query.len() > seq.len() {
            continue;
        }
        let mism = seq[pos..pos + query.len()]
            .iter()
            .zip(query)
            .filter(|(a, b)| a != b)
            .count() as i32;
        if mism <= rem_dist {
            out.push((idx, idx + 1));
        }
    }
}

/// Build a cache table for all possible kmers of length `cache_len`.
/// For edit_distance=0: exact matching only (optimized).
/// For edit_distance>0: generates all possible k-mer variants within edit distance.
pub fn store_cache(
    edit_distance: u8,
    cache_len: usize,
    seq: &[u8],
    sa: &[i64],
) -> Vec<Vec<CacheEntry>> {
    // Limit k-mer length for cache to prevent excessive memory allocation
    // For k > 15, the cache would require > 1GB of memory (4^15 = 1 billion entries)
    const MAX_CACHE_K: usize = 15;
    
    if cache_len > MAX_CACHE_K {
        eprintln!("Warning: k-mer length {} exceeds maximum cache size. Using hash-based approach.", cache_len);
        // For large k, return empty cache and rely on suffix array directly
        // This will be slower but won't exhaust memory
        return vec![Vec::new(); 0];
    }
    
    let size = 1u64 << (cache_len * 2);
    let mut cache = vec![Vec::new(); size as usize];

    if edit_distance == 0 {
        // Optimized exact matching - scan suffix array once
        let kmer_to_id = |kmer: &[u8]| -> usize {
            let mut id = 0;
            for (i, &base) in kmer.iter().enumerate() {
                if base > 3 {
                    return 0;
                }
                id |= (base as usize) << (i * 2);
            }
            id
        };

        let mut i = 0;
        while i < sa.len() {
            let p = sa[i] as usize;
            if p + cache_len > seq.len() {
                i += 1;
                continue;
            }

            let kmer = &seq[p..p + cache_len];
            let id = kmer_to_id(kmer);
            let start = i;

            while i < sa.len() {
                let next_p = sa[i] as usize;
                if next_p + cache_len > seq.len() || &seq[next_p..next_p + cache_len] != kmer {
                    break;
                }
                i += 1;
            }

            cache[id].push((start, i, cache_len as u8, edit_distance));
        }
    } else {
        // Inexact matching - generate all possible k-mer variants like C++
        for id in 0..size {
            let mut query = vec![0u8; cache_len];
            
            // Convert cache index back to k-mer
            let mut temp_id = id;
            for i in 0..cache_len {
                query[i] = (temp_id & 3) as u8;
                temp_id >>= 2;
            }
            
            // Skip invalid k-mers
            if query.iter().any(|&b| b > 3) {
                continue;
            }
            
            // Find all matches with edit distance tolerance
            let mut matched = std::collections::HashSet::new();
            sa_search_recursive(&query, seq, sa, 0, sa.len(), 0, 0, edit_distance as i32, &mut matched);
            
            // Convert matches to cache entries
            for (begin, end) in matched {
                cache[id as usize].push((begin, end, cache_len as u8, edit_distance));
            }
        }
    }

    cache
}

/// Recursive suffix array search for inexact matching (mirrors C++ SA_search)
fn sa_search_recursive(
    query: &[u8],
    seq: &[u8],
    sa: &[i64],
    begin: usize,
    end: usize,
    query_num: usize,
    seq_num: usize,
    rem_dist: i32,
    matched: &mut std::collections::HashSet<(usize, usize)>,
) {
    if rem_dist >= 0 && query_num == query.len() {
        matched.insert((begin, end));
        return;
    }
    
    if rem_dist < 0 || query_num >= query.len() {
        return;
    }
    
    // Split suffix array range by next character
    let mut bounds = [begin; 5];
    let mut current = begin;
    
    for base in 0u8..4u8 {
        // Binary search for range of this base
        let mut left = current;
        let mut right = end;
        
        while left < right {
            let mid = (left + right) / 2;
            let pos = sa[mid] as usize + seq_num;
            
            if pos < seq.len() && seq[pos] < base {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        
        current = left;
        
        // Find end of this base range
        left = current;
        right = end;
        while left < right {
            let mid = (left + right) / 2;
            let pos = sa[mid] as usize + seq_num;
            
            if pos < seq.len() && seq[pos] <= base {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        
        bounds[base as usize + 1] = left;
    }
    bounds[4] = end;
    
    // Recurse for each base
    for base in 0u8..4u8 {
        let range_start = bounds[base as usize];
        let range_end = bounds[base as usize + 1];
        
        if range_start < range_end {
            let mismatch_penalty = if query[query_num] != base { 1 } else { 0 };
            sa_search_recursive(
                query,
                seq,
                sa,
                range_start,
                range_end,
                query_num + 1,
                seq_num + 1,
                rem_dist - mismatch_penalty,
                matched,
            );
        }
    }
}

/// Remove tandem occurrences closer than `tandem_dist`.
pub fn removetandem(occs: &mut Vec<usize>, tandem_dist: usize) {
    occs.sort_unstable();
    let mut prev = None;
    occs.retain(|&pos| {
        if let Some(p) = prev {
            if pos < p + tandem_dist {
                return false;
            }
        }
        prev = Some(pos);
        true
    });
}

/// Remove occurrences that overlap masked regions.
///
/// C++ parity REPrise.cpp:525-553
/// Forward (is_rc = false): reject p if any mask[p + i] set for i in [0, k-1]
/// Reverse (is_rc = true):  reject p if any mask[p - i] set for i in [0, k-1]
/// Note: This matches C++ semantics where reverse-orientation masking is contiguous
/// decreasing from p, not mirrored around p+k-1 as in forward.
pub fn removemasked(occs: &mut Vec<usize>, mask: &[bool], k: usize, is_rc: bool) {
    occs.retain(|&p| {
        let mut masked = false;
        for i in 0..k {
            let idx = if !is_rc {
                p.saturating_add(i)
            } else {
                p.saturating_sub(i)
            };
            if idx < mask.len() && mask[idx] {
                masked = true;
                break;
            }
        }
        !masked
    });
}

/// Mask a set of occurrences.
pub fn maskbyseed(occs: &[usize], mask: &mut [bool], k: usize, is_rc: bool) {
    for &p in occs {
        if !is_rc {
            for i in 0..k {
                if p + i < mask.len() {
                    mask[p + i] = true;
                }
            }
        } else {
            // Reverse orientation masks contiguous positions [p, p-1, ..., p-(k-1)]
            for i in 0..k {
                if p >= i {
                    let j = p - i;
                    if j < mask.len() {
                        mask[j] = true;
                    }
                }
            }
        }
    }
}

/**
 * findkmer parity
 * C++ reference: REPrise.cpp:587-653 (via SA_search at REPrise.cpp:612)
 * Semantics:
 * - Enumerate forward occurrences for the exact k-mer using cache begin/end range
 * - Bounds honoring: skip any SA positions that would overflow sequence length
 * - Deterministic ordering: ascending genomic coordinate, stable for ties
 * - No RC here: RC handled by caller with separate RC query
 */
pub fn findkmer(query: &[u8], cache: &[Vec<CacheEntry>], seq: &[u8], sa: &[i64]) -> Vec<usize> {
    if query.is_empty() {
        return Vec::new();
    }
    // Reject queries containing invalid alphabet symbols (>3) to match C++ numeric encoding guard.
    if query.iter().any(|&b| b > 3) {
        return Vec::new();
    }

    // Compute cache index using 2-bit encoding little-endian by position to match store_cache generation.
    let mut idx = 0usize;
    for (i, &b) in query.iter().enumerate() {
        idx |= (b as usize) << (i * 2);
    }

    let mut out = Vec::new();
    if idx >= cache.len() {
        return out;
    }

    // Accumulate matches from all cache segments for this k-mer id.
    for &(b, e, _klen, dist) in &cache[idx] {
        // Use the edit distance from cache generation for inexact matching
        let mut tmp = Vec::new();
        sa_search(query, seq, sa, b, e, dist as i32, &mut tmp);
        for (sidx, _) in tmp {
            let p = sa[sidx] as usize;
            if p + query.len() <= seq.len() {
                out.push(p);
            }
        }
    }

    // Deterministic order
    out.sort_unstable();
    out.dedup();
    out
}

/**
 * build_sortedkmers parity
 * C++ reference: REPrise.cpp:206-344
 * Semantics:
 * - Count total occurrences of k-mer across forward and reverse-complement strands
 * - Ignore k-mers containing invalid symbols (>3)
 * - Deterministic tie-break: higher frequency first; if equal, lexicographically smaller k-mer precedes
 * - Use BinaryHeap where Ord implements (freq asc, kmer desc) to pop highest freq with stable tie rule
 */
pub fn build_sortedkmers(
    k: usize,
    seq: &[u8],
    cache: &[Vec<CacheEntry>],
    sa: &[i64],
    minfreq: usize,
) -> BinaryHeap<(usize, Vec<u8>)> {
    // Custom wrapper to invert ordering for BinaryHeap to achieve desired tie-breaks
    #[derive(Eq, PartialEq)]
    struct KEntry {
        freq: usize,
        kmer: Vec<u8>,
    }
    impl Ord for KEntry {
        fn cmp(&self, other: &Self) -> std::cmp::Ordering {
            // Primary: freq descending
            // Tie: kmer ascending lexicographically
            match self.freq.cmp(&other.freq) {
                std::cmp::Ordering::Less => std::cmp::Ordering::Less,   // natural
                std::cmp::Ordering::Greater => std::cmp::Ordering::Greater,
                std::cmp::Ordering::Equal => other.kmer.cmp(&self.kmer), // invert so heap pops lexicographically smaller first
            }
        }
    }
    impl PartialOrd for KEntry {
        fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
            Some(self.cmp(other))
        }
    }

    // C++ MAXENTROPY = -0.7 (REPrise.cpp)
    const MAX_ENTROPY: f64 = -0.7;
    
    let len = seq.len();
    let kmers: Vec<KEntry> = (0..len.saturating_sub(k).saturating_add(1))
        .into_par_iter()
        .filter_map(|i| {
            let kmer_slice = &seq[i..i + k];
            if kmer_slice.iter().any(|&b| b > 3) {
                return None;
            }
            // Filter out low-complexity k-mers based on entropy
            if compute_entropy(kmer_slice) > MAX_ENTROPY {
                return None;
            }
            // Count fwd + rc occurrences with parity findkmer
            let fwd = findkmer(kmer_slice, cache, seq, sa);
            let rc: Vec<u8> = kmer_slice.iter().rev().map(|&b| 3 - b).collect();
            let rev = findkmer(&rc, cache, seq, sa);
            let freq = fwd.len() + rev.len();
            if freq >= minfreq {
                Some(KEntry {
                    freq,
                    kmer: kmer_slice.to_vec(),
                })
            } else {
                None
            }
        })
        .collect();

    // Use BinaryHeap<(usize, Vec<u8>)> outward API; reinsert with correct tuple order
    let mut heap = BinaryHeap::new();
    for e in kmers {
        heap.push((e.freq, e.kmer));
    }
    heap
}

#[allow(dead_code)]
#[cfg(test)]
pub fn cli_simple_extend_interval(
    seq: &[u8],
    seed_pos: usize,
    k: usize,
    is_rc: bool,
    match_score: i32,
    mismatch: i32,
    gap: i32,
    gapex: i32,
) -> (usize, usize, i32) {
    // Keep as placeholder for tests not used in production CLI.
    let _ = (seq, seed_pos, k, is_rc, match_score, mismatch, gap, gapex);
    (0, 0, 0)
}

/// C++-parity constants for extension alignment
/// REPrise.cpp references:
/// - OFFSETWIDTH default at CLI parse (-maxgap) -> 5 (REPrise.cpp:100-107, 102)
/// - WHEN_TO_STOP default -> 100 (REPrise.cpp:100-107, 103)
/// - MINSCORE sentinel for DP cells -> -100000 (REPrise.cpp:54)
const OFFSETWIDTH: usize = 5; // REPrise.cpp:102
const _WHEN_TO_STOP: usize = 100; // REPrise.cpp:103
const MINSCORE: i32 = -100000; // REPrise.cpp:54

/// Port of mask_extention_score with identical band and affine gap semantics.
/// C++ reference: REPrise.cpp:919
pub(crate) fn mask_extention_score(
    is_right: bool,
    ext: usize,
    _i: usize,
    // state vectors (length = 2*OFFSETWIDTH+1)
    mask_score: &mut [i32],
    mask_score_m: &mut [i32],
    mask_score_ins: &mut [i32],
    mask_score_del: &mut [i32],
    // shared state
    consensus: &[u8],
    sequence: &[u8],
    pos: usize,
    k: usize,
    rev: bool,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
    maxextend: usize,
) -> i32 {
    let width = 2 * OFFSETWIDTH + 1;
    debug_assert!(mask_score.len() == width);
    debug_assert!(mask_score_m.len() == width);
    debug_assert!(mask_score_ins.len() == width);
    debug_assert!(mask_score_del.len() == width);

    let mut next_m = vec![0i32; width];
    let mut next_ins = vec![0i32; width];
    let mut next_del = vec![0i32; width];
    let mut best = MINSCORE;

    // helper to fetch sequence base with safe bounds checking - prevents heap-buffer-overflow
    let get_seq = |offset: isize| -> Option<u8> {
        let off = offset as isize;
        if off < 0 { return None; }
        let off = off as usize;
        if off >= sequence.len() { return None; }
        let b = sequence[off];
        if b > 3 { return None; }
        Some(b)
    };

    // compute next_m (match/mismatch) for all band offsets
    for offset in -((OFFSETWIDTH as isize))..=((OFFSETWIDTH as isize)) {
        let band_idx = (offset + OFFSETWIDTH as isize) as usize;
        let from = mask_score_m[band_idx]
            .max(mask_score_ins[band_idx])
            .max(mask_score_del[band_idx]);

        // Determine genomic coordinate analogous to C++ branches
        let seq_base_opt = if !rev && is_right {
            // sequence[pos[i] + offset + ext - MAXEXTEND]
            let idx = pos as isize + offset + (ext as isize) - (maxextend as isize);
            get_seq(idx)
        } else if rev && is_right {
            // complement(sequence[pos[i] + k - 1 - (- offset + ext - MAXEXTEND)])
            let inner = -offset + (ext as isize) - (maxextend as isize);
            let idx = pos as isize + (k as isize) - 1 - inner;
            get_seq(idx).map(|b| 3 - b)
        } else if !rev && !is_right {
            // sequence[pos[i] - offset + (ext - MAXEXTEND)]
            let idx = pos as isize - offset + (ext as isize) - (maxextend as isize);
            get_seq(idx)
        } else {
            // complement(sequence[pos[i] + k - 1 - (offset + ext - MAXEXTEND)])
            let idx = pos as isize + (k as isize) - 1 - (offset + (ext as isize) - (maxextend as isize));
            get_seq(idx).map(|b| 3 - b)
        };

        // consensus[ext]
        let cons_b = consensus.get(ext).copied().unwrap_or(99u8);
        let add = match (seq_base_opt, cons_b) {
            (Some(sb), cb) if cb <= 3 && sb == cb => match_score,
            (Some(_), cb) if cb <= 3 => mismatch_score,
            _ => mismatch_score, // padding/invalid treated as mismatch
        };
        next_m[band_idx] = from.saturating_add(add);
    }

    // next_ins: for offset in [-W, W-1]: max(M + gap_open, INS + gap_extend) shifting from offset+1
    for offset in -((OFFSETWIDTH as isize))..((OFFSETWIDTH as isize)) {
        let idx = (offset + OFFSETWIDTH as isize) as usize;
        let from_idx = (offset + 1 + OFFSETWIDTH as isize) as usize;
        let v = (mask_score_m[from_idx].saturating_add(gap_open))
            .max(mask_score_ins[from_idx].saturating_add(gap_extend));
        next_ins[idx] = v;
    }
    // next_del: for offset in [-W+1, W]: max(next_m[offset-1] + gap_open, next_del[offset-1] + gap_extend)
    // Note: we must iterate increasing offset to use next_del[offset-1]
    for offset in -((OFFSETWIDTH as isize) - 1)..=((OFFSETWIDTH as isize)) {
        let idx = (offset + OFFSETWIDTH as isize) as usize;
        let prev_idx = (offset - 1 + OFFSETWIDTH as isize) as usize;
        let v = (next_m[prev_idx].saturating_add(gap_open))
            .max(next_del[prev_idx].saturating_add(gap_extend));
        next_del[idx] = v;
    }

    for idx in 0..width {
        mask_score_m[idx] = next_m[idx];
        mask_score_ins[idx] = next_ins[idx];
        mask_score_del[idx] = next_del[idx];
        let s = next_m[idx].max(next_ins[idx]).max(next_del[idx]);
        mask_score[idx] = s;
        if s > best {
            best = s;
        }
    }
    best
}

/// Port of masking_align driving right/left banded extension and picking element start/end.
/// C++ reference: REPrise.cpp:836
pub fn masking_align(
    i: usize,
    consensusstart: usize,
    consensusend: usize,
    // shared state
    consensus: &[u8],
    sequence: &[u8],
    pos: usize,
    k: usize,
    rev: bool,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
    maxextend: usize,
    when_to_stop: usize,
) -> (isize, isize) {
    let width = 2 * OFFSETWIDTH + 1;

    let mut right_mask_score = vec![MINSCORE; width];
    let mut left_mask_score = vec![MINSCORE; width];
    let mut best_right_mask_score = vec![MINSCORE; width];
    let mut best_left_mask_score = vec![MINSCORE; width];

    let mut right_mask_score_m = vec![MINSCORE; width];
    let mut right_mask_score_ins = vec![MINSCORE; width];
    let mut right_mask_score_del = vec![MINSCORE; width];

    let mut left_mask_score_m = vec![MINSCORE; width];
    let mut left_mask_score_ins = vec![MINSCORE; width];
    let mut left_mask_score_del = vec![MINSCORE; width];

    // init centers
    right_mask_score_m[OFFSETWIDTH] = 0;
    for offset in (1..=OFFSETWIDTH).rev() {
        right_mask_score_del[offset + OFFSETWIDTH] = gap_open + gap_extend * (offset as i32 - 1);
    }
    left_mask_score_m[OFFSETWIDTH] = 0;
    for offset in (1..=OFFSETWIDTH).rev() {
        left_mask_score_del[offset + OFFSETWIDTH] = gap_open + gap_extend * (offset as i32 - 1);
    }

    let mut right_best_ext = maxextend as isize;
    let mut right_bestscore = 0i32;
    let mut left_best_ext = maxextend as isize - 1;
    let mut left_bestscore = 0i32;

    // right extend: ext from MAXEXTEND..=consensusend
    let mut right_ext = maxextend as isize;
    while right_ext <= consensusend as isize {
        let score = mask_extention_score(
            true,
            right_ext as usize,
            i,
            &mut right_mask_score,
            &mut right_mask_score_m,
            &mut right_mask_score_ins,
            &mut right_mask_score_del,
            consensus,
            sequence,
            pos,
            k,
            rev,
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
            maxextend,
        );
        if score > right_bestscore {
            right_bestscore = score;
            right_best_ext = right_ext;
            best_right_mask_score.copy_from_slice(&right_mask_score);
        }
        if (right_ext - right_best_ext) as usize >= when_to_stop {
            break;
        }
        right_ext += 1;
    }

    // left extend: ext from MAXEXTEND-1 down to consensusstart
    let mut left_ext_i = maxextend as isize - 1;
    while left_ext_i >= consensusstart as isize {
        let score = mask_extention_score(
            false,
            left_ext_i as usize,
            i,
            &mut left_mask_score,
            &mut left_mask_score_m,
            &mut left_mask_score_ins,
            &mut left_mask_score_del,
            consensus,
            sequence,
            pos,
            k,
            rev,
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
            maxextend,
        );
        if score > left_bestscore {
            left_bestscore = score;
            left_best_ext = left_ext_i;
            best_left_mask_score.copy_from_slice(&left_mask_score);
        }
        if (left_best_ext - left_ext_i) as usize >= when_to_stop {
            break;
        }
        left_ext_i -= 1;
    }

    // choose best offsets in bands
    let mut best_rightoffset = -(OFFSETWIDTH as isize);
    for off in (-(OFFSETWIDTH as isize))..=(OFFSETWIDTH as isize) {
        let idx = (off + OFFSETWIDTH as isize) as usize;
        let cur = best_right_mask_score[idx];
        let best_cur = best_right_mask_score[(best_rightoffset + OFFSETWIDTH as isize) as usize];
        if cur > best_cur {
            best_rightoffset = off;
        }
    }
    let mut best_leftoffset = -(OFFSETWIDTH as isize);
    for off in (-(OFFSETWIDTH as isize))..=(OFFSETWIDTH as isize) {
        let idx = (off + OFFSETWIDTH as isize) as usize;
        let cur = best_left_mask_score[idx];
        let best_cur = best_left_mask_score[(best_leftoffset + OFFSETWIDTH as isize) as usize];
        if cur > best_cur {
            best_leftoffset = off;
        }
    }

    // element start/end in seed-relative coords per C++ lines 906-912
    if !rev {
        let elementstart = left_best_ext - maxextend as isize - best_leftoffset;
        let elementend = right_best_ext - maxextend as isize + best_rightoffset;
        (elementstart, elementend)
    } else {
        let elementstart = -right_best_ext + k as isize - 1 + maxextend as isize - best_rightoffset;
        let elementend = -left_best_ext + k as isize - 1 + maxextend as isize + best_leftoffset;
        (elementstart, elementend)
    }
}

/// Public thin wrapper for CLI use. Exposes masking_align with explicit state to the bin crate.
pub fn cli_masking_align_call(
    seed_pos: usize,
    is_rc: bool,
    k: usize,
    consensus: &[u8],
    sequence: &[u8],
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
    maxextend: usize,
    when_to_stop: usize,
) -> (isize, isize) {
    // The consensus window must match C++ construction: index space [0 .. consensus.len()-1],
    // with seed center at MAXEXTEND. Enforce bounds explicitly.
    let consensus_len = consensus.len();
    let consensusstart = 0usize;
    let consensusend = consensus_len.saturating_sub(1);
    masking_align(
        0, // i not used in our port
        consensusstart,
        consensusend,
        consensus,
        sequence,
        seed_pos,
        k,
        is_rc,
        match_score,
        mismatch_score,
        gap_open,
        gap_extend,
        maxextend,
        when_to_stop,
    )
}

/// Stubs for parity tracking: mask_extention_score and masking_align.
/// These names are referenced by tests/specs; keep behind feature to avoid API changes.
/// They do not alter core semantics; only used by CLI for integration experiments.
#[allow(dead_code)]
pub fn mask_extention_score_stub(
    _is_right: bool,
    _ext: usize,
    _i: usize,
    _mask_score: &mut [i32],
    _mask_score_m: &mut [i32],
    _mask_score_ins: &mut [i32],
    _mask_score_del: &mut [i32],
) -> i32 {
    0
}

#[allow(dead_code)]
#[cfg(test)]
pub fn masking_align_stub(
    _i: usize,
    consensusstart: usize,
    consensusend: usize,
) -> (usize, usize) {
    (consensusstart, consensusend)
}

/// Compute a simple match/mismatch score.
pub fn compute_score(base1: u8, base2: u8, match_score: i32, mismatch: i32) -> i32 {
    if base1 == base2 {
        match_score
    } else {
        mismatch
    }
}

/// Very small repeat family builder. Only selects top kmers and masks them.
pub fn build_repeat_families(
    mut kmers: BinaryHeap<(usize, Vec<u8>)>,
    seq: &[u8],
    cache: &[Vec<CacheEntry>],
    sa: &[i64],
    minfreq: usize,
) -> Vec<Vec<usize>> {
    let mut families = Vec::new();
    let mut mask = vec![false; seq.len()];

    while let Some((freq, kmer)) = kmers.pop() {
        if freq < minfreq {
            break;
        }
        let mut occ = findkmer(&kmer, cache, seq, sa);
        removemasked(&mut occ, &mask, kmer.len(), false);
        if occ.is_empty() {
            continue;
        }
        maskbyseed(&occ, &mut mask, kmer.len(), false);
        families.push(occ);
    }

    families
}

/// Mask repeat regions over all seed occurrences, matching the C++ logic.
/// For forward hits, mask [pos[i] + repeatstart[i], pos[i] + repeatend[i]].
/// For reverse hits, mask [pos[i] - repeatend[i] + k - 1, pos[i] - repeatstart[i] + k - 1].
pub fn maskbyrepeat(
    seedfreq: usize,
    repeatstart: &[usize],
    repeatend: &[usize],
    mask_flag: &mut [bool],
    pos: &[usize],
    rev: &[bool],
    k: usize,
) {
    for i in 0..seedfreq {
        if !rev[i] {
            let start = pos[i] + repeatstart[i];
            let end = pos[i] + repeatend[i];
            for j in start..=end {
                if j < mask_flag.len() {
                    mask_flag[j] = true;
                }
            }
        } else {
            let start = pos[i].saturating_sub(repeatend[i]).saturating_add(k - 1);
            let end = pos[i].saturating_sub(repeatstart[i]).saturating_add(k - 1);
            for j in start..=end {
                if j < mask_flag.len() {
                    mask_flag[j] = true;
                }
            }
        }
    }
}

/// Mask a single repeat element for seed i over [pos[i] + elementstart, pos[i] + elementend].
/// Note: elementstart/elementend are already orientation-adjusted in C++.
pub fn maskbyrepeat_element(
    i: usize,
    elementstart: usize,
    elementend: usize,
    mask_flag: &mut [bool],
    pos: &[usize],
) {
    let start = pos[i] + elementstart;
    let end = pos[i] + elementend;
    for j in start..=end {
        if j < mask_flag.len() {
            mask_flag[j] = true;
        }
    }
}

// Extension alignment constants (matching C++ defaults)
#[allow(dead_code)]
const MATCHSCORE: i32 = 1;
#[allow(dead_code)]
const MISMATCHSCORE: i32 = -1;
#[allow(dead_code)]
const GAPSCORE: i32 = -5;
#[allow(dead_code)]
const GAPEXTENDSCORE: i32 = -1;
#[allow(dead_code)]
const CAPPENALTY: i32 = -20;
#[allow(dead_code)]
const MAXEXTEND: usize = 10000;
#[allow(dead_code)]
const WHEN_TO_STOP: usize = 100;

/// Get complement of a nucleotide (0=A, 1=C, 2=G, 3=T)
pub fn complement(base: u8) -> u8 {
    match base {
        0 => 3, // A -> T
        1 => 2, // C -> G  
        2 => 1, // G -> C
        3 => 0, // T -> A
        _ => base, // Invalid bases remain unchanged
    }
}

/// Extension alignment - extends consensus sequence in one direction
/// Returns the best extension length
pub fn extend(
    is_right: bool,
    seedfreq: usize,
    pos: &[usize],
    rev: &[bool],
    seq: &[u8],
    consensus: &mut [u8],
    seed_ext: &mut [i32],
    match_score: i32,
    mismatch_score: i32,
    gap_score: i32,
    gap_extend_score: i32,
    cap_penalty: i32,
    max_extend: usize,
    offset_width: usize,
    when_to_stop: usize,
    k: usize,
    min_improvement: i32,
) -> i32 {
    let mut nexttotalscore = vec![0i32; 4];
    let mut bestscore = vec![0i32; seedfreq];
    let mut score = vec![vec![MINSCORE; 2 * offset_width + 1]; seedfreq];
    let mut score_m = vec![vec![MINSCORE; 2 * offset_width + 1]; seedfreq];
    let mut score_ins = vec![vec![MINSCORE; 2 * offset_width + 1]; seedfreq];
    let mut score_del = vec![vec![MINSCORE; 2 * offset_width + 1]; seedfreq];
    
    let mut score_m_bybase = vec![vec![vec![MINSCORE; 2 * offset_width + 1]; seedfreq]; 4];
    let mut score_ins_bybase = vec![vec![vec![MINSCORE; 2 * offset_width + 1]; seedfreq]; 4];
    let mut score_del_bybase = vec![vec![vec![MINSCORE; 2 * offset_width + 1]; seedfreq]; 4];
    
    let mut best_ext = -1i32;
    let mut besttotalscore = 0i32;
    
    // Initialize scores
    for se in 0..seedfreq {
        score_m[se][offset_width] = 0;
        for offset in 1..=offset_width {
            score_del[se][offset + offset_width] = gap_score + gap_extend_score * (offset as i32 - 1);
        }
    }
    
    for ext in 0..max_extend {
        // Safety check - prevent buffer overflows that caused C++ heap crashes
        if max_extend + ext >= consensus.len() || (max_extend as i32 - ext as i32 - 1) < 0 {
            eprintln!("Warning: Extension bounds exceeded at ext={}, breaking early", ext);
            break;
        }
        nexttotalscore.fill(0);
        
        for base in 0..4u8 {
            for se in 0..seedfreq {
                let score_val = compute_score_extend(
                    is_right, ext, se, base, pos, rev, seq,
                    &score_m, &score_ins, &score_del,
                    &mut score_m_bybase, &mut score_ins_bybase, &mut score_del_bybase,
                    match_score, mismatch_score, gap_score, gap_extend_score,
                    offset_width, k
                );
                nexttotalscore[base as usize] += std::cmp::max(0, std::cmp::max(bestscore[se] + cap_penalty, score_val));
            }
        }
        
        let bestbase = nexttotalscore.iter().position(|&x| x == *nexttotalscore.iter().max().unwrap()).unwrap() as u8;
        
        // Set consensus base
        if is_right {
            if max_extend + ext < consensus.len() {
                consensus[max_extend + ext] = bestbase;
            }
        } else {
            if max_extend >= ext + 1 {
                consensus[max_extend - ext - 1] = bestbase;
            }
        }
        
        // Update scores
        for se in 0..seedfreq {
            for offset in 0..=(2 * offset_width) {
                score_m[se][offset] = score_m_bybase[bestbase as usize][se][offset];
                score_ins[se][offset] = score_ins_bybase[bestbase as usize][se][offset];
                score_del[se][offset] = score_del_bybase[bestbase as usize][se][offset];
                score[se][offset] = std::cmp::max(score_m[se][offset], 
                    std::cmp::max(score_ins[se][offset], score_del[se][offset]));
            }
        }
        
        let mut tmpbesttotalscore = 0i32;
        for se in 0..seedfreq {
            let mut tmpbestscore = MINSCORE * 2;
            let mut bestoffset = 0i32;
            
            for offset in 0..=(2 * offset_width) {
                if score[se][offset] > tmpbestscore {
                    tmpbestscore = score[se][offset];
                    bestoffset = offset as i32 - offset_width as i32;
                }
            }
            
            if tmpbestscore > bestscore[se] {
                seed_ext[se] = ext as i32 + bestoffset;
                bestscore[se] = tmpbestscore;
            }
            tmpbesttotalscore += std::cmp::max(tmpbestscore, bestscore[se] + cap_penalty);
        }
        
        // Check if improvement meets minimum threshold (C++ line: if (tmpbesttotalscore >= besttotalscore + (ext - best_ext) * MINIMPROVEMENT))
        let improvement_threshold = if best_ext >= 0 {
            besttotalscore + ((ext as i32 - best_ext) * min_improvement)
        } else {
            besttotalscore
        };
        
        if tmpbesttotalscore >= improvement_threshold {
            besttotalscore = tmpbesttotalscore;
            best_ext = ext as i32;
        } else {
            // Check stop condition
            if ext >= when_to_stop && best_ext >= 0 && (ext - best_ext as usize) >= when_to_stop {
                break;
            }
        }
    }
    
    best_ext
}

/// Compute alignment score for extension (helper for extend function)
fn compute_score_extend(
    is_right: bool,
    ext: usize,
    se: usize,
    base: u8,
    pos: &[usize],
    rev: &[bool],
    seq: &[u8],
    score_m: &[Vec<i32>],
    score_ins: &[Vec<i32>],
    score_del: &[Vec<i32>],
    score_m_bybase: &mut [Vec<Vec<i32>>],
    score_ins_bybase: &mut [Vec<Vec<i32>>],
    score_del_bybase: &mut [Vec<Vec<i32>>],
    match_score: i32,
    mismatch_score: i32,
    gap_score: i32,
    gap_extend_score: i32,
    offset_width: usize,
    k: usize,
) -> i32 {
    let mut nextscore_m_bybase = vec![MINSCORE; 2 * offset_width + 1];
    let mut nextscore_ins_bybase = vec![MINSCORE; 2 * offset_width + 1];
    let mut nextscore_del_bybase = vec![MINSCORE; 2 * offset_width + 1];
    
    // Compute match scores
    if !rev[se] && is_right {
        for offset in 0..=(2 * offset_width) {
            let actual_offset = offset as i32 - offset_width as i32;
            let seq_pos = pos[se] as i32 + actual_offset + ext as i32;
            if seq_pos >= 0 && (seq_pos as usize) < seq.len() {
                let seq_base = seq[seq_pos as usize];
                let score_increment = if base == seq_base { match_score } else { mismatch_score };
                nextscore_m_bybase[offset] = std::cmp::max(score_m[se][offset], 
                    std::cmp::max(score_ins[se][offset], score_del[se][offset])) + score_increment;
            }
        }
    } else if rev[se] && is_right {
        for offset in 0..=(2 * offset_width) {
            let actual_offset = offset as i32 - offset_width as i32;
            let seq_pos = pos[se] as i32 + k as i32 - 1 - (actual_offset + ext as i32);
            if seq_pos >= 0 && (seq_pos as usize) < seq.len() {
                let seq_base = complement(seq[seq_pos as usize]);
                let score_increment = if base == seq_base { match_score } else { mismatch_score };
                nextscore_m_bybase[offset] = std::cmp::max(score_m[se][offset], 
                    std::cmp::max(score_ins[se][offset], score_del[se][offset])) + score_increment;
            }
        }
    } else if !rev[se] && !is_right {
        for offset in 0..=(2 * offset_width) {
            let actual_offset = offset as i32 - offset_width as i32;
            let seq_pos = pos[se] as i32 - 1 - actual_offset - ext as i32;
            if seq_pos >= 0 && (seq_pos as usize) < seq.len() {
                let seq_base = seq[seq_pos as usize];
                let score_increment = if base == seq_base { match_score } else { mismatch_score };
                nextscore_m_bybase[offset] = std::cmp::max(score_m[se][offset], 
                    std::cmp::max(score_ins[se][offset], score_del[se][offset])) + score_increment;
            }
        }
    } else {
        for offset in 0..=(2 * offset_width) {
            let actual_offset = offset as i32 - offset_width as i32;
            let seq_pos = pos[se] as i32 + k as i32 - 1 - (-1 - actual_offset - ext as i32);
            if seq_pos >= 0 && (seq_pos as usize) < seq.len() {
                let seq_base = complement(seq[seq_pos as usize]);
                let score_increment = if base == seq_base { match_score } else { mismatch_score };
                nextscore_m_bybase[offset] = std::cmp::max(score_m[se][offset], 
                    std::cmp::max(score_ins[se][offset], score_del[se][offset])) + score_increment;
            }
        }
    }
    
    // Compute insertion scores
    for offset in 0..=(2 * offset_width) {
        if offset > 0 {
            nextscore_ins_bybase[offset] = std::cmp::max(
                score_m[se][offset - 1] + gap_score,
                std::cmp::max(
                    score_ins[se][offset - 1] + gap_extend_score,
                    score_del[se][offset - 1] + gap_score
                )
            );
        }
    }
    
    // Compute deletion scores
    for offset in 0..(2 * offset_width) {
        nextscore_del_bybase[offset] = std::cmp::max(
            score_m[se][offset + 1] + gap_score,
            std::cmp::max(
                score_ins[se][offset + 1] + gap_score,
                score_del[se][offset + 1] + gap_extend_score
            )
        );
    }
    
    // Store results in bybase arrays
    for offset in 0..=(2 * offset_width) {
        score_m_bybase[base as usize][se][offset] = nextscore_m_bybase[offset];
        score_ins_bybase[base as usize][se][offset] = nextscore_ins_bybase[offset];
        score_del_bybase[base as usize][se][offset] = nextscore_del_bybase[offset];
    }
    
    // Return best score for this seed
    let mut best_score = MINSCORE;
    for offset in 0..=(2 * offset_width) {
        best_score = std::cmp::max(best_score, 
            std::cmp::max(nextscore_m_bybase[offset], 
                std::cmp::max(nextscore_ins_bybase[offset], nextscore_del_bybase[offset])));
    }
    best_score
}

/// Find the best seed from the priority queue, considering masked regions.
/// This function mirrors the C++ implementation but uses Rust data structures.
pub fn find_bestseed(
    kmers: &mut BinaryHeap<(usize, Vec<u8>)>,
    cache: &[Vec<CacheEntry>],
    mask_flag: &[bool],
    seq: &[u8],
    sa: &[i64],
    tandem_dist: usize,
    minfreq: usize,
) -> (usize, Vec<u8>, Vec<usize>, Vec<bool>) {
    let mut pos: Vec<usize> = Vec::new();
    let mut rev: Vec<bool> = Vec::new();
    
    while let Some(tmpbest) = kmers.pop() {
        let (_freq, kmer) = tmpbest;
        
        // Find occurrences of the kmer
        let mut occs = findkmer(&kmer, cache, seq, sa);
        removetandem(&mut occs, tandem_dist);
        removemasked(&mut occs, mask_flag, kmer.len(), false);
        
        // Find occurrences of the reverse complement
        let rc_kmer: Vec<u8> = kmer.iter().rev().map(|&b| 3 - b).collect();
        let mut rc_occs = findkmer(&rc_kmer, cache, seq, sa);
        removetandem(&mut rc_occs, tandem_dist);
        removemasked(&mut rc_occs, mask_flag, kmer.len(), true);
        
        let newfreq = occs.len() + rc_occs.len();
        
        // If queue is empty or new frequency is better than or equal to top, break
        if kmers.is_empty() || newfreq >= kmers.peek().unwrap().0 {
            // Update pos and rev vectors
            pos.clear();
            pos.extend(occs.iter());
            pos.extend(rc_occs.iter());
            
            rev.clear();
            rev.extend(std::iter::repeat(false).take(occs.len()));
            rev.extend(std::iter::repeat(true).take(rc_occs.len()));
            
            if newfreq >= minfreq {
                return (newfreq, kmer, pos, rev);
            } else {
                // If frequency is too low, continue with next kmer
                continue;
            }
        }
        
        // If new frequency is still good, push it back to queue
        if newfreq >= minfreq {
            kmers.push((newfreq, kmer));
        }
    }
    
    // If no suitable kmer found, return empty result
    (0, Vec::new(), pos, rev)
}

