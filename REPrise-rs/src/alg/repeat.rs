use rayon::prelude::*;
use std::collections::BinaryHeap;
use std::sync::{Arc, Mutex};

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
/// The returned table mirrors the C++ `store_cache` output but uses
/// Rust data structures and parallel iteration.
pub fn store_cache(
    edit_distance: u8,
    cache_len: usize,
    seq: &[u8],
    sa: &[i64],
) -> Vec<Vec<CacheEntry>> {
    let size = 1usize << (cache_len * 2);
    let cache: Vec<_> = (0..size).map(|_| Mutex::new(Vec::new())).collect();
    let cache = Arc::new(cache);

    (0..size).into_par_iter().for_each(|id| {
        let mut kmer = vec![0u8; cache_len];
        for i in 0..cache_len {
            kmer[cache_len - 1 - i] = ((id >> (i * 2)) & 3) as u8;
        }
        let mut begin = None;
        let mut end = None;
        for (i, &p) in sa.iter().enumerate() {
            let p = p as usize;
            if p + cache_len > seq.len() {
                continue;
            }
            if &seq[p..p + cache_len] == &kmer[..] {
                if begin.is_none() {
                    begin = Some(i);
                }
                end = Some(i + 1);
            } else if begin.is_some() {
                break;
            }
        }
        if let (Some(b), Some(e)) = (begin, end) {
            let mut slot = cache[id].lock().unwrap();
            slot.push((b, e, cache_len as u8, edit_distance));
        }
    });

    Arc::try_unwrap(cache)
        .unwrap()
        .into_iter()
        .map(|m| m.into_inner().unwrap())
        .collect()
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
pub fn removemasked(occs: &mut Vec<usize>, mask: &[bool], k: usize, is_rc: bool) {
    occs.retain(|&p| {
        let range = if !is_rc {
            p..p + k
        } else {
            p.saturating_sub(k - 1)..=p
        };
        !range.into_iter().any(|i| i < mask.len() && mask[i])
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
            for i in 0..k {
                if p + k - 1 >= i {
                    mask[p + k - 1 - i] = true;
                }
            }
        }
    }
}

/// Very naive kmer search using the suffix array cache table.
pub fn findkmer(query: &[u8], cache: &[Vec<CacheEntry>], seq: &[u8], sa: &[i64]) -> Vec<usize> {
    let mut idx = 0usize;
    for (i, &b) in query.iter().enumerate() {
        idx |= (b as usize) << (i * 2);
    }
    let mut out = Vec::new();
    if idx < cache.len() {
        for &(b, e, _, dist) in &cache[idx] {
            let mut tmp = Vec::new();
            sa_search(query, seq, sa, b, e, dist as i32, &mut tmp);
            for (s, _) in tmp {
                out.push(sa[s] as usize);
            }
        }
    }
    out
}

/// Build sorted kmers together with their frequencies.
/// This function is heavily simplified compared to the C++ version.
pub fn build_sortedkmers(
    k: usize,
    seq: &[u8],
    cache: &[Vec<CacheEntry>],
    sa: &[i64],
    minfreq: usize,
) -> BinaryHeap<(usize, Vec<u8>)> {
    let len = seq.len();
    let results: Vec<_> = (0..len.saturating_sub(k))
        .into_par_iter()
        .filter_map(|i| {
            let kmer = &seq[i..i + k];
            if kmer.iter().any(|&b| b > 3) {
                return None;
            }
            let mut freq = findkmer(kmer, cache, seq, sa).len();
            let rc: Vec<u8> = kmer
                .iter()
                .rev()
                .map(|&b| 3 - b)
                .collect();
            freq += findkmer(&rc, cache, seq, sa).len();
            if freq >= minfreq {
                Some((freq, kmer.to_vec()))
            } else {
                None
            }
        })
        .collect();

    let mut heap = BinaryHeap::new();
    for r in results {
        heap.push(r);
    }
    heap
}

/// Placeholder extension routine returning zero extension.
pub fn extend(_is_right: bool, seedfreq: usize, seed_ext: &mut Vec<usize>) -> i32 {
    for s in seed_ext.iter_mut().take(seedfreq) {
        *s = 0;
    }
    0
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

