use std::fs::File;
use std::io::{self, BufRead, BufReader};

use crate::utils::{char_to_num, complement};

const PAD_LENGTH: usize = 11000;

#[derive(Debug)]
pub struct Genome {
    pub sequence: Vec<u8>,
    pub chrtable: Vec<(String, usize)>,
}

pub fn read_fasta(path: &str) -> io::Result<Genome> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut sequence = Vec::new();
    let mut chrtable = Vec::new();
    chrtable.push(("unknown".to_string(), 0));

    // initial padding
    sequence.extend(std::iter::repeat(99u8).take(PAD_LENGTH));

    let mut length = PAD_LENGTH;
    let mut name: Option<String> = None;

    for line_res in reader.lines() {
        let line = line_res?;
        if line.starts_with('>') {
            // finalize previous chr
            if name.is_some() {
                chrtable.push(("padding".to_string(), length));
                sequence.extend(std::iter::repeat(99u8).take(PAD_LENGTH));
                length += PAD_LENGTH;
            }
            let chrname = line[1..].split_whitespace().next().unwrap_or("").to_string();
            chrtable.push((chrname.clone(), length));
            name = Some(chrname);
        } else {
            for c in line.trim().chars() {
                sequence.push(char_to_num(c));
                length += 1;
            }
        }
    }
    // final padding
    chrtable.push(("padding".to_string(), length));
    sequence.extend(std::iter::repeat(99u8).take(PAD_LENGTH));

    Ok(Genome { sequence, chrtable })
}

/// Compute log-base-e entropy of a kmer represented as numeric bases (0..3).
pub fn compute_entropy(kmer: &[u8]) -> f64 {
    let mut count = [0usize; 4];
    for &b in kmer {
        if b < 4 {
            count[b as usize] += 1;
        }
    }
    let k = kmer.len() as f64;
    let mut ans = 0.0f64;
    for &c in &count {
        if c == 0 {
            continue;
        }
        let y = c as f64 / k;
        ans += y * y.ln();
    }
    ans
}

/// Compute the default k‑mer length based on genome length and allowed edits.
pub fn default_k(len: usize, kmerdist: usize) -> usize {
    const MAX: usize = 40;
    let mut v = vec![vec![0u128; MAX]; MAX];
    for i in 0..MAX {
        v[i][0] = 1;
        v[i][i] = 1;
    }
    for kk in 1..MAX {
        for j in 1..kk {
            v[kk][j] = v[kk - 1][j - 1] + v[kk - 1][j];
        }
    }

    let mut l = 2usize;
    while l <= MAX {
        let mut comb = 0f64;
        for d in 0..=kmerdist {
            if d <= l {
                comb += v[l][d] as f64 * 3f64.powi(d as i32);
            }
        }
        let e = len as f64 * comb / 4f64.powi(l as i32);
        if e < 1.0 {
            break;
        }
        l += 1;
    }
    l + 1
}

/// Naive suffix array construction used while porting the algorithm.
pub fn build_suffix_array(seq: &[u8]) -> Vec<usize> {
    let mut sa: Vec<usize> = (0..seq.len()).collect();
    sa.sort_by(|&a, &b| seq[a..].cmp(&seq[b..]));
    sa
}

/// Find exact kmer occurrences using a suffix array.
pub fn find_kmer(seq: &[u8], sa: &[usize], query: &[u8]) -> Vec<usize> {
    let qlen = query.len();
    sa.iter()
        .copied()
        .filter(|&pos| pos + qlen <= seq.len() && &seq[pos..pos + qlen] == query)
        .collect()
}

use std::collections::HashMap;

/// Count the occurrences of every k-mer in the sequence.
pub fn kmer_frequencies(seq: &[u8], k: usize) -> HashMap<Vec<u8>, usize> {
    let mut counts: HashMap<Vec<u8>, usize> = HashMap::new();
    if k == 0 || seq.len() < k {
        return counts;
    }
    for window in seq.windows(k) {
        counts.entry(window.to_vec()).and_modify(|c| *c += 1).or_insert(1);
    }
    counts
}

/// Find the k-mer with the highest frequency.
pub fn best_kmer(seq: &[u8], k: usize) -> Option<(Vec<u8>, usize)> {
    kmer_frequencies(seq, k)
        .into_iter()
        .max_by_key(|(_, c)| *c)
}

/// Remove occurrences that are within `dist` bases of the previous one.
pub fn remove_tandem(occs: &mut Vec<usize>, dist: usize) {
    occs.sort_unstable();
    let mut prev: Option<usize> = None;
    occs.retain(|&occ| {
        if let Some(p) = prev {
            if occ < p + dist {
                false
            } else {
                prev = Some(occ);
                true
            }
        } else {
            prev = Some(occ);
            true
        }
    });
}

/// Remove occurrences overlapping masked regions.
pub fn remove_masked(occs: &mut Vec<usize>, mask: &[bool], is_rc: bool, k: usize) {
    occs.retain(|&occ| {
        if !is_rc {
            (0..k).all(|i| !mask[occ + i])
        } else {
            (0..k).all(|i| !mask[occ + k - 1 - i])
        }
    });
}

/// Mark seed positions as masked.
pub fn mask_by_seed(occs: &[usize], mask: &mut [bool], is_rc: bool, k: usize) {
    if !is_rc {
        for &o in occs {
            for i in 0..k {
                mask[o + i] = true;
            }
        }
    } else {
        for &o in occs {
            for i in 0..k {
                mask[o + k - 1 - i] = true;
            }
        }
    }
}

/// Return the chromosome name and start index for a sequence position.
pub fn chrtracer(table: &[(String, usize)], pos: usize) -> (String, usize) {
    let mut result = table.last().cloned().unwrap_or_else(|| ("".into(), 0));
    for i in 0..table.len() {
        if pos < table[i].1 {
            if i > 0 {
                result = table[i - 1].clone();
            }
            break;
        }
    }
    result
}
