pub mod alg;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// Number used to pad sequences between contigs, mirroring the C++ implementation.
pub const PADLENGTH: usize = 11000;

/// Mapping from nucleotide characters to numeric codes used in the C++ version.
fn char_to_num(c: char) -> u8 {
    match c {
        'A' | 'a' => 0,
        'C' | 'c' => 1,
        'G' | 'g' => 2,
        'T' | 't' => 3,
        'N' | 'n' | 'x'
        | 'R' | 'r' | 'Y' | 'y' | 'M' | 'm'
        | 'K' | 'k' | 'W' | 'w' | 'S' | 's'
        | 'B' | 'b' | 'D' | 'd' | 'H' | 'h' | 'V' | 'v' => 99,
        _ => 99,
    }
}

/// Sequence and chromosome table returned from `build_sequence`.
pub struct SequenceData {
    pub sequence: Vec<u8>,
    pub chrtable: Vec<(String, usize)>,
}

/// Read a FASTA file and generate the padded numeric sequence.
/// This follows the logic of `build_sequence` from `REPrise.cpp`.
pub fn build_sequence<P: AsRef<Path>>(path: P) -> io::Result<SequenceData> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut sequence: Vec<u8> = Vec::new();
    let mut chrtable: Vec<(String, usize)> = Vec::new();

    // Initial padding and unknown entry.
    chrtable.push(("unknown".to_string(), 0));
    sequence.extend(std::iter::repeat(99u8).take(PADLENGTH));

    for line_res in reader.lines() {
        let line = line_res?;
        if line.starts_with('>') {
            // Padding between contigs
            chrtable.push(("padding".to_string(), sequence.len()));
            sequence.extend(std::iter::repeat(99u8).take(PADLENGTH));

            let mut name = line[1..].trim().to_string();
            name = name
                .chars()
                .map(|c| if c == ' ' || c == '\t' { '_' } else { c })
                .collect();
            chrtable.push((name, sequence.len()));
        } else {
            for ch in line.trim().chars() {
                if (ch as u32) > 64 {
                    sequence.push(char_to_num(ch));
                }
            }
        }
    }

    // Trailing padding at EOF.
    sequence.extend(std::iter::repeat(99u8).take(PADLENGTH));

    Ok(SequenceData { sequence, chrtable })
}

/// Very simple (and slow) suffix array implementation used for testing.
pub fn suffix_array(seq: &[u8]) -> Vec<i64> {
    let mut sa: Vec<i64> = (0..seq.len() as i64).collect();
    sa.sort_by(|&a, &b| seq[a as usize..].cmp(&seq[b as usize..]));
    sa
}

/// Convert numeric nucleotide code to character (mirrors C++ num_to_char).
pub fn num_to_char(z: u8) -> char {
    match z {
        0 => 'A',
        1 => 'C', 
        2 => 'G',
        3 => 'T',
        _ => 'N',
    }
}

/// Get complement of nucleotide (mirrors C++ complement).
pub fn complement(c: u8) -> u8 {
    match c {
        0 => 3, // A -> T
        1 => 2, // C -> G  
        2 => 1, // G -> C
        3 => 0, // T -> A
        _ => 99, // Invalid -> N
    }
}

/// Get reverse complement of sequence (mirrors C++ reverse_complement).
pub fn reverse_complement(query: &[u8]) -> Vec<u8> {
    query.iter().rev().map(|&b| complement(b)).collect()
}

/// Compute entropy of k-mer (mirrors C++ compute_entropy).
pub fn compute_entropy(kmer: &[u8]) -> f64 {
    let mut count = [0; 4];
    for &base in kmer {
        if base <= 3 {
            count[base as usize] += 1;
        }
    }
    
    let mut answer = 0.0;
    let k = kmer.len() as f64;
    for &cnt in &count {
        if cnt > 0 {
            let y = cnt as f64 / k;
            answer += y * y.ln();
        }
    }
    answer
}

/// Map absolute sequence position to chromosome name and offset (mirrors C++ chrtracer).
pub fn chrtracer(chrtable: &[(String, usize)], stringpos: usize) -> (String, usize) {
    let mut result = &chrtable[chrtable.len() - 1];
    for entry in chrtable {
        if stringpos < entry.1 {
            // Found first entry with start > pos, so use previous entry
            let idx = chrtable.iter().position(|x| x.1 == entry.1).unwrap();
            if idx > 0 {
                result = &chrtable[idx - 1];
            }
            break;
        }
        result = entry;
    }
    (result.0.clone(), result.1)
}

/// Calculate default k-mer length (mirrors C++ default_k).
pub fn default_k(len: usize, kmer_dist: usize) -> usize {
    // Build Pascal's triangle for combinations
    let mut v = vec![vec![0i64; 40]; 40];
    for i in 0..v.len() {
        v[i][0] = 1;
        v[i][i] = 1;
    }
    for kk in 1..v.len() {
        for j in 1..kk {
            v[kk][j] = v[kk-1][j-1] + v[kk-1][j];
        }
    }

    for l in 2..=40 {
        let mut comb = 0.0;
        for d in 0..=kmer_dist {
            if d < v[l].len() {
                comb += v[l][d] as f64 * 3.0_f64.powi(d as i32);
            }
        }
        let e = len as f64 * comb / 4.0_f64.powi(l as i32);
        if e < 1.0 {
            return l + 1;
        }
    }
    41
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_and_sa() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../test/tst.fa");
        let SequenceData { sequence, .. } = build_sequence(&path).expect("load test fasta");
        assert_eq!(sequence.len(), 46000);

        let sa = suffix_array(&sequence);
        assert_eq!(sa.len(), sequence.len());

        for i in 1..sa.len() {
            let a = &sequence[sa[i - 1] as usize..];
            let b = &sequence[sa[i] as usize..];
            assert!(a <= b);
        }
    }
}
