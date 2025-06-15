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
