pub fn char_to_num(c: char) -> u8 {
    match c {
        'A' | 'a' => 0,
        'C' | 'c' => 1,
        'G' | 'g' => 2,
        'T' | 't' => 3,
        'N' | 'n' => 99,
        'R' | 'r' | 'Y' | 'y' | 'M' | 'm' |
        'K' | 'k' | 'W' | 'w' | 'S' | 's' |
        'B' | 'b' | 'D' | 'd' | 'H' | 'h' |
        'V' | 'v' => 99,
        _ => 99,
    }
}

pub fn num_to_char(n: u8) -> char {
    match n {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        _ => 'N',
    }
}

pub fn complement(n: u8) -> u8 {
    match n {
        0 => 3,
        1 => 2,
        2 => 1,
        3 => 0,
        _ => 99,
    }
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}
