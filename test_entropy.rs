fn compute_entropy(kmer: &[u8]) -> f64 {
    let mut count = [0; 4];
    for &base in kmer {
        if base <= 3 {
            count[base as usize] += 1;
        }
    }
    let len = kmer.len() as f64;
    let mut entropy = 0.0;
    for &c in &count {
        if c > 0 {
            let p = c as f64 / len;
            entropy -= p * p.log2();
        }
    }
    entropy
}

fn main() {
    // Test k-mer "TAGAAATAT" from C++ output
    let kmer = vec![3, 0, 2, 0, 0, 0, 3, 0, 3]; // T=3, A=0, G=2
    println!("Entropy of TAGAAATAT: {}", compute_entropy(&kmer));
    
    // Test k-mer "AAAGTAGCC" from Rust output  
    let kmer2 = vec![0, 0, 0, 2, 3, 0, 2, 1, 1]; // A=0, G=2, T=3, C=1
    println!("Entropy of AAAGTAGCC: {}", compute_entropy(&kmer2));
    
    println!("MAX_ENTROPY threshold: -0.7");
}