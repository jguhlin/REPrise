fn compute_entropy(kmer: &[u8]) -> f64 {
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
            answer += y * y.ln();  // Natural log like C++
        }
    }
    answer
}

fn main() {
    // Test k-mer "TAGAAATAT" from C++ output
    let kmer1 = vec![3, 0, 2, 0, 0, 0, 3, 0, 3]; // T=3, A=0, G=2
    let e1 = compute_entropy(&kmer1);
    println!("Entropy of TAGAAATAT: {}", e1);
    println!("Passes filter (> -0.7)? {}", e1 > -0.7);
    
    // Test k-mer "AAAGTAGCC" from Rust output  
    let kmer2 = vec![0, 0, 0, 2, 3, 0, 2, 1, 1]; // A=0, G=2, T=3, C=1
    let e2 = compute_entropy(&kmer2);
    println!("Entropy of AAAGTAGCC: {}", e2);
    println!("Passes filter (> -0.7)? {}", e2 > -0.7);
    
    println!("\nNote: Higher entropy (closer to 0) = MORE complex");
    println!("Filter removes k-mers with entropy > -0.7 (low complexity)");
}