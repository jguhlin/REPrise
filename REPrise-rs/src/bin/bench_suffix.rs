// Benchmark script to compare suffix array implementations
use reprise::{build_sequence, suffix_array};
use std::time::Instant;

fn naive_suffix_array(seq: &[u8]) -> Vec<i64> {
    let mut sa: Vec<i64> = (0..seq.len() as i64).collect();
    sa.sort_by(|&a, &b| seq[a as usize..].cmp(&seq[b as usize..]));
    sa
}

fn main() {
    println!("Suffix Array Performance Benchmark");
    println!("==================================");

    // Test with the standard test file
    let test_file = "../test/tst.fa";
    
    println!("Loading sequence from {}", test_file);
    let sequence_data = build_sequence(test_file).expect("Failed to load test sequence");
    let seq = &sequence_data.sequence;
    
    println!("Sequence length: {}", seq.len());
    println!();

    // Test naive implementation
    println!("Testing naive O(n² log n) implementation:");
    let start = Instant::now();
    let _naive_sa = naive_suffix_array(seq);
    let naive_time = start.elapsed();
    println!("Time: {:?}", naive_time);
    println!();

    // Test suffix crate implementation  
    println!("Testing suffix crate (SAIS) implementation:");
    let start = Instant::now();
    let _suffix_sa = suffix_array(seq);
    let suffix_time = start.elapsed();
    println!("Time: {:?}", suffix_time);
    println!();

    // Compare results
    println!("Performance Comparison:");
    println!("=======================");
    println!("Naive implementation:   {:?}", naive_time);
    println!("Suffix crate:          {:?}", suffix_time);
    
    if naive_time > suffix_time {
        let speedup = naive_time.as_secs_f64() / suffix_time.as_secs_f64();
        println!("Speedup: {:.2}x faster", speedup);
    } else {
        let slowdown = suffix_time.as_secs_f64() / naive_time.as_secs_f64();
        println!("Slowdown: {:.2}x slower", slowdown);
    }

    // Test with larger sequence if available
    let large_file = "../data/ecoli.fasta";
    if std::path::Path::new(large_file).exists() {
        println!();
        println!("Testing with larger E. coli sequence:");
        println!("=====================================");
        
        let large_data = build_sequence(large_file).expect("Failed to load E. coli sequence");
        let large_seq = &large_data.sequence;
        println!("E. coli sequence length: {}", large_seq.len());

        // Only test suffix crate on large sequence (naive would be too slow)
        println!("Testing suffix crate on large sequence:");
        let start = Instant::now();
        let _large_sa = suffix_array(large_seq);
        let large_time = start.elapsed();
        println!("Time: {:?}", large_time);
    }
}