// Benchmark just the suffix array construction on E. coli
use reprise::{build_sequence, suffix_array};
use std::time::Instant;

fn main() {
    println!("E. coli Suffix Array Benchmark");
    println!("===============================");

    let ecoli_file = "../data/ecoli.fasta";
    
    println!("Loading E. coli sequence...");
    let sequence_data = build_sequence(ecoli_file).expect("Failed to load E. coli sequence");
    let seq = &sequence_data.sequence;
    
    println!("E. coli sequence length: {}", seq.len());
    println!();

    // Test suffix array construction only
    println!("Testing suffix array construction with suffix crate:");
    let start = Instant::now();
    let _sa = suffix_array(seq);
    let duration = start.elapsed();
    println!("Time: {:?}", duration);
    println!("Rate: {:.2} bases/ms", seq.len() as f64 / duration.as_millis() as f64);
}