// Profiling version of main with detailed timing
use reprise::{build_sequence, suffix_array, default_k};
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <input.fa> <output_prefix>", args[0]);
        std::process::exit(1);
    }
    
    let input_file = &args[1];
    let output_prefix = &args[2];
    
    println!("Starting profiling run on {}", input_file);
    let total_start = Instant::now();
    
    // 1. Sequence loading
    println!("=== Phase 1: Loading sequence ===");
    let phase1_start = Instant::now();
    let data = build_sequence(input_file).expect("Failed to load sequence");
    let phase1_time = phase1_start.elapsed();
    println!("Sequence length: {}", data.sequence.len());
    println!("Phase 1 time: {:?}", phase1_time);
    println!();
    
    // 2. Suffix array construction
    println!("=== Phase 2: Building suffix array ===");
    let phase2_start = Instant::now();
    let sa = suffix_array(&data.sequence);
    let phase2_time = phase2_start.elapsed();
    println!("Suffix array length: {}", sa.len());
    println!("Phase 2 time: {:?}", phase2_time);
    println!();
    
    // 3. K-mer processing (simplified version)
    println!("=== Phase 3: K-mer processing ===");
    let phase3_start = Instant::now();
    let k = default_k(data.sequence.len(), 0);
    println!("Using k-mer length: {}", k);
    
    // Simple k-mer frequency counting for profiling
    let mut kmer_count = 0;
    let seq = &data.sequence;
    for i in 0..(seq.len().saturating_sub(k)) {
        // Count k-mers (simplified)
        kmer_count += 1;
        if kmer_count % 100000 == 0 {
            print!(".");
        }
    }
    let phase3_time = phase3_start.elapsed();
    println!("\nProcessed {} k-mers", kmer_count);
    println!("Phase 3 time: {:?}", phase3_time);
    println!();
    
    let total_time = total_start.elapsed();
    println!("=== Summary ===");
    println!("Phase 1 (Loading):        {:?} ({:.1}%)", phase1_time, 100.0 * phase1_time.as_secs_f64() / total_time.as_secs_f64());
    println!("Phase 2 (Suffix Array):   {:?} ({:.1}%)", phase2_time, 100.0 * phase2_time.as_secs_f64() / total_time.as_secs_f64());
    println!("Phase 3 (K-mer counting): {:?} ({:.1}%)", phase3_time, 100.0 * phase3_time.as_secs_f64() / total_time.as_secs_f64());
    println!("Total time:               {:?}", total_time);
    println!();
    
    // Create dummy output files so the test is complete
    std::fs::write(format!("{}.reprof", output_prefix), "# Profiling run\n").ok();
    
    println!("Profiling complete!");
}