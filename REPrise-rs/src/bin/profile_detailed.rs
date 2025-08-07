// Detailed profiling of the main algorithm steps
use reprise::{build_sequence, suffix_array, default_k, reverse_complement, num_to_char};
use std::time::Instant;
use std::io::{Write, BufWriter};
use std::fs::File;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <input.fa> <output_prefix>", args[0]);
        std::process::exit(1);
    }
    
    let input_file = &args[1];
    let output_prefix = &args[2];
    
    println!("Starting detailed profiling run on {}", input_file);
    let total_start = Instant::now();
    
    // Phase 1: Loading
    println!("=== Phase 1: Loading sequence ===");
    let phase1_start = Instant::now();
    let data = build_sequence(input_file).expect("failed to read FASTA");
    let phase1_time = phase1_start.elapsed();
    println!("Sequence length: {}", data.sequence.len());
    println!("Phase 1 time: {:?}", phase1_time);
    
    // Phase 2: Suffix array
    println!("=== Phase 2: Building suffix array ===");
    let phase2_start = Instant::now();
    let sa = suffix_array(&data.sequence);
    let phase2_time = phase2_start.elapsed();
    println!("Suffix array length: {}", sa.len());
    println!("Phase 2 time: {:?}", phase2_time);
    
    // Phase 3: Parameters
    let k = default_k(data.sequence.len(), 0);
    let min_freq = 3;
    println!("Using k-mer length: {}", k);
    
    // Phase 4: Cache building
    println!("=== Phase 3: Building k-mer cache ===");
    let phase3_start = Instant::now();
    let cache = reprise::alg::repeat::store_cache(0, k, &data.sequence, &sa);
    let phase3_time = phase3_start.elapsed();
    println!("Cache size: {}", cache.len());
    println!("Phase 3 time: {:?}", phase3_time);
    
    // Phase 5: Sorted k-mers
    println!("=== Phase 4: Building sorted k-mers ===");
    let phase4_start = Instant::now();
    let mut kmers = reprise::alg::repeat::build_sortedkmers(
        k,
        &data.sequence,
        &cache,
        &sa,
        min_freq,
    );
    let phase4_time = phase4_start.elapsed();
    println!("K-mers heap size: {}", kmers.len());
    println!("Phase 4 time: {:?}", phase4_time);
    
    // Phase 6: Process a few k-mers
    println!("=== Phase 5: Processing k-mers (first 10) ===");
    let phase5_start = Instant::now();
    let mut mask = vec![false; data.sequence.len()];
    let tandem_dist = 500;
    let mut processed = 0;
    
    while let Some((freq, kmer)) = kmers.pop() {
        if freq < min_freq { break; }
        if processed >= 10 { break; } // Only process first 10 for profiling
        
        let kmer_start = Instant::now();
        
        // Forward occurrences
        let mut occ_f = reprise::alg::repeat::findkmer(&kmer, &cache, &data.sequence, &sa);
        reprise::alg::repeat::removetandem(&mut occ_f, tandem_dist);
        reprise::alg::repeat::removemasked(&mut occ_f, &mask, k, false);

        // Reverse complement occurrences
        let rc_kmer = reverse_complement(&kmer);
        let mut occ_r = reprise::alg::repeat::findkmer(&rc_kmer, &cache, &data.sequence, &sa);
        reprise::alg::repeat::removetandem(&mut occ_r, tandem_dist);
        reprise::alg::repeat::removemasked(&mut occ_r, &mask, k, true);

        let total_freq = occ_f.len() + occ_r.len();
        if total_freq < min_freq { continue; }

        // Mask the used positions
        reprise::alg::repeat::maskbyseed(&occ_f, &mut mask, k, false);
        reprise::alg::repeat::maskbyseed(&occ_r, &mut mask, k, true);
        
        let kmer_time = kmer_start.elapsed();
        println!("  K-mer {} (freq={}): {:?}", processed, freq, kmer_time);
        
        processed += 1;
    }
    
    let phase5_time = phase5_start.elapsed();
    println!("Processed {} k-mers", processed);
    println!("Phase 5 time: {:?}", phase5_time);
    
    let total_time = total_start.elapsed();
    println!();
    println!("=== Summary ===");
    println!("Phase 1 (Loading):         {:?} ({:.1}%)", phase1_time, 100.0 * phase1_time.as_secs_f64() / total_time.as_secs_f64());
    println!("Phase 2 (Suffix Array):    {:?} ({:.1}%)", phase2_time, 100.0 * phase2_time.as_secs_f64() / total_time.as_secs_f64());
    println!("Phase 3 (Cache):           {:?} ({:.1}%)", phase3_time, 100.0 * phase3_time.as_secs_f64() / total_time.as_secs_f64());
    println!("Phase 4 (Sorted K-mers):   {:?} ({:.1}%)", phase4_time, 100.0 * phase4_time.as_secs_f64() / total_time.as_secs_f64());
    println!("Phase 5 (K-mer processing): {:?} ({:.1}%)", phase5_time, 100.0 * phase5_time.as_secs_f64() / total_time.as_secs_f64());
    println!("Total time:                {:?}", total_time);
    
    // Create dummy output
    std::fs::write(format!("{}.reprof", output_prefix), "# Detailed profiling run\n").ok();
    
    Ok(())
}