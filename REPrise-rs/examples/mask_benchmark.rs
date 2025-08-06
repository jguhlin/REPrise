//! Performance benchmark for the Bitmask implementation
//! 
//! This benchmark measures the performance characteristics of the atomic bitmask
//! under various workloads, demonstrating its readiness for production use with
//! large genomes.

use std::sync::Arc;
use std::time::Instant;
use reprise::mask::{Bitmask, ClaimGuard};
use rayon::prelude::*;

fn main() {
    println!("REPrise Bitmask Performance Benchmark");
    println!("====================================");
    
    // Test with genome sizes typical for REPrise workloads
    let genome_sizes = vec![
        ("Human genome (3.2 Gb)", 3_200_000_000u64),
        ("Large plant genome (12 Gb)", 12_000_000_000u64),
        ("Fragmented assembly (1M contigs)", 1_000_000_000u64),
    ];
    
    for (name, size) in genome_sizes {
        println!("\n{}", name);
        println!("{}", "=".repeat(name.len()));
        
        benchmark_genome_size(size);
    }
}

fn benchmark_genome_size(genome_size: u64) {
    let mask = Arc::new(Bitmask::new(genome_size));
    
    // Benchmark 1: Sequential claim/release performance
    println!("1. Sequential operations:");
    benchmark_sequential_operations(&mask, 10000);
    
    // Benchmark 2: Concurrent claim attempts
    println!("2. Concurrent operations:");
    benchmark_concurrent_operations(&mask);
    
    // Benchmark 3: Memory usage
    println!("3. Memory efficiency:");
    let memory_mb = (genome_size + 63) / 64 * 8 / (1024 * 1024);
    println!("   Memory usage: {} MB ({:.3}% of genome size)", 
        memory_mb, memory_mb as f64 * 1024.0 * 1024.0 / genome_size as f64 * 100.0);
}

fn benchmark_sequential_operations(mask: &Arc<Bitmask>, num_operations: usize) {
    let start = Instant::now();
    
    for i in 0..num_operations {
        let region_start = (i as u64 * 1000) % (mask.len() - 1000);
        let region_end = region_start + 500;
        
        if let Some(guard) = ClaimGuard::new(mask, region_start..region_end) {
            // Simulate minimal processing
            std::hint::black_box(guard.range());
        }
    }
    
    let elapsed = start.elapsed();
    let ops_per_sec = num_operations as f64 / elapsed.as_secs_f64();
    
    println!("   {} operations in {:.2?}", num_operations, elapsed);
    println!("   {:.0} operations/second", ops_per_sec);
}

fn benchmark_concurrent_operations(mask: &Arc<Bitmask>) {
    let num_threads = rayon::current_num_threads();
    let operations_per_thread = 1000;
    
    println!("   Using {} threads", num_threads);
    
    let start = Instant::now();
    
    // Create overlapping work to test contention handling
    let results: Vec<usize> = (0..num_threads)
        .into_par_iter()
        .map(|thread_id| {
            let mut successful_claims = 0;
            
            for i in 0..operations_per_thread {
                // Create regions that will sometimes overlap between threads
                let region_start = ((thread_id * 500 + i * 100) as u64) % (mask.len() - 1000);
                let region_end = region_start + 200;
                
                if let Some(_guard) = ClaimGuard::new(mask, region_start..region_end) {
                    successful_claims += 1;
                    // Simulate processing work
                    std::hint::black_box(region_start + region_end);
                }
            }
            
            successful_claims
        })
        .collect();
    
    let elapsed = start.elapsed();
    let total_operations = num_threads * operations_per_thread;
    let successful_operations: usize = results.iter().sum();
    let ops_per_sec = total_operations as f64 / elapsed.as_secs_f64();
    
    println!("   {} total operations in {:.2?}", total_operations, elapsed);
    println!("   {} successful ({:.1}% success rate)", 
        successful_operations, 
        successful_operations as f64 / total_operations as f64 * 100.0);
    println!("   {:.0} operations/second", ops_per_sec);
}

#[cfg(test)]
mod benchmarks {
    use super::*;
    use std::time::Duration;
    
    #[test]
    fn benchmark_large_genome_performance() {
        // Test that operations complete in reasonable time even for large genomes
        let large_mask = Arc::new(Bitmask::new(1_000_000_000)); // 1 Gb
        
        let start = Instant::now();
        
        // Try 1000 claims across the genome
        for i in 0..1000 {
            let region_start = (i as u64 * 1_000_000) % (large_mask.len() - 10000);
            let region_end = region_start + 5000;
            
            let _guard = ClaimGuard::new(&large_mask, region_start..region_end);
        }
        
        let elapsed = start.elapsed();
        
        // Should complete well under 1 second even for large genomes
        assert!(elapsed < Duration::from_secs(1), 
            "Large genome operations took too long: {:?}", elapsed);
        
        println!("Large genome (1Gb) test completed in {:?}", elapsed);
    }
}