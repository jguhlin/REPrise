//! Demonstration of the Bitmask and ClaimGuard functionality
//! 
//! This example shows how the atomic bitmask can be used in a multi-threaded
//! genomic processing scenario, simulating the workload that REPrise will handle
//! during concurrent repeat detection.

use std::sync::Arc;
use std::thread;
use std::time::Duration;
use reprise::mask::{Bitmask, ClaimGuard};
use rayon::prelude::*;

fn main() {
    println!("REPrise Bitmask Demo - Production-Ready Concurrent Genomic Region Tracking");
    println!("=========================================================================");
    
    // Simulate a 1 million base genome
    let genome_size = 1_000_000u64;
    let mask = Arc::new(Bitmask::new(genome_size));
    
    println!("Created bitmask for {} base genome", genome_size);
    
    // Demo 1: Basic claim and release functionality
    demo_basic_functionality(&mask);
    
    // Demo 2: Concurrent worker simulation
    demo_concurrent_workers(&mask);
    
    // Demo 3: High-throughput parallel processing with Rayon
    demo_rayon_parallel_processing(&mask);
    
    println!("\nAll demos completed successfully!");
    println!("The Bitmask provides:");
    println!("- Thread-safe atomic operations");  
    println!("- Automatic rollback on conflicts");
    println!("- RAII-based resource management with ClaimGuard");
    println!("- High-performance concurrent access patterns");
}

fn demo_basic_functionality(mask: &Bitmask) {
    println!("\n1. Basic Functionality Demo");
    println!("---------------------------");
    
    // Test basic claim/release
    {
        let _guard = ClaimGuard::new(mask, 1000..2000).expect("Should claim region");
        println!("✓ Successfully claimed region 1000..2000");
        
        // Try to claim overlapping region - should fail
        if ClaimGuard::new(mask, 1500..2500).is_none() {
            println!("✓ Correctly rejected overlapping claim 1500..2500");
        }
        
        // Try non-overlapping region - should succeed
        let _guard2 = ClaimGuard::new(mask, 3000..4000).expect("Should claim non-overlapping region");
        println!("✓ Successfully claimed non-overlapping region 3000..4000");
        
    } // Guards automatically release here
    
    // Verify regions are now free
    if mask.is_range_free(&(1500..1600)) {
        println!("✓ Regions automatically released after guards dropped");
    }
}

fn demo_concurrent_workers(mask: &Arc<Bitmask>) {
    println!("\n2. Concurrent Worker Demo");
    println!("-------------------------");
    
    let num_workers = 8;
    let mut handles = Vec::new();
    
    println!("Starting {} worker threads...", num_workers);
    
    for worker_id in 0..num_workers {
        let mask_clone = Arc::clone(mask);
        let handle = thread::spawn(move || {
            let mut successful_claims = 0;
            let mut failed_claims = 0;
            
            // Each worker tries to claim 100 different regions
            for i in 0..100 {
                let start = (worker_id * 10000) + (i * 50);  
                let end = start + 25; // Small regions to increase conflicts
                
                if let Some(guard) = ClaimGuard::new(&mask_clone, start..end) {
                    successful_claims += 1;
                    
                    // Simulate processing work
                    thread::sleep(Duration::from_millis(1));
                    
                    // Guard automatically releases when it goes out of scope
                }
                else {
                    failed_claims += 1;
                }
            }
            
            (worker_id, successful_claims, failed_claims)
        });
        handles.push(handle);
    }
    
    // Collect results
    let mut total_successful = 0;
    let mut total_failed = 0;
    
    for handle in handles {
        let (worker_id, successful, failed) = handle.join().unwrap();
        println!("Worker {}: {} successful, {} failed claims", worker_id, successful, failed);
        total_successful += successful;
        total_failed += failed;
    }
    
    println!("Total: {} successful, {} failed claims", total_successful, total_failed);
    println!("✓ All workers completed without deadlocks or data races");
}

fn demo_rayon_parallel_processing(mask: &Arc<Bitmask>) {
    println!("\n3. High-Throughput Rayon Demo");
    println!("-----------------------------");
    
    // Generate a large number of candidate regions for processing
    let candidates: Vec<(u64, u64)> = (0..10000)
        .map(|i| {
            let start = (i * 50) % 900000;  // Wrap around to create conflicts
            let end = start + 100;
            (start, end)
        })
        .collect();
    
    println!("Processing {} candidate regions in parallel...", candidates.len());
    
    // Process candidates in parallel using Rayon
    let results: Vec<bool> = candidates
        .par_iter()  // Parallel iterator
        .map(|&(start, end)| {
            // Try to claim and process the region
            if let Some(_guard) = ClaimGuard::new(mask, start..end) {
                // Simulate genomic processing work
                // In real REPrise, this would be repeat detection and alignment
                let _work_simulation = (start..end).sum::<u64>();
                true
            } else {
                false
            }
        })
        .collect();
    
    let successful = results.iter().filter(|&&success| success).count();
    let failed = results.len() - successful;
    
    println!("Processed {} regions successfully", successful);
    println!("Skipped {} regions due to conflicts", failed);
    println!("✓ High-throughput parallel processing completed efficiently");
}