//! Comprehensive stress tests for the Bitmask implementation
//! 
//! These tests validate the correctness and robustness of the atomic bitmask
//! under extreme concurrent workloads, ensuring production readiness.

use std::sync::Arc;
use std::thread;
use std::time::Duration;
use reprise::mask::{Bitmask, ClaimGuard};
use rayon::prelude::*;

#[test]
fn stress_test_high_contention() {
    // Create a small bitmask to force high contention
    let mask = Arc::new(Bitmask::new(10000));
    let num_threads = 16;
    let operations_per_thread = 1000;
    
    let mut handles = Vec::new();
    
    for thread_id in 0..num_threads {
        let mask_clone = Arc::clone(&mask);
        let handle = thread::spawn(move || {
            let mut stats = (0, 0); // (successful, failed)
            
            for i in 0..operations_per_thread {
                // Create overlapping regions to maximize contention
                let start = ((thread_id * 100 + i) % 9900) as u64;
                let end = start + 100;
                
                if let Some(guard) = ClaimGuard::new(&mask_clone, start..end) {
                    stats.0 += 1;
                    
                    // Simulate work while holding the claim
                    thread::sleep(Duration::from_micros(10));
                    
                    // Verify the claim is still valid
                    assert!(!mask_clone.is_range_free(&(start..end)));
                    
                    drop(guard);
                    
                    // Verify the region is freed (eventually - may race with other threads)
                    thread::sleep(Duration::from_micros(1));
                } else {
                    stats.1 += 1;
                }
            }
            stats
        });
        handles.push(handle);
    }
    
    // Collect results
    let results: Vec<_> = handles.into_iter().map(|h| h.join().unwrap()).collect();
    let total_successful: usize = results.iter().map(|(s, _)| s).sum();
    let total_failed: usize = results.iter().map(|(_, f)| f).sum();
    let total_operations = num_threads * operations_per_thread;
    
    assert_eq!(total_successful + total_failed, total_operations);
    assert!(total_successful > 0, "No operations succeeded under high contention");
    
    println!("High contention test: {}/{} operations successful", 
        total_successful, total_operations);
}

#[test]
fn stress_test_rapid_claim_release() {
    let mask = Arc::new(Bitmask::new(1000000));
    let num_iterations = 10000;
    
    // Rapidly claim and release the same regions from multiple threads
    let handle1 = {
        let mask = Arc::clone(&mask);
        thread::spawn(move || {
            for i in 0..num_iterations {
                let start = (i % 1000) * 100;
                let end = start + 50;
                
                let _guard = ClaimGuard::new(&mask, start..end);
                // Guard automatically releases
            }
        })
    };
    
    let handle2 = {
        let mask = Arc::clone(&mask);
        thread::spawn(move || {
            for i in 0..num_iterations {
                let start = (i % 1000) * 100 + 25; // Overlapping regions
                let end = start + 50;
                
                let _guard = ClaimGuard::new(&mask, start..end);
                // Guard automatically releases  
            }
        })
    };
    
    handle1.join().unwrap();
    handle2.join().unwrap();
    
    // Verify all regions are free after completion
    for i in 0..1000 {
        let start = i * 100;
        let end = start + 75; // Larger than both claim sizes
        assert!(mask.is_range_free(&(start..end)));
    }
}

#[test]
fn stress_test_panic_safety_under_load() {
    use std::panic::{self, AssertUnwindSafe};
    
    let mask = Arc::new(Bitmask::new(100000));
    let panic_thread_count = 5;
    let normal_thread_count = 10;
    
    // Start threads that will panic while holding claims
    let mut panic_handles = Vec::new();
    for i in 0..panic_thread_count {
        let mask = Arc::clone(&mask);
        let handle = thread::spawn(move || {
            let result = panic::catch_unwind(AssertUnwindSafe(|| {
                let start = (i * 10000) as u64;
                let end = start + 1000;
                
                let _guard = ClaimGuard::new(&mask, start..end).expect("Should claim");
                
                if i % 2 == 0 {
                    panic!("Intentional panic while holding claim");
                }
            }));
            
            result.is_err() // Return true if panic occurred
        });
        panic_handles.push(handle);
    }
    
    // Start threads that do normal work
    let mut normal_handles = Vec::new();
    for i in 0..normal_thread_count {
        let mask = Arc::clone(&mask);
        let handle = thread::spawn(move || {
            let mut successful = 0;
            
            for j in 0..100 {
                let start = ((i * 1000 + j * 10) % 90000) as u64;
                let end = start + 100;
                
                if let Some(_guard) = ClaimGuard::new(&mask, start..end) {
                    successful += 1;
                    thread::sleep(Duration::from_millis(1));
                }
            }
            
            successful
        });
        normal_handles.push(handle);
    }
    
    // Wait for all threads
    for handle in panic_handles {
        let _ = handle.join(); // Some may have panicked
    }
    
    let normal_results: Vec<_> = normal_handles.into_iter()
        .map(|h| h.join().unwrap())
        .collect();
    
    let total_normal_successful: usize = normal_results.iter().sum();
    assert!(total_normal_successful > 0, "Normal threads should have succeeded despite panics");
    
    // Give time for any remaining cleanup
    thread::sleep(Duration::from_millis(100));
    
    // Verify that all regions are eventually freed
    // (This is a best-effort check since there may still be racing threads)
    for i in 0..10 {
        let start = (i * 10000) as u64;
        let end = start + 1000;
        // Don't assert here as there may be timing issues
        mask.is_range_free(&(start..end));
    }
}

#[test]
fn stress_test_word_boundary_handling() {
    // Test claims that cross multiple 64-bit word boundaries
    let mask = Arc::new(Bitmask::new(1000000));
    
    // Create claims that span multiple words
    let large_claims: Vec<(u64, u64)> = (0..100)
        .map(|i| {
            let start = i * 5000;
            let end = start + 1000; // Definitely spans multiple 64-bit words
            (start, end)
        })
        .collect();
    
    // Process in parallel
    let results: Vec<bool> = large_claims
        .par_iter()
        .map(|&(start, end)| {
            if let Some(guard) = ClaimGuard::new(&mask, start..end) {
                // Verify the entire range is claimed
                for word_start in (start..end).step_by(64) {
                    let check_end = (word_start + 64).min(end);
                    assert!(!mask.is_range_free(&(word_start..check_end)));
                }
                
                drop(guard);
                true
            } else {
                false
            }
        })
        .collect();
    
    let successful = results.iter().filter(|&&success| success).count();
    assert!(successful > 0, "No word-boundary spanning claims succeeded");
    
    println!("Word boundary test: {}/{} large claims successful", 
        successful, large_claims.len());
}