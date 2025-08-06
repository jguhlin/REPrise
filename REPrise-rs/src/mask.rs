//! Atomic bitmask for thread-safe genomic region tracking
//! 
//! This module provides the core concurrency primitives for safe parallel repeat detection.
//! The `Bitmask` uses atomic operations to enable lock-free claiming and releasing of genomic
//! regions, while `ClaimGuard` provides RAII-based automatic cleanup with panic safety.

use std::sync::atomic::{AtomicU64, Ordering};
use std::ops::Range;

/// A memory-efficient, thread-safe bitmask for tracking claimed genomic regions.
/// 
/// The `Bitmask` enables multiple threads to atomically claim non-overlapping regions
/// of a genome for repeat detection processing. It uses 64-bit atomic integers for
/// efficient concurrent access and provides rollback capabilities for failed claims.
pub struct Bitmask {
    data: Vec<AtomicU64>,
    len_bases: u64,
}

impl Bitmask {
    const BITS_PER_WORD: u64 = 64;

    /// Creates a new bitmask capable of tracking `len_bases` genomic positions.
    /// 
    /// # Arguments
    /// * `len_bases` - The total number of genomic positions to track
    /// 
    /// # Returns
    /// A new `Bitmask` with all positions initially unclaimed
    pub fn new(len_bases: u64) -> Self {
        let num_words = (len_bases + Self::BITS_PER_WORD - 1) / Self::BITS_PER_WORD;
        let data = (0..num_words).map(|_| AtomicU64::new(0)).collect();
        Self { data, len_bases }
    }

    /// Attempts to atomically claim all bits in a range.
    /// This is an all-or-nothing operation with internal rollback on conflict.
    /// 
    /// # Arguments
    /// * `range` - The genomic range to claim (half-open interval [start, end))
    /// 
    /// # Returns
    /// * `true` if the entire range was successfully claimed
    /// * `false` if any part of the range was already claimed, or if the range is invalid
    /// 
    /// # Behavior
    /// If a conflict is detected partway through claiming a range, all previously
    /// claimed parts are automatically rolled back, ensuring atomic semantics.
    pub fn claim_range(&self, range: &Range<u64>) -> bool {
        if range.start >= range.end || range.end > self.len_bases { 
            return false; 
        }

        let mut claimed_words: Vec<(usize, u64)> = Vec::new();
        
        // Process each 64-bit word that overlaps with the range
        for word_idx in (range.start / Self::BITS_PER_WORD)..=((range.end - 1) / Self::BITS_PER_WORD) {
            let word_start = word_idx * Self::BITS_PER_WORD;
            let word_end = word_start + Self::BITS_PER_WORD;

            let claim_start = range.start.max(word_start);
            let claim_end = range.end.min(word_end);

            // Build bitmask for the portion of this word we want to claim
            let mut bitmask_for_word = 0u64;
            for i in claim_start..claim_end {
                bitmask_for_word |= 1u64 << (i % Self::BITS_PER_WORD);
            }

            let atomic_word = &self.data[word_idx as usize];
            let mut current_val = atomic_word.load(Ordering::Acquire);
            
            // Retry loop for compare-and-swap
            loop {
                // Check for conflict inside the loop to avoid race conditions
                if (current_val & bitmask_for_word) != 0 {
                    // Conflict detected. Roll back all previous claims.
                    for &(idx, mask) in &claimed_words {
                        self.data[idx].fetch_and(!mask, Ordering::Release);
                    }
                    return false;
                }
                
                match atomic_word.compare_exchange_weak(
                    current_val, 
                    current_val | bitmask_for_word,
                    Ordering::AcqRel, 
                    Ordering::Acquire,
                ) {
                    Ok(_) => {
                        claimed_words.push((word_idx as usize, bitmask_for_word));
                        break;
                    }
                    Err(new_current_val) => current_val = new_current_val,
                }
            }
        }
        true
    }

    /// Releases all bits in a range. This operation is idempotent.
    /// 
    /// # Arguments
    /// * `range` - The genomic range to release (half-open interval [start, end))
    /// 
    /// # Behavior
    /// Clearing already-clear bits is safe and has no effect. Invalid ranges are ignored.
    pub fn release_range(&self, range: &Range<u64>) {
        if range.start >= range.end || range.end > self.len_bases { 
            return; 
        }
        
        // Process each 64-bit word that overlaps with the range
        for word_idx in (range.start / Self::BITS_PER_WORD)..=((range.end - 1) / Self::BITS_PER_WORD) {
            let word_start = word_idx * Self::BITS_PER_WORD;
            let claim_start = range.start.max(word_start);
            let claim_end = range.end.min(word_start + Self::BITS_PER_WORD);
            
            // Build bitmask for the portion of this word we want to clear
            let mut bitmask = 0u64;
            for i in claim_start..claim_end {
                bitmask |= 1u64 << (i % Self::BITS_PER_WORD);
            }
            
            // Atomically clear the bits using fetch_and with inverted mask
            self.data[word_idx as usize].fetch_and(!bitmask, Ordering::Release);
        }
    }

    /// Performs a cheap, non-atomic check to see if a range appears to be free.
    /// This is intended for pre-filtering to avoid expensive atomic operations.
    /// 
    /// # Arguments
    /// * `range` - The genomic range to check (half-open interval [start, end))
    /// 
    /// # Returns
    /// * `true` if the range appears to be free (but may become claimed before use)
    /// * `false` if any part of the range is definitely claimed, or if the range is invalid
    /// 
    /// # Note
    /// This uses `Ordering::Relaxed` for performance. The result may be stale by the time
    /// it's used, so callers should still be prepared to handle failed claims.
    pub fn is_range_free(&self, range: &Range<u64>) -> bool {
        if range.start >= range.end || range.end > self.len_bases { 
            return false; 
        }
        
        // Check each 64-bit word that overlaps with the range
        for word_idx in (range.start / Self::BITS_PER_WORD)..=((range.end - 1) / Self::BITS_PER_WORD) {
            let word_start = word_idx * Self::BITS_PER_WORD;
            let claim_start = range.start.max(word_start);
            let claim_end = range.end.min(word_start + Self::BITS_PER_WORD);
            
            // Build bitmask for the portion of this word we're checking
            let mut bitmask = 0u64;
            for i in claim_start..claim_end {
                bitmask |= 1u64 << (i % Self::BITS_PER_WORD);
            }
            
            // Use relaxed ordering for this non-critical check
            if (self.data[word_idx as usize].load(Ordering::Relaxed) & bitmask) != 0 {
                return false;
            }
        }
        true
    }

    /// Returns the total number of bases this bitmask can track.
    pub fn len(&self) -> u64 {
        self.len_bases
    }

    /// Returns true if this bitmask tracks zero bases.
    pub fn is_empty(&self) -> bool {
        self.len_bases == 0
    }
}

/// A panic-safe RAII guard for atomically claiming a genomic region.
/// 
/// The `ClaimGuard` ensures that claimed regions are automatically released
/// when the guard goes out of scope, even in the presence of panics. This
/// prevents deadlocks and resource leaks in concurrent repeat detection.
#[must_use = "Claim must be held by a variable to ensure its lifetime"]
pub struct ClaimGuard<'a> {
    mask: &'a Bitmask,
    claimed_range: Option<Range<u64>>,
}

impl<'a> Drop for ClaimGuard<'a> {
    fn drop(&mut self) {
        if let Some(range) = self.claimed_range.take() {
            self.mask.release_range(&range);
        }
    }
}

impl<'a> ClaimGuard<'a> {
    /// Attempts to create a new claim guard by atomically claiming the specified range.
    /// 
    /// # Arguments
    /// * `mask` - The bitmask to claim from
    /// * `range` - The genomic range to claim
    /// 
    /// # Returns
    /// * `Some(ClaimGuard)` if the range was successfully claimed
    /// * `None` if the range could not be claimed (already in use or invalid)
    /// 
    /// # Behavior
    /// The returned guard will automatically release the claim when dropped,
    /// providing exception-safe resource management.
    pub fn new(mask: &'a Bitmask, range: Range<u64>) -> Option<Self> {
        if mask.claim_range(&range) {
            Some(Self { 
                mask, 
                claimed_range: Some(range) 
            })
        } else {
            None
        }
    }

    /// Returns the range that this guard has claimed, if any.
    pub fn range(&self) -> Option<&Range<u64>> {
        self.claimed_range.as_ref()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use std::thread;

    #[test]
    fn test_bitmask_basic_operations() {
        let mask = Bitmask::new(1000);
        
        // Test successful claim
        assert!(mask.claim_range(&(100..200)));
        assert!(!mask.is_range_free(&(150..160)));
        
        // Test overlapping claim fails
        assert!(!mask.claim_range(&(150..250)));
        
        // Test non-overlapping claim succeeds
        assert!(mask.claim_range(&(300..400)));
        
        // Test release
        mask.release_range(&(100..200));
        assert!(mask.is_range_free(&(150..160)));
        
        // Test can reclaim after release
        assert!(mask.claim_range(&(100..200)));
    }

    #[test]
    fn test_bitmask_edge_cases() {
        let mask = Bitmask::new(100);
        
        // Test empty range
        assert!(!mask.claim_range(&(50..50)));
        
        // Test invalid range (start >= end)
        assert!(!mask.claim_range(&(60..50)));
        
        // Test out of bounds
        assert!(!mask.claim_range(&(50..150)));
        
        // Test single position
        assert!(mask.claim_range(&(50..51)));
        assert!(!mask.is_range_free(&(50..51)));
        
        // Test full range
        let small_mask = Bitmask::new(64);
        assert!(small_mask.claim_range(&(0..64)));
        assert!(!small_mask.is_range_free(&(0..1)));
    }

    #[test]
    fn test_word_boundaries() {
        let mask = Bitmask::new(200);
        
        // Test claim spanning word boundary (around bit 64)
        assert!(mask.claim_range(&(60..68)));
        assert!(!mask.is_range_free(&(63..65)));
        
        // Test claim exactly at word boundary
        assert!(mask.claim_range(&(128..132)));
        assert!(!mask.is_range_free(&(128..129)));
    }

    #[test]
    fn test_claim_guard_basic() {
        let mask = Bitmask::new(1000);
        
        {
            let _guard = ClaimGuard::new(&mask, 100..200).expect("Should claim successfully");
            assert!(!mask.is_range_free(&(150..160)));
        } // guard drops here
        
        // Range should be free after guard drops
        assert!(mask.is_range_free(&(150..160)));
    }

    #[test]
    fn test_claim_guard_failure() {
        let mask = Bitmask::new(1000);
        
        let _guard1 = ClaimGuard::new(&mask, 100..200).expect("Should claim successfully");
        let guard2 = ClaimGuard::new(&mask, 150..250);
        
        assert!(guard2.is_none(), "Overlapping claim should fail");
    }

    #[test]
    fn test_concurrent_claims() {
        let mask = Arc::new(Bitmask::new(10000));
        let mut handles = Vec::new();
        
        // Spawn multiple threads trying to claim different ranges
        for i in 0..10 {
            let mask_clone = Arc::clone(&mask);
            let handle = thread::spawn(move || {
                let start = i * 1000;
                let end = start + 500;
                mask_clone.claim_range(&(start..end))
            });
            handles.push(handle);
        }
        
        // All threads should succeed (non-overlapping ranges)
        let results: Vec<_> = handles.into_iter().map(|h| h.join().unwrap()).collect();
        assert!(results.iter().all(|&success| success));
    }

    #[test]
    fn test_concurrent_conflicting_claims() {
        let mask = Arc::new(Bitmask::new(1000));
        let mut handles = Vec::new();
        
        // Spawn multiple threads trying to claim the same range
        for _ in 0..10 {
            let mask_clone = Arc::clone(&mask);
            let handle = thread::spawn(move || {
                mask_clone.claim_range(&(100..200))
            });
            handles.push(handle);
        }
        
        // Exactly one thread should succeed
        let results: Vec<_> = handles.into_iter().map(|h| h.join().unwrap()).collect();
        let success_count = results.iter().filter(|&&success| success).count();
        assert_eq!(success_count, 1);
    }

    #[test]
    fn test_panic_safety() {
        use std::panic;
        
        let mask = Arc::new(Bitmask::new(1000));
        let mask_clone = Arc::clone(&mask);
        
        let result = panic::catch_unwind(move || {
            let _guard = ClaimGuard::new(&mask_clone, 100..200).expect("Should claim");
            panic!("Simulated panic");
        });
        
        assert!(result.is_err());
        
        // Range should be free after panic (guard should have been dropped)
        assert!(mask.is_range_free(&(150..160)));
        
        // Should be able to claim again
        assert!(mask.claim_range(&(100..200)));
    }

    #[test]
    fn test_rollback_on_partial_failure() {
        let mask = Bitmask::new(200);
        
        // Claim a range that will conflict with our multi-word claim
        assert!(mask.claim_range(&(120..130)));
        
        // Try to claim a range spanning multiple words that includes the conflict
        // This should fail and roll back any partial progress
        assert!(!mask.claim_range(&(60..140)));
        
        // Verify that the non-conflicting parts weren't claimed
        assert!(mask.is_range_free(&(60..120)));
        assert!(mask.is_range_free(&(130..140)));
        
        // Original claim should still be intact
        assert!(!mask.is_range_free(&(125..126)));
    }
}

/// Loom-based stress tests for concurrency validation
#[cfg(test)]
#[cfg(loom)]
mod loom_tests {
    use super::*;
    use loom::sync::Arc;
    use loom::thread;

    #[test]
    fn loom_concurrent_claims() {
        loom::model(|| {
            let mask = Arc::new(Bitmask::new(100));
            
            let mask1 = mask.clone();
            let mask2 = mask.clone();
            
            let h1 = thread::spawn(move || {
                mask1.claim_range(&(10..20))
            });
            
            let h2 = thread::spawn(move || {
                mask2.claim_range(&(10..20))
            });
            
            let r1 = h1.join().unwrap();
            let r2 = h2.join().unwrap();
            
            // Exactly one should succeed
            assert_ne!(r1, r2);
        });
    }

    #[test]
    fn loom_claim_guard_drop() {
        loom::model(|| {
            let mask = Arc::new(Bitmask::new(100));
            
            let mask1 = mask.clone();
            let mask2 = mask.clone();
            
            let h1 = thread::spawn(move || {
                let _guard = ClaimGuard::new(&mask1, 10..20);
                // Guard drops at end of scope
            });
            
            let h2 = thread::spawn(move || {
                // This may succeed or fail depending on timing
                ClaimGuard::new(&mask2, 10..20).is_some()
            });
            
            h1.join().unwrap();
            let _result = h2.join().unwrap();
            
            // After both threads finish, range should be free
            assert!(mask.is_range_free(&(10..20)));
        });
    }
}