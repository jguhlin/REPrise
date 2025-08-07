// In src/mask.rs
use std::sync::atomic::{AtomicU64, Ordering};
use std::ops::Range;

// A memory-efficient, thread-safe bitmask for tracking claimed genomic regions.
pub struct Bitmask {
    data: Vec<AtomicU64>,
    len_bases: u64,
}

impl Bitmask {
    const BITS_PER_WORD: u64 = 64;

    pub fn new(len_bases: u64) -> Self {
        let num_words = (len_bases + Self::BITS_PER_WORD - 1) / Self::BITS_PER_WORD;
        let data = (0..num_words).map(|_| AtomicU64::new(0)).collect();
        Self { data, len_bases }
    }

    /// Attempts to atomically claim all bits in a range.
    /// This is an all-or-nothing operation with internal rollback.
    pub fn claim_range(&self, range: &Range<u64>) -> bool {
        if range.start >= range.end || range.end > self.len_bases { return false; }

        let mut claimed_words: Vec<(usize, u64)> = Vec::new();
        for word_idx in (range.start / Self::BITS_PER_WORD)..=((range.end - 1) / Self::BITS_PER_WORD) {
            let word_start = word_idx * Self::BITS_PER_WORD;
            let word_end = word_start + Self::BITS_PER_WORD;

            let claim_start = range.start.max(word_start);
            let claim_end = range.end.min(word_end);

            let mut bitmask_for_word = 0u64;
            for i in claim_start..claim_end {
                bitmask_for_word |= 1u64 << (i % Self::BITS_PER_WORD);
            }

            let atomic_word = &self.data[word_idx as usize];
            let mut current_val = atomic_word.load(Ordering::Acquire);
            loop {
                // Check for conflict inside the loop to avoid race conditions
                if (current_val & bitmask_for_word) != 0 {
                    // Conflict detected. Roll back all previous claims.
                    for (idx, mask) in &claimed_words {
                        self.data[*idx].fetch_and(!mask, Ordering::Release);
                    }
                    return false;
                }
                
                match atomic_word.compare_exchange_weak(
                    current_val, current_val | bitmask_for_word,
                    Ordering::AcqRel, Ordering::Acquire,
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

    /// Releases all bits in a range. Idempotent.
    pub fn release_range(&self, range: &Range<u64>) {
        if range.start >= range.end || range.end > self.len_bases { return; }
        
        for word_idx in (range.start / Self::BITS_PER_WORD)..=((range.end - 1) / Self::BITS_PER_WORD) {
            let word_start = word_idx * Self::BITS_PER_WORD;
            let claim_start = range.start.max(word_start);
            let claim_end = range.end.min(word_start + Self::BITS_PER_WORD);
            
            let mut bitmask = 0u64;
            for i in claim_start..claim_end {
                bitmask |= 1u64 << (i % Self::BITS_PER_WORD);
            }
            
            self.data[word_idx as usize].fetch_and(!bitmask, Ordering::Release);
        }
    }

    /// A cheap, non-atomic check for pre-filtering.
    pub fn is_range_free(&self, range: &Range<u64>) -> bool {
        if range.start >= range.end || range.end > self.len_bases { return false; }
        
        for word_idx in (range.start / Self::BITS_PER_WORD)..=((range.end - 1) / Self::BITS_PER_WORD) {
            let word_start = word_idx * Self::BITS_PER_WORD;
            let claim_start = range.start.max(word_start);
            let claim_end = range.end.min(word_start + Self::BITS_PER_WORD);
            
            let mut bitmask = 0u64;
            for i in claim_start..claim_end {
                bitmask |= 1u64 << (i % Self::BITS_PER_WORD);
            }
            
            if (self.data[word_idx as usize].load(Ordering::Relaxed) & bitmask) != 0 {
                return false;
            }
        }
        true
    }
}

/// A panic-safe RAII guard for atomically claiming a region.
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
    pub fn new(mask: &'a Bitmask, range: Range<u64>) -> Option<Self> {
        if mask.claim_range(&range) {
            Some(Self { mask, claimed_range: Some(range) })
        } else {
            None
        }
    }
}