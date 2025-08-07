//! Thread-local object pooling for memory-efficient temporary buffer management
//!
//! This module provides thread-local object pools to reuse temporary buffers
//! and reduce allocation overhead during intensive genomic processing operations.

use std::cell::RefCell;
use std::collections::VecDeque;

/// Maximum number of objects to keep in each pool to prevent unbounded growth
const MAX_POOL_SIZE: usize = 64;

/// Thread-local pool for reusable Vec<u8> buffers
thread_local! {
    static BYTE_BUFFER_POOL: RefCell<VecDeque<Vec<u8>>> = RefCell::new(VecDeque::new());
    static POSITION_BUFFER_POOL: RefCell<VecDeque<Vec<u64>>> = RefCell::new(VecDeque::new());
    static KMER_BUFFER_POOL: RefCell<VecDeque<Vec<crate::kmer::Kmer>>> = RefCell::new(VecDeque::new());
}

/// RAII wrapper for pooled byte buffers that automatically returns to pool on drop
pub struct PooledByteBuffer {
    buffer: Option<Vec<u8>>,
}

impl PooledByteBuffer {
    /// Get a buffer from the thread-local pool, or create a new one
    pub fn new() -> Self {
        let buffer = BYTE_BUFFER_POOL.with(|pool| {
            pool.borrow_mut().pop_front()
        }).unwrap_or_else(Vec::new);
        
        Self {
            buffer: Some(buffer),
        }
    }
    
    /// Get a buffer with pre-allocated capacity
    pub fn with_capacity(capacity: usize) -> Self {
        let mut buffer = BYTE_BUFFER_POOL.with(|pool| {
            pool.borrow_mut().pop_front()
        }).unwrap_or_else(|| Vec::with_capacity(capacity));
        
        // Ensure capacity is at least what was requested
        if buffer.capacity() < capacity {
            buffer.reserve(capacity - buffer.len());
        }
        
        Self {
            buffer: Some(buffer),
        }
    }
    
    /// Get mutable reference to the underlying buffer
    pub fn as_mut(&mut self) -> &mut Vec<u8> {
        self.buffer.as_mut().expect("Buffer should always be present")
    }
    
    /// Get immutable reference to the underlying buffer
    pub fn as_ref(&self) -> &Vec<u8> {
        self.buffer.as_ref().expect("Buffer should always be present")
    }
    
    /// Clear the buffer while preserving capacity
    pub fn clear(&mut self) {
        self.as_mut().clear();
    }
    
    /// Get the current length of the buffer
    pub fn len(&self) -> usize {
        self.as_ref().len()
    }
    
    /// Check if the buffer is empty
    pub fn is_empty(&self) -> bool {
        self.as_ref().is_empty()
    }
    
    /// Get the capacity of the buffer
    pub fn capacity(&self) -> usize {
        self.as_ref().capacity()
    }
    
    /// Reserve additional capacity
    pub fn reserve(&mut self, additional: usize) {
        self.as_mut().reserve(additional);
    }
    
    /// Extend the buffer with data
    pub fn extend_from_slice(&mut self, data: &[u8]) {
        self.as_mut().extend_from_slice(data);
    }
    
    /// Push a single byte
    pub fn push(&mut self, byte: u8) {
        self.as_mut().push(byte);
    }
}

impl Drop for PooledByteBuffer {
    fn drop(&mut self) {
        if let Some(mut buffer) = self.buffer.take() {
            // Clear the buffer and return to pool if pool isn't full
            buffer.clear();
            
            BYTE_BUFFER_POOL.with(|pool| {
                let mut pool = pool.borrow_mut();
                if pool.len() < MAX_POOL_SIZE {
                    pool.push_back(buffer);
                }
                // Otherwise, let it drop to prevent unbounded growth
            });
        }
    }
}

/// RAII wrapper for pooled position buffers (Vec<u64>)
pub struct PooledPositionBuffer {
    buffer: Option<Vec<u64>>,
}

impl PooledPositionBuffer {
    /// Get a position buffer from the thread-local pool
    pub fn new() -> Self {
        let buffer = POSITION_BUFFER_POOL.with(|pool| {
            pool.borrow_mut().pop_front()
        }).unwrap_or_else(Vec::new);
        
        Self {
            buffer: Some(buffer),
        }
    }
    
    /// Get a position buffer with pre-allocated capacity
    pub fn with_capacity(capacity: usize) -> Self {
        let mut buffer = POSITION_BUFFER_POOL.with(|pool| {
            pool.borrow_mut().pop_front()
        }).unwrap_or_else(|| Vec::with_capacity(capacity));
        
        if buffer.capacity() < capacity {
            buffer.reserve(capacity - buffer.len());
        }
        
        Self {
            buffer: Some(buffer),
        }
    }
    
    /// Get mutable reference to the underlying buffer
    pub fn as_mut(&mut self) -> &mut Vec<u64> {
        self.buffer.as_mut().expect("Buffer should always be present")
    }
    
    /// Get immutable reference to the underlying buffer
    pub fn as_ref(&self) -> &Vec<u64> {
        self.buffer.as_ref().expect("Buffer should always be present")
    }
    
    /// Clear the buffer
    pub fn clear(&mut self) {
        self.as_mut().clear();
    }
    
    /// Get the length
    pub fn len(&self) -> usize {
        self.as_ref().len()
    }
    
    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.as_ref().is_empty()
    }
    
    /// Push a position
    pub fn push(&mut self, position: u64) {
        self.as_mut().push(position);
    }
    
    /// Extend with positions
    pub fn extend<I: IntoIterator<Item = u64>>(&mut self, iter: I) {
        self.as_mut().extend(iter);
    }
    
    /// Get the capacity
    pub fn capacity(&self) -> usize {
        self.as_ref().capacity()
    }
}

impl Drop for PooledPositionBuffer {
    fn drop(&mut self) {
        if let Some(mut buffer) = self.buffer.take() {
            buffer.clear();
            
            POSITION_BUFFER_POOL.with(|pool| {
                let mut pool = pool.borrow_mut();
                if pool.len() < MAX_POOL_SIZE {
                    pool.push_back(buffer);
                }
            });
        }
    }
}

/// RAII wrapper for pooled k-mer buffers
pub struct PooledKmerBuffer {
    buffer: Option<Vec<crate::kmer::Kmer>>,
}

impl PooledKmerBuffer {
    /// Get a k-mer buffer from the thread-local pool
    pub fn new() -> Self {
        let buffer = KMER_BUFFER_POOL.with(|pool| {
            pool.borrow_mut().pop_front()
        }).unwrap_or_else(Vec::new);
        
        Self {
            buffer: Some(buffer),
        }
    }
    
    /// Get a k-mer buffer with pre-allocated capacity
    pub fn with_capacity(capacity: usize) -> Self {
        let mut buffer = KMER_BUFFER_POOL.with(|pool| {
            pool.borrow_mut().pop_front()
        }).unwrap_or_else(|| Vec::with_capacity(capacity));
        
        if buffer.capacity() < capacity {
            buffer.reserve(capacity - buffer.len());
        }
        
        Self {
            buffer: Some(buffer),
        }
    }
    
    /// Get mutable reference to the underlying buffer
    pub fn as_mut(&mut self) -> &mut Vec<crate::kmer::Kmer> {
        self.buffer.as_mut().expect("Buffer should always be present")
    }
    
    /// Get immutable reference to the underlying buffer  
    pub fn as_ref(&self) -> &Vec<crate::kmer::Kmer> {
        self.buffer.as_ref().expect("Buffer should always be present")
    }
    
    /// Clear the buffer
    pub fn clear(&mut self) {
        self.as_mut().clear();
    }
    
    /// Get the length
    pub fn len(&self) -> usize {
        self.as_ref().len()
    }
    
    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.as_ref().is_empty()
    }
    
    /// Push a k-mer
    pub fn push(&mut self, kmer: crate::kmer::Kmer) {
        self.as_mut().push(kmer);
    }
    
    /// Extend with k-mers
    pub fn extend<I: IntoIterator<Item = crate::kmer::Kmer>>(&mut self, iter: I) {
        self.as_mut().extend(iter);
    }
}

impl Drop for PooledKmerBuffer {
    fn drop(&mut self) {
        if let Some(mut buffer) = self.buffer.take() {
            buffer.clear();
            
            KMER_BUFFER_POOL.with(|pool| {
                let mut pool = pool.borrow_mut();
                if pool.len() < MAX_POOL_SIZE {
                    pool.push_back(buffer);
                }
            });
        }
    }
}

/// Memory pool statistics for monitoring and debugging
#[derive(Debug, Clone)]
pub struct PoolStats {
    pub byte_buffers_pooled: usize,
    pub position_buffers_pooled: usize,
    pub kmer_buffers_pooled: usize,
}

impl PoolStats {
    /// Get current pool statistics for this thread
    pub fn current() -> Self {
        let byte_buffers_pooled = BYTE_BUFFER_POOL.with(|pool| pool.borrow().len());
        let position_buffers_pooled = POSITION_BUFFER_POOL.with(|pool| pool.borrow().len());
        let kmer_buffers_pooled = KMER_BUFFER_POOL.with(|pool| pool.borrow().len());
        
        Self {
            byte_buffers_pooled,
            position_buffers_pooled,
            kmer_buffers_pooled,
        }
    }
    
    /// Clear all thread-local pools (useful for testing or cleanup)
    pub fn clear_all_pools() {
        BYTE_BUFFER_POOL.with(|pool| pool.borrow_mut().clear());
        POSITION_BUFFER_POOL.with(|pool| pool.borrow_mut().clear());
        KMER_BUFFER_POOL.with(|pool| pool.borrow_mut().clear());
    }
}

/// Utility functions for efficient memory usage patterns
pub mod utils {
    use super::*;
    
    /// Create a byte buffer optimized for k-mer extraction
    /// 
    /// Pre-allocates based on typical k-mer extraction patterns
    pub fn kmer_extraction_buffer(sequence_len: usize, k: usize) -> PooledByteBuffer {
        // Estimate buffer size needed for k-mer extraction operations
        let estimated_size = (sequence_len.saturating_sub(k - 1)) * 4; // Conservative estimate
        PooledByteBuffer::with_capacity(estimated_size.min(64 * 1024)) // Cap at 64KB
    }
    
    /// Create a position buffer for genomic ranges
    pub fn position_buffer(estimated_positions: usize) -> PooledPositionBuffer {
        PooledPositionBuffer::with_capacity(estimated_positions)
    }
    
    /// Create a k-mer buffer for sequence processing
    pub fn kmer_buffer(estimated_kmers: usize) -> PooledKmerBuffer {
        PooledKmerBuffer::with_capacity(estimated_kmers)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_byte_buffer_pooling() {
        // Clear pools to start fresh
        PoolStats::clear_all_pools();
        
        let initial_stats = PoolStats::current();
        assert_eq!(initial_stats.byte_buffers_pooled, 0);
        
        // Create and use a buffer
        {
            let mut buffer = PooledByteBuffer::with_capacity(1024);
            buffer.push(42);
            buffer.extend_from_slice(&[1, 2, 3, 4]);
            assert_eq!(buffer.len(), 5);
            assert!(buffer.capacity() >= 1024);
        } // Buffer should return to pool here
        
        let after_stats = PoolStats::current();
        assert_eq!(after_stats.byte_buffers_pooled, 1);
        
        // Reuse the buffer
        {
            let buffer = PooledByteBuffer::new();
            assert!(buffer.capacity() >= 1024); // Should reuse the same buffer
            assert_eq!(buffer.len(), 0); // Should be cleared
        }
    }
    
    #[test]
    fn test_position_buffer_pooling() {
        PoolStats::clear_all_pools();
        
        {
            let mut buffer = PooledPositionBuffer::with_capacity(100);
            buffer.push(12345);
            buffer.extend([67890, 11111]);
            assert_eq!(buffer.len(), 3);
        }
        
        let stats = PoolStats::current();
        assert_eq!(stats.position_buffers_pooled, 1);
    }
    
    #[test]
    fn test_pool_size_limits() {
        PoolStats::clear_all_pools();
        
        // Create more buffers than the pool limit
        let mut buffers = Vec::new();
        for _ in 0..MAX_POOL_SIZE + 10 {
            buffers.push(PooledByteBuffer::new());
        }
        
        // Drop all buffers
        drop(buffers);
        
        let stats = PoolStats::current();
        // Should not exceed MAX_POOL_SIZE
        assert!(stats.byte_buffers_pooled <= MAX_POOL_SIZE);
    }
    
    #[test]
    fn test_utility_functions() {
        PoolStats::clear_all_pools();
        
        let buffer = utils::kmer_extraction_buffer(1000, 13);
        assert!(buffer.capacity() > 0);
        
        let pos_buffer = utils::position_buffer(500);
        assert!(pos_buffer.capacity() >= 500);
    }
}