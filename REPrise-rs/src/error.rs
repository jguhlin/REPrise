//! Error handling for REPrise scaling architecture
//! 
//! This module provides comprehensive error types for all REPrise operations,
//! enabling clear error reporting and proper error handling throughout the
//! concurrent genomic processing pipeline.

use thiserror::Error;
use crossbeam::channel::SendError;
use rayon::ThreadPoolBuildError;
use crate::pipeline::CandidatePair;

/// Comprehensive error type for all REPrise operations
#[derive(Error, Debug)]
pub enum REPriseError {
    /// I/O errors (file operations, network, etc.)
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    
    /// Invalid k-mer containing non-ACGT bases
    #[error("Invalid k-mer containing non-ACGT base at position {0}")]
    InvalidKmer(u64),
    
    /// Alignment attempted across contig boundary
    #[error("Alignment attempted across contig boundary")]
    BoundaryViolation,
    
    /// External tool execution failed
    #[error("External tool '{0}' failed: {1}")]
    ExternalTool(String, String),
    
    /// Invalid FASTA format
    #[error("Invalid FASTA format at line {line}: {message}")]
    InvalidFasta { line: usize, message: String },
    
    /// Genome too large for current architecture
    #[error("Genome size {0} exceeds maximum supported size {1}")]
    GenomeTooLarge(u64, u64),
    
    /// Invalid k-mer length
    #[error("K-mer length {k} is invalid (must be between {min} and {max})")]
    InvalidKmerLength { k: usize, min: usize, max: usize },
    
    /// Contig not found
    #[error("Contig with ID {0} not found")]
    ContigNotFound(u32),
    
    /// Invalid genomic range
    #[error("Invalid genomic range {start}..{end} (start >= end or out of bounds)")]
    InvalidRange { start: u64, end: u64 },
    
    /// Memory allocation failed
    #[error("Memory allocation failed: requested {0} bytes")]
    OutOfMemory(u64),
    
    /// Concurrency error (e.g., channel closed unexpectedly)
    #[error("Concurrency error: {0}")]
    Concurrency(String),
    
    /// Configuration error
    #[error("Configuration error: {0}")]
    Config(String),
    
    /// Parse error for numeric or other structured data
    #[error("Parse error: {0}")]
    Parse(String),

    #[error("Send error: {0}")]
    SendError(#[from] SendError<CandidatePair>),

    #[error("Thread pool build error: {0}")]
    ThreadPoolBuildError(#[from] ThreadPoolBuildError),
    
    /// Thread join error
    #[error("Thread error: {0}")]
    ThreadError(String),
    
    /// Needletail parser error
    #[error("FASTA/Q parsing error: {0}")]
    NeedletailError(#[from] needletail::errors::ParseError),
    
    /// Generic anyhow error for complex nested errors
    #[error(transparent)]
    Other(#[from] anyhow::Error),
}

impl REPriseError {
    /// Create an InvalidFasta error with line number and message
    pub fn invalid_fasta(line: usize, message: impl Into<String>) -> Self {
        Self::InvalidFasta {
            line,
            message: message.into(),
        }
    }
    
    /// Create an InvalidKmerLength error with current and valid ranges
    pub fn invalid_kmer_length(k: usize, min: usize, max: usize) -> Self {
        Self::InvalidKmerLength { k, min, max }
    }
    
    /// Create an InvalidRange error
    pub fn invalid_range(start: u64, end: u64) -> Self {
        Self::InvalidRange { start, end }
    }
    
    /// Create a Concurrency error
    pub fn concurrency(message: impl Into<String>) -> Self {
        Self::Concurrency(message.into())
    }
    
    /// Create a Config error
    pub fn config(message: impl Into<String>) -> Self {
        Self::Config(message.into())
    }
    
    /// Create a Parse error
    pub fn parse(message: impl Into<String>) -> Self {
        Self::Parse(message.into())
    }
    
    /// Create an I/O error
    pub fn io_error(message: impl Into<String>) -> Self {
        Self::Io(std::io::Error::new(std::io::ErrorKind::Other, message.into()))
    }
}

/// Result type alias for REPrise operations
pub type Result<T> = std::result::Result<T, REPriseError>;

#[cfg(test)]
mod tests {
    use super::*;
    use std::io;
    
    #[test]
    fn test_error_display() {
        let err = REPriseError::InvalidKmer(12345);
        assert_eq!(err.to_string(), "Invalid k-mer containing non-ACGT base at position 12345");
        
        let err = REPriseError::invalid_fasta(42, "Missing sequence header");
        assert_eq!(err.to_string(), "Invalid FASTA format at line 42: Missing sequence header");
        
        let err = REPriseError::invalid_kmer_length(50, 4, 32);
        assert_eq!(err.to_string(), "K-mer length 50 is invalid (must be between 4 and 32)");
    }
    
    #[test]
    fn test_error_from_io() {
        let io_err = io::Error::new(io::ErrorKind::NotFound, "File not found");
        let reprise_err: REPriseError = io_err.into();
        
        match reprise_err {
            REPriseError::Io(_) => (),
            _ => panic!("Expected Io error"),
        }
    }
    
    #[test]
    fn test_helper_methods() {
        let err = REPriseError::concurrency("Channel closed");
        assert!(err.to_string().contains("Channel closed"));
        
        let err = REPriseError::config("Invalid thread count");
        assert!(err.to_string().contains("Invalid thread count"));
        
        let err = REPriseError::parse("Invalid number format");
        assert!(err.to_string().contains("Invalid number format"));
    }
}