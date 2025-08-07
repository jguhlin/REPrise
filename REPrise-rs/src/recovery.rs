//! Error recovery and resilience system for REPrise
//!
//! Provides checkpoint/resume functionality, automatic retry mechanisms,
//! and graceful error recovery for production deployments.

use crate::config::RecoverySettings;
use crate::error::{RepriseError, Result};
use crate::pipeline::{DetectedRepeat, PipelineStats};
use bincode;
use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};
use std::collections::VecDeque;
use std::fs::{File, OpenOptions};
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, RwLock};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};
use tracing::{debug, error, info, instrument, warn};
use uuid::Uuid;

/// Checkpoint data for resuming interrupted runs
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Checkpoint {
    /// Unique checkpoint identifier
    pub id: Uuid,
    /// Timestamp when checkpoint was created
    pub timestamp: DateTime<Utc>,
    /// Run identifier this checkpoint belongs to
    pub run_id: Uuid,
    /// Current processing state
    pub state: ProcessingState,
    /// Detected repeats so far
    pub detected_repeats: Vec<DetectedRepeat>,
    /// Pipeline statistics at checkpoint time
    pub pipeline_stats: CheckpointStats,
    /// Configuration hash to verify compatibility
    pub config_hash: u64,
    /// Input file hash to verify consistency
    pub input_hash: Option<String>,
    /// Progress information
    pub progress: ProgressInfo,
}

/// Processing state information for checkpoints
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProcessingState {
    /// Index of the last processed k-mer
    pub last_kmer_index: usize,
    /// Number of candidates generated so far
    pub candidates_generated: u64,
    /// Number of candidates processed so far
    pub candidates_processed: u64,
    /// Current pipeline phase
    pub phase: ProcessingPhase,
    /// Additional state data
    pub custom_data: serde_json::Value,
}

/// Current phase of pipeline processing
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum ProcessingPhase {
    /// Building genome index
    IndexBuilding,
    /// Generating candidate pairs
    CandidateGeneration,
    /// Processing candidates with alignment
    CandidateProcessing,
    /// Writing output files
    OutputWriting,
    /// Processing completed
    Completed,
}

/// Pipeline statistics for checkpoints (simplified version)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CheckpointStats {
    pub candidates_generated: u64,
    pub candidates_processed: u64,
    pub candidates_skipped: u64,
    pub repeats_detected: u64,
    pub processing_time_ms: u64,
    pub memory_usage: u64,
    pub peak_memory: u64,
}

impl From<&PipelineStats> for CheckpointStats {
    fn from(stats: &PipelineStats) -> Self {
        Self {
            candidates_generated: stats.candidates_generated(),
            candidates_processed: stats.candidates_processed(),
            candidates_skipped: stats.candidates_skipped(),
            repeats_detected: stats.repeats_detected(),
            processing_time_ms: stats.total_processing_time.load(Ordering::Relaxed),
            memory_usage: stats.memory_usage(),
            peak_memory: stats.peak_memory_usage(),
        }
    }
}

/// Progress information for user feedback
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProgressInfo {
    pub current_step: String,
    pub steps_completed: u64,
    pub total_steps: Option<u64>,
    pub percentage_complete: Option<f64>,
    pub estimated_time_remaining: Option<Duration>,
}

/// Retry configuration and state
#[derive(Debug, Clone)]
pub struct RetryState {
    pub attempts: u32,
    pub last_error: Option<String>,
    pub backoff_duration: Duration,
    pub next_retry: Option<Instant>,
}

/// Transient error types that can be retried
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TransientError {
    /// Temporary I/O error
    IoError,
    /// Memory allocation failure
    MemoryError,
    /// Network-related error
    NetworkError,
    /// Resource temporarily unavailable
    ResourceUnavailable,
    /// System overload
    SystemOverload,
}

/// Error recovery manager
pub struct RecoveryManager {
    settings: RecoverySettings,
    run_id: Uuid,
    checkpoint_counter: AtomicU64,
    last_checkpoint: RwLock<Option<Checkpoint>>,
    retry_states: RwLock<std::collections::HashMap<String, RetryState>>,
    shutdown_requested: Arc<AtomicBool>,
}

impl RecoveryManager {
    /// Create a new recovery manager
    pub fn new(settings: RecoverySettings) -> Self {
        Self {
            settings,
            run_id: Uuid::new_v4(),
            checkpoint_counter: AtomicU64::new(0),
            last_checkpoint: RwLock::new(None),
            retry_states: RwLock::new(std::collections::HashMap::new()),
            shutdown_requested: Arc::new(AtomicBool::new(false)),
        }
    }

    /// Get the run ID for this session
    pub fn run_id(&self) -> Uuid {
        self.run_id
    }

    /// Check if graceful shutdown was requested
    pub fn shutdown_requested(&self) -> bool {
        self.shutdown_requested.load(Ordering::Relaxed)
    }

    /// Request graceful shutdown
    pub fn request_shutdown(&self) {
        self.shutdown_requested.store(true, Ordering::Relaxed);
        info!("Graceful shutdown requested");
    }

    /// Create a checkpoint of current processing state
    #[instrument(skip(self, detected_repeats, pipeline_stats))]
    pub fn create_checkpoint(
        &self,
        state: ProcessingState,
        detected_repeats: &[DetectedRepeat],
        pipeline_stats: &PipelineStats,
        config_hash: u64,
        input_hash: Option<String>,
        progress: ProgressInfo,
    ) -> Result<()> {
        if !self.settings.enable_checkpoints {
            return Ok(());
        }

        let checkpoint_id = self.checkpoint_counter.fetch_add(1, Ordering::Relaxed);
        
        let checkpoint = Checkpoint {
            id: Uuid::new_v4(),
            timestamp: Utc::now(),
            run_id: self.run_id,
            state,
            detected_repeats: detected_repeats.to_vec(),
            pipeline_stats: CheckpointStats::from(pipeline_stats),
            config_hash,
            input_hash,
            progress,
        };

        // Write checkpoint to file
        let checkpoint_path = self.get_checkpoint_path(checkpoint_id)?;
        self.write_checkpoint(&checkpoint, &checkpoint_path)?;

        // Update last checkpoint
        *self.last_checkpoint.write().unwrap() = Some(checkpoint.clone());

        // Clean up old checkpoints (keep only the last few)
        self.cleanup_old_checkpoints(checkpoint_id)?;

        info!(
            checkpoint_id = %checkpoint.id,
            checkpoint_path = %checkpoint_path.display(),
            phase = ?checkpoint.state.phase,
            repeats_found = checkpoint.detected_repeats.len(),
            "Checkpoint created"
        );

        Ok(())
    }

    /// Load the most recent checkpoint for resuming
    pub fn load_latest_checkpoint(&self) -> Result<Option<Checkpoint>> {
        if !self.settings.enable_checkpoints {
            return Ok(None);
        }

        let checkpoint_dir = self.get_checkpoint_dir()?;
        if !checkpoint_dir.exists() {
            return Ok(None);
        }

        // Find the most recent checkpoint file
        let mut checkpoint_files: Vec<_> = std::fs::read_dir(&checkpoint_dir)?
            .filter_map(|entry| entry.ok())
            .filter(|entry| {
                entry.path().extension()
                    .and_then(|ext| ext.to_str()) == Some("checkpoint")
            })
            .collect();

        if checkpoint_files.is_empty() {
            return Ok(None);
        }

        // Sort by modification time (most recent first)
        checkpoint_files.sort_by_key(|entry| {
            entry.metadata()
                .and_then(|m| m.modified())
                .unwrap_or(SystemTime::UNIX_EPOCH)
        });
        checkpoint_files.reverse();

        // Try to load the most recent valid checkpoint
        for entry in checkpoint_files {
            match self.load_checkpoint(&entry.path()) {
                Ok(checkpoint) => {
                    info!(
                        checkpoint_id = %checkpoint.id,
                        timestamp = %checkpoint.timestamp,
                        phase = ?checkpoint.state.phase,
                        repeats_found = checkpoint.detected_repeats.len(),
                        "Loaded checkpoint for resume"
                    );
                    return Ok(Some(checkpoint));
                }
                Err(e) => {
                    warn!(
                        checkpoint_file = %entry.path().display(),
                        error = %e,
                        "Failed to load checkpoint, trying next"
                    );
                    continue;
                }
            }
        }

        warn!("No valid checkpoints found for resume");
        Ok(None)
    }

    /// Execute an operation with automatic retry on transient errors
    #[instrument(skip(self, operation))]
    pub fn with_retry<T, F>(&self, operation_name: &str, operation: F) -> Result<T>
    where
        F: Fn() -> Result<T>,
    {
        let mut retry_state = self.get_or_create_retry_state(operation_name);
        
        loop {
            // Check if we should wait before retrying
            if let Some(next_retry) = retry_state.next_retry {
                let now = Instant::now();
                if now < next_retry {
                    let wait_time = next_retry - now;
                    debug!(
                        operation = operation_name,
                        wait_time = ?wait_time,
                        attempt = retry_state.attempts + 1,
                        "Waiting before retry"
                    );
                    std::thread::sleep(wait_time);
                }
            }

            // Attempt the operation
            match operation() {
                Ok(result) => {
                    if retry_state.attempts > 0 {
                        info!(
                            operation = operation_name,
                            attempts = retry_state.attempts + 1,
                            "Operation succeeded after retry"
                        );
                    }
                    
                    // Reset retry state on success
                    self.reset_retry_state(operation_name);
                    return Ok(result);
                }
                Err(e) => {
                    // Check if this is a transient error
                    let is_transient = self.is_transient_error(&e);
                    retry_state.attempts += 1;
                    retry_state.last_error = Some(e.to_string());

                    if !is_transient || retry_state.attempts >= self.settings.max_retry_attempts {
                        error!(
                            operation = operation_name,
                            attempts = retry_state.attempts,
                            error = %e,
                            is_transient = is_transient,
                            "Operation failed after maximum retries"
                        );
                        return Err(e);
                    }

                    // Calculate backoff delay
                    let backoff = self.calculate_backoff(retry_state.attempts);
                    retry_state.backoff_duration = backoff;
                    retry_state.next_retry = Some(Instant::now() + backoff);

                    warn!(
                        operation = operation_name,
                        attempt = retry_state.attempts,
                        backoff = ?backoff,
                        error = %e,
                        "Transient error, will retry"
                    );

                    // Update retry state
                    self.update_retry_state(operation_name, retry_state.clone());
                }
            }
        }
    }

    /// Check if we can resume from a checkpoint
    pub fn can_resume_from_checkpoint(
        &self,
        checkpoint: &Checkpoint,
        config_hash: u64,
        input_hash: Option<&str>,
    ) -> bool {
        // Check configuration compatibility
        if checkpoint.config_hash != config_hash {
            warn!(
                checkpoint_config_hash = checkpoint.config_hash,
                current_config_hash = config_hash,
                "Configuration has changed, cannot resume from checkpoint"
            );
            return false;
        }

        // Check input file consistency if hash is available
        if let (Some(checkpoint_hash), Some(current_hash)) = 
            (checkpoint.input_hash.as_ref(), input_hash) {
            if checkpoint_hash != current_hash {
                warn!(
                    checkpoint_input_hash = checkpoint_hash,
                    current_input_hash = current_hash,
                    "Input file has changed, cannot resume from checkpoint"
                );
                return false;
            }
        }

        // Check if checkpoint is not too old
        let age = Utc::now().signed_duration_since(checkpoint.timestamp);
        let max_age = chrono::Duration::hours(24); // Don't resume from checkpoints older than 24h
        if age > max_age {
            warn!(
                checkpoint_age = ?age,
                max_age = ?max_age,
                "Checkpoint is too old, cannot resume"
            );
            return false;
        }

        true
    }

    /// Get checkpoint directory
    fn get_checkpoint_dir(&self) -> Result<PathBuf> {
        let dir = if let Some(checkpoint_dir) = &self.settings.checkpoint_dir {
            checkpoint_dir.clone()
        } else {
            std::env::temp_dir().join("reprise_checkpoints")
        };

        std::fs::create_dir_all(&dir)
            .map_err(|e| RepriseError::io_error(format!("Failed to create checkpoint directory: {}", e)))?;

        Ok(dir)
    }

    /// Get path for a specific checkpoint
    fn get_checkpoint_path(&self, checkpoint_id: u64) -> Result<PathBuf> {
        let dir = self.get_checkpoint_dir()?;
        let filename = format!("checkpoint_{}_{:06}.checkpoint", self.run_id, checkpoint_id);
        Ok(dir.join(filename))
    }

    /// Write checkpoint to file
    fn write_checkpoint(&self, checkpoint: &Checkpoint, path: &Path) -> Result<()> {
        let file = OpenOptions::new()
            .create(true)
            .truncate(true)
            .write(true)
            .open(path)
            .map_err(|e| RepriseError::io_error(format!("Failed to create checkpoint file: {}", e)))?;

        let writer = BufWriter::new(file);
        bincode::serialize_into(writer, checkpoint)
            .map_err(|e| RepriseError::config(format!("Failed to serialize checkpoint: {}", e)))?;

        Ok(())
    }

    /// Load checkpoint from file
    fn load_checkpoint(&self, path: &Path) -> Result<Checkpoint> {
        let file = File::open(path)
            .map_err(|e| RepriseError::io_error(format!("Failed to open checkpoint file: {}", e)))?;

        let reader = BufReader::new(file);
        bincode::deserialize_from(reader)
            .map_err(|e| RepriseError::config(format!("Failed to deserialize checkpoint: {}", e)))
    }

    /// Clean up old checkpoint files
    fn cleanup_old_checkpoints(&self, current_checkpoint_id: u64) -> Result<()> {
        let checkpoint_dir = self.get_checkpoint_dir()?;
        let keep_count = 5; // Keep last 5 checkpoints

        let mut checkpoint_files: Vec<_> = std::fs::read_dir(&checkpoint_dir)?
            .filter_map(|entry| entry.ok())
            .filter(|entry| {
                entry.file_name()
                    .to_str()
                    .map(|name| name.starts_with(&format!("checkpoint_{}", self.run_id)))
                    .unwrap_or(false)
            })
            .collect();

        // Sort by checkpoint ID (embedded in filename)
        checkpoint_files.sort_by_key(|entry| {
            entry.file_name()
                .to_str()
                .and_then(|name| name.split('_').nth(2))
                .and_then(|id_str| id_str.strip_suffix(".checkpoint"))
                .and_then(|id_str| id_str.parse::<u64>().ok())
                .unwrap_or(0)
        });

        // Remove old checkpoints, keeping only the most recent ones
        if checkpoint_files.len() > keep_count {
            let to_remove = checkpoint_files.len() - keep_count;
            for entry in &checkpoint_files[..to_remove] {
                if let Err(e) = std::fs::remove_file(entry.path()) {
                    warn!(
                        file = %entry.path().display(),
                        error = %e,
                        "Failed to remove old checkpoint file"
                    );
                }
            }
        }

        Ok(())
    }

    /// Get or create retry state for an operation
    fn get_or_create_retry_state(&self, operation_name: &str) -> RetryState {
        let retry_states = self.retry_states.read().unwrap();
        retry_states.get(operation_name).cloned().unwrap_or_else(|| {
            RetryState {
                attempts: 0,
                last_error: None,
                backoff_duration: self.settings.retry_base_delay,
                next_retry: None,
            }
        })
    }

    /// Update retry state for an operation
    fn update_retry_state(&self, operation_name: &str, state: RetryState) {
        let mut retry_states = self.retry_states.write().unwrap();
        retry_states.insert(operation_name.to_string(), state);
    }

    /// Reset retry state for an operation
    fn reset_retry_state(&self, operation_name: &str) {
        let mut retry_states = self.retry_states.write().unwrap();
        retry_states.remove(operation_name);
    }

    /// Check if an error is transient and should be retried
    fn is_transient_error(&self, error: &RepriseError) -> bool {
        match error {
            RepriseError::Io(_) => true,
            RepriseError::OutOfMemory(_) => true,
            RepriseError::Concurrency(_) => true,
            // Add more transient error patterns as needed
            _ => false,
        }
    }

    /// Calculate exponential backoff delay
    fn calculate_backoff(&self, attempt: u32) -> Duration {
        let base_delay = self.settings.retry_base_delay;
        let max_delay = self.settings.retry_max_delay;
        
        // Exponential backoff with jitter
        let exponential_delay = base_delay * (2_u32.pow(attempt.saturating_sub(1)));
        let jittered_delay = exponential_delay + Duration::from_millis(
            (rand::random::<u64>() % 1000) // Add up to 1 second of jitter
        );
        
        jittered_delay.min(max_delay)
    }
}

/// Signal handler for graceful shutdown
pub struct SignalHandler {
    recovery_manager: Arc<RecoveryManager>,
}

impl SignalHandler {
    /// Create a new signal handler
    pub fn new(recovery_manager: Arc<RecoveryManager>) -> Self {
        Self { recovery_manager }
    }

    /// Install signal handlers for graceful shutdown
    pub fn install_handlers(&self) -> Result<()> {
        if !self.recovery_manager.settings.enable_signal_handling {
            return Ok(());
        }

        use signal_hook::{consts::SIGTERM, consts::SIGINT, iterator::Signals};
        
        let recovery_manager = Arc::clone(&self.recovery_manager);
        let mut signals = Signals::new(&[SIGINT, SIGTERM])
            .map_err(|e| RepriseError::config(format!("Failed to register signal handlers: {}", e)))?;

        // Spawn signal handling thread
        std::thread::spawn(move || {
            for signal in signals.forever() {
                match signal {
                    SIGINT | SIGTERM => {
                        info!(signal = signal, "Received shutdown signal");
                        recovery_manager.request_shutdown();
                        break;
                    }
                    _ => {}
                }
            }
        });

        info!("Signal handlers installed for graceful shutdown");
        Ok(())
    }
}

// Add rand dependency placeholder (would normally be in Cargo.toml)
mod rand {
    pub fn random<T: Default>() -> T {
        // Simplified random number generation for testing
        // In production, use a proper random number generator
        T::default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pipeline::{GenomicRegion, PipelineStats};
    use crate::kmer::Kmer;
    use std::time::Duration;
    use tempfile::TempDir;

    fn create_test_settings() -> RecoverySettings {
        RecoverySettings {
            enable_checkpoints: true,
            checkpoint_interval: Duration::from_secs(60),
            checkpoint_dir: None,
            max_retry_attempts: 3,
            retry_base_delay: Duration::from_millis(100),
            retry_max_delay: Duration::from_secs(5),
            enable_signal_handling: false,
            shutdown_timeout: Duration::from_secs(10),
        }
    }

    fn create_test_checkpoint() -> Checkpoint {
        let region1 = GenomicRegion::new(1000, 2000, 0, 500, 1500);
        let region2 = GenomicRegion::new(3000, 4000, 1, 0, 1000);
        let kmer = Kmer::from_sequence(&[0, 1, 2, 3], 4).unwrap();
        
        let repeat = DetectedRepeat {
            region1,
            region2,
            score: 50,
            identity: 0.75,
            length: 1000,
            seed_kmer: kmer,
        };

        Checkpoint {
            id: Uuid::new_v4(),
            timestamp: Utc::now(),
            run_id: Uuid::new_v4(),
            state: ProcessingState {
                last_kmer_index: 100,
                candidates_generated: 1000,
                candidates_processed: 500,
                phase: ProcessingPhase::CandidateProcessing,
                custom_data: serde_json::json!({}),
            },
            detected_repeats: vec![repeat],
            pipeline_stats: CheckpointStats {
                candidates_generated: 1000,
                candidates_processed: 500,
                candidates_skipped: 50,
                repeats_detected: 1,
                processing_time_ms: 60000,
                memory_usage: 100000000,
                peak_memory: 150000000,
            },
            config_hash: 12345,
            input_hash: Some("abcdef123456".to_string()),
            progress: ProgressInfo {
                current_step: "Processing candidates".to_string(),
                steps_completed: 500,
                total_steps: Some(1000),
                percentage_complete: Some(50.0),
                estimated_time_remaining: Some(Duration::from_secs(60)),
            },
        }
    }

    #[test]
    fn test_recovery_manager_creation() {
        let settings = create_test_settings();
        let manager = RecoveryManager::new(settings);
        
        assert!(!manager.shutdown_requested());
        assert!(manager.run_id() != Uuid::nil());
    }

    #[test]
    fn test_checkpoint_creation_and_loading() {
        let temp_dir = TempDir::new().unwrap();
        let mut settings = create_test_settings();
        settings.checkpoint_dir = Some(temp_dir.path().to_path_buf());
        
        let manager = RecoveryManager::new(settings);
        let pipeline_stats = PipelineStats::new();
        
        let state = ProcessingState {
            last_kmer_index: 50,
            candidates_generated: 100,
            candidates_processed: 50,
            phase: ProcessingPhase::CandidateGeneration,
            custom_data: serde_json::json!({"test": "data"}),
        };
        
        let progress = ProgressInfo {
            current_step: "Testing".to_string(),
            steps_completed: 50,
            total_steps: Some(100),
            percentage_complete: Some(50.0),
            estimated_time_remaining: Some(Duration::from_secs(30)),
        };
        
        // Create checkpoint
        manager.create_checkpoint(
            state,
            &[],
            &pipeline_stats,
            12345,
            Some("test_hash".to_string()),
            progress,
        ).unwrap();
        
        // Load checkpoint
        let loaded = manager.load_latest_checkpoint().unwrap();
        assert!(loaded.is_some());
        
        let checkpoint = loaded.unwrap();
        assert_eq!(checkpoint.state.last_kmer_index, 50);
        assert_eq!(checkpoint.state.phase, ProcessingPhase::CandidateGeneration);
        assert_eq!(checkpoint.config_hash, 12345);
        assert_eq!(checkpoint.input_hash, Some("test_hash".to_string()));
    }

    #[test]
    fn test_checkpoint_compatibility() {
        let settings = create_test_settings();
        let manager = RecoveryManager::new(settings);
        let checkpoint = create_test_checkpoint();
        
        // Compatible checkpoint
        assert!(manager.can_resume_from_checkpoint(&checkpoint, 12345, Some("abcdef123456")));
        
        // Incompatible config hash
        assert!(!manager.can_resume_from_checkpoint(&checkpoint, 54321, Some("abcdef123456")));
        
        // Incompatible input hash
        assert!(!manager.can_resume_from_checkpoint(&checkpoint, 12345, Some("different_hash")));
    }

    #[test]
    fn test_retry_mechanism() {
        let settings = create_test_settings();
        let manager = RecoveryManager::new(settings);
        
        let mut attempt_count = 0;
        let result = manager.with_retry("test_operation", || {
            attempt_count += 1;
            if attempt_count < 3 {
                Err(RepriseError::io_error("Transient error".to_string()))
            } else {
                Ok(42)
            }
        });
        
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), 42);
        assert_eq!(attempt_count, 3);
    }

    #[test]
    fn test_retry_exhaustion() {
        let mut settings = create_test_settings();
        settings.max_retry_attempts = 2;
        let manager = RecoveryManager::new(settings);
        
        let mut attempt_count = 0;
        let result = manager.with_retry("test_operation", || {
            attempt_count += 1;
            Err(RepriseError::io_error("Persistent error".to_string()))
        });
        
        assert!(result.is_err());
        assert_eq!(attempt_count, 2); // Should try max_retry_attempts times
    }

    #[test]
    fn test_processing_phase_serialization() {
        let phase = ProcessingPhase::CandidateProcessing;
        let serialized = serde_json::to_string(&phase).unwrap();
        let deserialized: ProcessingPhase = serde_json::from_str(&serialized).unwrap();
        assert_eq!(phase, deserialized);
    }

    #[test]
    fn test_checkpoint_serialization() {
        let checkpoint = create_test_checkpoint();
        let serialized = bincode::serialize(&checkpoint).unwrap();
        let deserialized: Checkpoint = bincode::deserialize(&serialized).unwrap();
        
        assert_eq!(checkpoint.id, deserialized.id);
        assert_eq!(checkpoint.state.phase, deserialized.state.phase);
        assert_eq!(checkpoint.detected_repeats.len(), deserialized.detected_repeats.len());
    }

    #[test]
    fn test_backoff_calculation() {
        let settings = create_test_settings();
        let manager = RecoveryManager::new(settings);
        
        let backoff1 = manager.calculate_backoff(1);
        let backoff2 = manager.calculate_backoff(2);
        let backoff3 = manager.calculate_backoff(3);
        
        // Exponential backoff should increase
        assert!(backoff1 <= backoff2);
        assert!(backoff2 <= backoff3);
        
        // Should not exceed max delay
        let max_backoff = manager.calculate_backoff(10);
        assert!(max_backoff <= manager.settings.retry_max_delay);
    }

    #[test]
    fn test_shutdown_request() {
        let settings = create_test_settings();
        let manager = RecoveryManager::new(settings);
        
        assert!(!manager.shutdown_requested());
        manager.request_shutdown();
        assert!(manager.shutdown_requested());
    }
}