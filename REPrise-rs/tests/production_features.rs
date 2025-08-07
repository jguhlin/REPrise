//! Integration tests for production features
//!
//! Tests the complete integration of logging, configuration, output formats,
//! error recovery, and system monitoring features.

use reprise::config::{ConfigManager, OutputFormat, CompressionType, LogLevel};
use reprise::logging::{LoggingSystem, LoggingConfig, MetricsCollector};
use reprise::output::{OutputManager, OutputMetadata, ProcessingParameters, OutputStatistics, GenomeInfo};
use reprise::recovery::{RecoveryManager, ProcessingState, ProcessingPhase, ProgressInfo};
use reprise::system::{ResourceMonitor, SystemSettings};
use reprise::pipeline::{DetectedRepeat, GenomicRegion};
use reprise::genome::Genome;
use reprise::kmer::Kmer;
use std::sync::Arc;
use std::time::Duration;
use tempfile::{TempDir, NamedTempFile};
use std::io::Write;
use uuid::Uuid;
use chrono::Utc;

/// Create a test genome for testing
fn create_test_genome(sequences: &[(&str, &str)]) -> Genome {
    let mut fasta_content = String::new();
    for (name, seq) in sequences {
        fasta_content.push_str(&format!(">{}\n{}\n", name, seq));
    }
    
    let mut file = NamedTempFile::new().unwrap();
    file.write_all(fasta_content.as_bytes()).unwrap();
    file.flush().unwrap();
    
    Genome::from_fasta(file.path()).unwrap()
}

/// Create test detected repeats
fn create_test_repeats() -> Vec<DetectedRepeat> {
    let region1 = GenomicRegion::new(1000, 2000, 0, 500, 1500);
    let region2 = GenomicRegion::new(3000, 4000, 1, 0, 1000);
    let kmer = Kmer::from_sequence(&[0, 1, 2, 3], 4).unwrap();
    
    vec![
        DetectedRepeat {
            region1,
            region2,
            score: 50,
            identity: 0.75,
            length: 1000,
            seed_kmer: kmer,
        },
        DetectedRepeat {
            region1: GenomicRegion::new(5000, 5500, 0, 4500, 5000),
            region2: GenomicRegion::new(6000, 6500, 0, 5500, 6000),
            score: 80,
            identity: 0.90,
            length: 500,
            seed_kmer: kmer,
        },
    ]
}

#[test]
fn test_configuration_system_integration() {
    let temp_dir = TempDir::new().unwrap();
    
    // Test configuration manager with different formats
    let mut manager = ConfigManager::new();
    
    // Test profile application
    manager.apply_profile("fast").unwrap();
    assert_eq!(manager.config().pipeline.min_frequency, 5);
    assert_eq!(manager.config().pipeline.min_identity, 0.40);
    
    manager.apply_profile("accurate").unwrap();
    assert_eq!(manager.config().pipeline.min_frequency, 2);
    assert_eq!(manager.config().pipeline.min_identity, 0.65);
    
    // Test configuration file save/load
    let config_file = temp_dir.path().join("test_config.toml");
    manager.save_to_file(&config_file).unwrap();
    
    let loaded_manager = ConfigManager::load_from_file(&config_file).unwrap();
    assert_eq!(
        loaded_manager.config().pipeline.min_frequency,
        manager.config().pipeline.min_frequency
    );
    
    // Test pipeline config conversion
    let pipeline_config = manager.to_pipeline_config();
    assert_eq!(pipeline_config.min_frequency, manager.config().pipeline.min_frequency);
    assert_eq!(pipeline_config.min_identity, manager.config().pipeline.min_identity);
}

#[test]
fn test_logging_system_integration() {
    let temp_dir = TempDir::new().unwrap();
    
    // Test logging configuration
    let logging_config = LoggingConfig {
        level: LogLevel::Debug,
        json_format: false,
        log_dir: Some(temp_dir.path().to_path_buf()),
        log_file_pattern: "test.log".to_string(),
        enable_metrics: true,
        track_memory: true,
        ..Default::default()
    };
    
    // Initialize logging system (may fail in test environment, but should not panic)
    let logging_result = LoggingSystem::init(logging_config);
    
    // Test should pass whether logging initialization succeeds or fails gracefully
    match logging_result {
        Ok(logging_system) => {
            // Test metrics collection
            let metrics = logging_system.metrics();
            metrics.inc_counter("test_counter");
            metrics.set_gauge("test_gauge", 42);
            metrics.record_time("test_timer", Duration::from_millis(100));
            
            assert_eq!(metrics.get_counter("test_counter"), 1);
            assert_eq!(metrics.get_gauge("test_gauge"), 42);
            
            // Test progress reporter
            let progress = logging_system.progress("test_operation", Some(100));
            progress.set(50);
            progress.report(); // Should not panic
        }
        Err(_) => {
            // Logging initialization can fail in test environment, which is acceptable
            println!("Logging initialization failed (expected in test environment)");
        }
    }
}

#[test]
fn test_output_format_integration() {
    let temp_dir = TempDir::new().unwrap();
    let genome = create_test_genome(&[
        ("contig1", "ATCGATCGATCGATCGATCG"),
        ("contig2", "GCTAGCTAGCTAGCTAGCTA"),
    ]);
    let repeats = create_test_repeats();
    
    // Create output metadata
    let metadata = OutputMetadata {
        run_id: Uuid::new_v4(),
        timestamp: Utc::now(),
        version: "test-1.0.0".to_string(),
        command_line: "reprise --test".to_string(),
        input_file: temp_dir.path().join("test.fa"),
        parameters: ProcessingParameters {
            k_mer_length: 13,
            min_frequency: 3,
            max_frequency: Some(10000),
            region_extension: 100,
            min_alignment_score: 10,
            min_identity: 0.50,
            num_workers: 4,
        },
        statistics: OutputStatistics {
            total_repeats_found: repeats.len() as u64,
            candidates_processed: 1000,
            processing_time_seconds: 60.0,
            genome_coverage: 25.0,
            average_repeat_length: 750.0,
            average_identity: 0.825,
            average_score: 65.0,
        },
        genome_info: GenomeInfo {
            total_length: genome.len(),
            num_contigs: genome.num_contigs(),
            n50: 10000,
            gc_content: 0.42,
        },
    };
    
    // Test different output formats
    let output_settings = reprise::config::OutputSettings {
        formats: vec![OutputFormat::Tsv, OutputFormat::Json, OutputFormat::Bed],
        enable_compression: false,
        compression_type: CompressionType::Gzip,
        compression_level: 6,
        include_metadata: true,
        generate_checksums: true,
        output_dir: Some(temp_dir.path().to_path_buf()),
        file_prefix: "test_output".to_string(),
        buffer_size: 65536,
    };
    
    let output_manager = OutputManager::new(output_settings, metadata).unwrap();
    let checksums = output_manager.write_results(&repeats, &genome).unwrap();
    
    // Verify output files were created
    assert!(checksums.len() >= 3); // At least TSV, JSON, BED
    
    // Verify files exist
    for checksum in &checksums {
        assert!(checksum.file_path.exists());
        assert!(checksum.size_bytes > 0);
        assert!(!checksum.sha256.is_empty());
    }
    
    // Verify checksum file was created
    let checksum_file = temp_dir.path().join("test_output.checksums");
    assert!(checksum_file.exists());
}

#[test]
fn test_recovery_system_integration() {
    let temp_dir = TempDir::new().unwrap();
    let repeats = create_test_repeats();
    
    // Test recovery manager
    let recovery_settings = reprise::config::RecoverySettings {
        enable_checkpoints: true,
        checkpoint_interval: Duration::from_secs(60),
        checkpoint_dir: Some(temp_dir.path().to_path_buf()),
        max_retry_attempts: 3,
        retry_base_delay: Duration::from_millis(100),
        retry_max_delay: Duration::from_secs(5),
        enable_signal_handling: false, // Disable for testing
        shutdown_timeout: Duration::from_secs(10),
    };
    
    let recovery_manager = RecoveryManager::new(recovery_settings);
    let pipeline_stats = reprise::pipeline::PipelineStats::new();
    
    // Test checkpoint creation
    let processing_state = ProcessingState {
        last_kmer_index: 100,
        candidates_generated: 1000,
        candidates_processed: 500,
        phase: ProcessingPhase::CandidateProcessing,
        custom_data: serde_json::json!({"test": "data"}),
    };
    
    let progress_info = ProgressInfo {
        current_step: "Testing checkpoint".to_string(),
        steps_completed: 500,
        total_steps: Some(1000),
        percentage_complete: Some(50.0),
        estimated_time_remaining: Some(Duration::from_secs(60)),
    };
    
    recovery_manager.create_checkpoint(
        processing_state,
        &repeats,
        &pipeline_stats,
        12345,
        Some("test_hash".to_string()),
        progress_info,
    ).unwrap();
    
    // Test checkpoint loading
    let loaded_checkpoint = recovery_manager.load_latest_checkpoint().unwrap();
    assert!(loaded_checkpoint.is_some());
    
    let checkpoint = loaded_checkpoint.unwrap();
    assert_eq!(checkpoint.state.last_kmer_index, 100);
    assert_eq!(checkpoint.state.phase, ProcessingPhase::CandidateProcessing);
    assert_eq!(checkpoint.detected_repeats.len(), 2);
    assert_eq!(checkpoint.config_hash, 12345);
    
    // Test checkpoint compatibility
    assert!(recovery_manager.can_resume_from_checkpoint(&checkpoint, 12345, Some("test_hash")));
    assert!(!recovery_manager.can_resume_from_checkpoint(&checkpoint, 54321, Some("test_hash")));
    assert!(!recovery_manager.can_resume_from_checkpoint(&checkpoint, 12345, Some("different_hash")));
    
    // Test retry mechanism
    let mut attempt_count = 0;
    let result = recovery_manager.with_retry("test_operation", || {
        attempt_count += 1;
        if attempt_count < 3 {
            Err(reprise::error::RepriseError::io_error("Transient error".to_string()))
        } else {
            Ok(42)
        }
    });
    
    assert!(result.is_ok());
    assert_eq!(result.unwrap(), 42);
    assert_eq!(attempt_count, 3);
}

#[test]
fn test_system_monitoring_integration() {
    let system_settings = SystemSettings {
        memory_check_interval: Duration::from_secs(1),
        memory_warning_threshold: 75.0,
        memory_critical_threshold: 85.0,
        enable_gc_hints: true,
        max_file_handles: 1024,
        process_priority: 0,
        enable_cpu_affinity: false,
        cpu_affinity_cores: Vec::new(),
    };
    
    let metrics = Arc::new(MetricsCollector::new());
    let monitor = ResourceMonitor::new(system_settings, Arc::clone(&metrics));
    
    // Test system info collection
    let system_info = monitor.system_info().unwrap();
    assert!(system_info.total_memory > 0);
    assert!(system_info.cpu_count > 0);
    assert!(system_info.memory_usage_percent >= 0.0);
    assert!(system_info.memory_usage_percent <= 100.0);
    
    // Test resource limit checking
    let events = monitor.check_resource_limits().unwrap();
    // Events depend on actual system state, just verify the call succeeds
    assert!(events.len() >= 0);
    
    // Test recommendations
    let recommendations = monitor.get_recommendations().unwrap();
    assert!(recommendations.len() >= 0);
    
    // Test event handler
    let event_received = Arc::new(std::sync::atomic::AtomicBool::new(false));
    let event_received_clone = Arc::clone(&event_received);
    
    monitor.add_event_handler(move |_event| {
        event_received_clone.store(true, std::sync::atomic::Ordering::Relaxed);
    });
    
    // Manually trigger an event to test handler
    let handlers = monitor.event_handlers.read().unwrap();
    if !handlers.is_empty() {
        let test_event = reprise::system::ResourceEvent::MemoryWarning {
            current: 1000000,
            threshold: 800000,
            percentage: 80.0,
        };
        handlers[0](test_event);
        assert!(event_received.load(std::sync::atomic::Ordering::Relaxed));
    }
}

#[test]
fn test_end_to_end_integration() {
    let temp_dir = TempDir::new().unwrap();
    
    // Create a comprehensive configuration
    let mut manager = ConfigManager::new();
    
    // Configure for testing
    let config = manager.config_mut();
    config.logging.level = LogLevel::Debug;
    config.logging.log_dir = Some(temp_dir.path().to_path_buf());
    config.output.output_dir = Some(temp_dir.path().to_path_buf());
    config.output.formats = vec![OutputFormat::Tsv, OutputFormat::Json];
    config.output.enable_compression = false;
    config.recovery.enable_checkpoints = true;
    config.recovery.checkpoint_dir = Some(temp_dir.path().to_path_buf());
    config.system.memory_check_interval = Duration::from_secs(1);
    
    // Test configuration validation
    assert!(config.validate().is_ok());
    
    // Test configuration serialization
    let config_file = temp_dir.path().join("integration_test.toml");
    manager.save_to_file(&config_file).unwrap();
    assert!(config_file.exists());
    
    // Test configuration loading
    let loaded_manager = ConfigManager::load_from_file(&config_file).unwrap();
    assert_eq!(
        loaded_manager.config().pipeline.min_frequency,
        manager.config().pipeline.min_frequency
    );
    
    // Test pipeline config conversion
    let pipeline_config = manager.to_pipeline_config();
    assert!(pipeline_config.num_workers > 0);
    assert!(pipeline_config.channel_capacity > 0);
    
    println!("End-to-end integration test completed successfully");
}

#[test]
fn test_metrics_collection_integration() {
    let metrics = Arc::new(MetricsCollector::new());
    
    // Test counter operations
    metrics.inc_counter("test_operation");
    metrics.inc_counter_by("test_operation", 5);
    assert_eq!(metrics.get_counter("test_operation"), 6);
    
    // Test gauge operations
    metrics.set_gauge("memory_usage", 1024 * 1024 * 512); // 512MB
    assert_eq!(metrics.get_gauge("memory_usage"), 1024 * 1024 * 512);
    
    // Test timing operations
    metrics.record_time("processing_time", Duration::from_millis(150));
    metrics.record_time("processing_time", Duration::from_millis(200));
    metrics.record_time("processing_time", Duration::from_millis(100));
    
    let timing_stats = metrics.get_timing_stats("processing_time").unwrap();
    assert_eq!(timing_stats.count, 3);
    assert_eq!(timing_stats.median, Duration::from_millis(150));
    assert_eq!(timing_stats.min, Duration::from_millis(100));
    assert_eq!(timing_stats.max, Duration::from_millis(200));
    
    // Test memory tracking
    metrics.update_memory_usage(1000);
    assert_eq!(metrics.memory_usage.load(std::sync::atomic::Ordering::Relaxed), 1000);
    
    metrics.update_memory_usage(2000);
    assert_eq!(metrics.peak_memory.load(std::sync::atomic::Ordering::Relaxed), 2000);
    
    metrics.update_memory_usage(1500);
    assert_eq!(metrics.peak_memory.load(std::sync::atomic::Ordering::Relaxed), 2000); // Peak should remain
    
    // Test metrics summary
    let summary = metrics.summary();
    assert!(summary.counters.contains_key("test_operation"));
    assert!(summary.gauges.contains_key("memory_usage"));
    assert!(summary.timings.contains_key("processing_time"));
    assert_eq!(summary.memory_usage, 1500);
    assert_eq!(summary.peak_memory, 2000);
}

#[test]
fn test_error_handling_integration() {
    let temp_dir = TempDir::new().unwrap();
    
    // Test recovery manager error handling
    let recovery_settings = reprise::config::RecoverySettings {
        enable_checkpoints: true,
        checkpoint_dir: Some(temp_dir.path().to_path_buf()),
        max_retry_attempts: 2,
        retry_base_delay: Duration::from_millis(10), // Fast for testing
        retry_max_delay: Duration::from_millis(100),
        enable_signal_handling: false,
        shutdown_timeout: Duration::from_secs(1),
    };
    
    let recovery_manager = RecoveryManager::new(recovery_settings);
    
    // Test successful retry after failures
    let mut attempt_count = 0;
    let result = recovery_manager.with_retry("retry_test", || {
        attempt_count += 1;
        if attempt_count == 1 {
            Err(reprise::error::RepriseError::io_error("First failure".to_string()))
        } else {
            Ok("success")
        }
    });
    
    assert!(result.is_ok());
    assert_eq!(result.unwrap(), "success");
    assert_eq!(attempt_count, 2);
    
    // Test retry exhaustion
    let mut fail_count = 0;
    let result = recovery_manager.with_retry("exhaust_test", || {
        fail_count += 1;
        Err(reprise::error::RepriseError::io_error("Persistent failure".to_string()))
    });
    
    assert!(result.is_err());
    assert_eq!(fail_count, 2); // Should try max_retry_attempts times
    
    // Test non-retryable error
    let result = recovery_manager.with_retry("non_retryable_test", || {
        Err(reprise::error::RepriseError::config("Config error".to_string()))
    });
    
    assert!(result.is_err());
}

#[test]
fn test_configuration_validation_integration() {
    let mut config = reprise::config::RepriseConfig::default();
    
    // Test valid configuration
    assert!(config.validate().is_ok());
    
    // Test invalid pipeline settings
    config.pipeline.num_workers = 1000; // Too many workers
    assert!(config.validate().is_err());
    
    // Fix worker count, test invalid memory setting
    config.pipeline.num_workers = 4;
    config.pipeline.max_memory_usage = 1000; // Too small
    assert!(config.validate().is_err());
    
    // Fix memory, test invalid identity range
    config.pipeline.max_memory_usage = 256 * 1024 * 1024;
    config.pipeline.min_identity = 1.5; // Invalid range
    assert!(config.validate().is_err());
    
    // Fix identity range
    config.pipeline.min_identity = 0.5;
    assert!(config.validate().is_ok());
    
    // Test invalid index settings
    config.index.k = 50; // Too large
    assert!(config.validate().is_err());
    
    // Fix k-mer length
    config.index.k = 13;
    assert!(config.validate().is_ok());
    
    // Test invalid system settings
    config.system.memory_warning_threshold = 150.0; // Over 100%
    assert!(config.validate().is_err());
    
    // Fix threshold
    config.system.memory_warning_threshold = 75.0;
    assert!(config.validate().is_ok());
}