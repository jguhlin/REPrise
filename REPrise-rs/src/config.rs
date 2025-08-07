//! Configuration management system for REPrise
//!
//! Provides comprehensive configuration loading, validation, and management
//! with support for YAML/TOML files, environment variables, and preset profiles.

use crate::error::{RepriseError, Result};
use crate::logging::LoggingConfig;
use crate::pipeline::PipelineConfig;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::env;
use std::path::{Path, PathBuf};
use std::time::Duration;
use validator::{Validate, ValidationError};

/// Main configuration structure for REPrise
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct RepriseConfig {
    /// Logging configuration
    #[serde(default)]
    pub logging: LoggingConfig,

    /// Pipeline processing configuration
    #[serde(default)]
    #[validate(nested)]
    pub pipeline: PipelineSettings,

    /// Index building configuration
    #[serde(default)]
    #[validate(nested)]
    pub index: IndexSettings,

    /// Output configuration
    #[serde(default)]
    #[validate(nested)]
    pub output: OutputSettings,

    /// System resource configuration
    #[serde(default)]
    #[validate(nested)]
    pub system: SystemSettings,

    /// Error recovery and resilience settings
    #[serde(default)]
    #[validate(nested)]
    pub recovery: RecoverySettings,
}

/// Pipeline processing settings
#[derive(Debug, Clone, Serialize, Deserialize, Validate)]
pub struct PipelineSettings {
    /// Number of worker threads (0 = auto-detect)
    #[validate(range(min = 0, max = 256))]
    pub num_workers: usize,

    /// Channel capacity for bounded channels
    #[validate(range(min = 100, max = 10000000))]
    pub channel_capacity: usize,

    /// Batch size for candidate processing
    #[validate(range(min = 1, max = 10000))]
    pub batch_size: usize,

    /// Maximum memory usage in bytes
    #[validate(range(min = 100000000, max = 1000000000000))] // 100MB - 1TB
    pub max_memory_usage: u64,

    /// Minimum k-mer frequency for processing
    #[validate(range(min = 1, max = 1000))]
    pub min_frequency: u32,

    /// Maximum k-mer frequency for processing
    pub max_frequency: Option<u32>,

    /// Region extension size around k-mer positions
    #[validate(range(min = 10, max = 100000))]
    pub region_extension: u64,

    /// Maximum region size to process
    #[validate(range(min = 100, max = 1000000))]
    pub max_region_size: u64,

    /// Enable backpressure control
    pub enable_backpressure: bool,

    /// Backpressure threshold (0.0-1.0)
    #[validate(range(min = 0.1, max = 1.0))]
    pub backpressure_threshold: f64,

    /// Minimum alignment score for repeat detection
    #[validate(range(min = 1, max = 10000))]
    pub min_alignment_score: i32,

    /// Minimum percent identity for repeat detection (0.0-1.0)
    #[validate(range(min = 0.1, max = 1.0))]
    pub min_identity: f64,
}

impl Default for PipelineSettings {
    fn default() -> Self {
        Self {
            num_workers: 0,
            channel_capacity: 10000,
            batch_size: 100,
            max_memory_usage: 256 * 1024 * 1024, // 256MB
            min_frequency: 3,
            max_frequency: Some(10000),
            region_extension: 100,
            max_region_size: 10000,
            enable_backpressure: true,
            backpressure_threshold: 0.8,
            min_alignment_score: 10,
            min_identity: 0.50,
        }
    }
}

/// Index building settings
#[derive(Debug, Clone, Serialize, Deserialize, Validate)]
pub struct IndexSettings {
    /// K-mer length (0 = auto-calculate)
    #[validate(range(min = 0, max = 32))]
    pub k: usize,

    /// Minimum k-mer frequency for indexing
    #[validate(range(min = 1, max = 1000))]
    pub min_frequency: u32,

    /// Maximum k-mer frequency for indexing
    pub max_frequency: Option<u32>,

    /// Enable parallel index building
    pub parallel: bool,

    /// Maximum positions per k-mer to store
    #[validate(range(min = 10, max = 100000))]
    pub max_positions_per_kmer: u32,

    /// Expected memory usage hint in bytes
    pub expected_memory_usage: Option<u64>,
}

impl Default for IndexSettings {
    fn default() -> Self {
        Self {
            k: 0, // Auto-calculate
            min_frequency: 3,
            max_frequency: Some(10000),
            parallel: true,
            max_positions_per_kmer: 10000,
            expected_memory_usage: None,
        }
    }
}

/// Output format configuration
#[derive(Debug, Clone, Serialize, Deserialize, Validate)]
pub struct OutputSettings {
    /// Available output formats
    pub formats: Vec<OutputFormat>,

    /// Enable output compression
    pub enable_compression: bool,

    /// Compression type
    pub compression_type: CompressionType,

    /// Compression level (1-9, higher = better compression, slower)
    #[validate(range(min = 1, max = 9))]
    pub compression_level: u32,

    /// Include metadata headers in output files
    pub include_metadata: bool,

    /// Generate checksums for output files
    pub generate_checksums: bool,

    /// Output directory
    pub output_dir: Option<PathBuf>,

    /// Output file prefix
    pub file_prefix: String,

    /// Buffer size for output writing
    #[validate(range(min = 4096, max = 10485760))] // 4KB - 10MB
    pub buffer_size: usize,
}

impl Default for OutputSettings {
    fn default() -> Self {
        Self {
            formats: vec![OutputFormat::Tsv, OutputFormat::Bed],
            enable_compression: false,
            compression_type: CompressionType::Gzip,
            compression_level: 6,
            include_metadata: true,
            generate_checksums: true,
            output_dir: None,
            file_prefix: "reprise_output".to_string(),
            buffer_size: 65536, // 64KB
        }
    }
}

/// Supported output formats
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum OutputFormat {
    /// Tab-separated values
    Tsv,
    /// Browser Extensible Data format
    Bed,
    /// JavaScript Object Notation
    Json,
    /// General Feature Format version 3
    Gff3,
    /// Custom REPrise format
    Reprise,
}

impl std::fmt::Display for OutputFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            OutputFormat::Tsv => write!(f, "tsv"),
            OutputFormat::Bed => write!(f, "bed"),
            OutputFormat::Json => write!(f, "json"),
            OutputFormat::Gff3 => write!(f, "gff3"),
            OutputFormat::Reprise => write!(f, "reprof"),
        }
    }
}

/// Supported compression types
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum CompressionType {
    /// No compression
    None,
    /// Gzip compression
    Gzip,
    /// Bzip2 compression
    Bzip2,
}

/// System resource settings
#[derive(Debug, Clone, Serialize, Deserialize, Validate)]
pub struct SystemSettings {
    /// Memory usage monitoring interval
    #[serde(with = "duration_serde")]
    pub memory_check_interval: Duration,

    /// Memory usage warning threshold (percentage)
    #[validate(range(min = 50.0, max = 100.0))]
    pub memory_warning_threshold: f64,

    /// Memory usage critical threshold (percentage)
    #[validate(range(min = 60.0, max = 100.0))]
    pub memory_critical_threshold: f64,

    /// Enable automatic garbage collection
    pub enable_gc_hints: bool,

    /// Maximum number of file handles to use
    #[validate(range(min = 10, max = 65536))]
    pub max_file_handles: usize,

    /// Process priority (nice value on Unix)
    #[validate(range(min = -20, max = 19))]
    pub process_priority: i8,

    /// Enable CPU affinity setting
    pub enable_cpu_affinity: bool,

    /// CPU cores to bind to (empty = no affinity)
    pub cpu_affinity_cores: Vec<usize>,
}

impl Default for SystemSettings {
    fn default() -> Self {
        Self {
            memory_check_interval: Duration::from_secs(10),
            memory_warning_threshold: 80.0,
            memory_critical_threshold: 90.0,
            enable_gc_hints: true,
            max_file_handles: 1024,
            process_priority: 0,
            enable_cpu_affinity: false,
            cpu_affinity_cores: Vec::new(),
        }
    }
}

/// Error recovery and resilience settings
#[derive(Debug, Clone, Serialize, Deserialize, Validate)]
pub struct RecoverySettings {
    /// Enable checkpoint/resume functionality
    pub enable_checkpoints: bool,

    /// Checkpoint save interval
    #[serde(with = "duration_serde")]
    pub checkpoint_interval: Duration,

    /// Checkpoint directory
    pub checkpoint_dir: Option<PathBuf>,

    /// Maximum number of retry attempts for transient errors
    #[validate(range(min = 0, max = 10))]
    pub max_retry_attempts: u32,

    /// Base delay for retry backoff
    #[serde(with = "duration_serde")]
    pub retry_base_delay: Duration,

    /// Maximum delay for retry backoff
    #[serde(with = "duration_serde")]
    pub retry_max_delay: Duration,

    /// Enable graceful shutdown on signals
    pub enable_signal_handling: bool,

    /// Shutdown timeout
    #[serde(with = "duration_serde")]
    pub shutdown_timeout: Duration,
}

impl Default for RecoverySettings {
    fn default() -> Self {
        Self {
            enable_checkpoints: false,
            checkpoint_interval: Duration::from_secs(300), // 5 minutes
            checkpoint_dir: None,
            max_retry_attempts: 3,
            retry_base_delay: Duration::from_millis(100),
            retry_max_delay: Duration::from_secs(30),
            enable_signal_handling: true,
            shutdown_timeout: Duration::from_secs(30),
        }
    }
}

impl Default for RepriseConfig {
    fn default() -> Self {
        Self {
            logging: LoggingConfig::default(),
            pipeline: PipelineSettings::default(),
            index: IndexSettings::default(),
            output: OutputSettings::default(),
            system: SystemSettings::default(),
            recovery: RecoverySettings::default(),
        }
    }
}

impl Validate for RepriseConfig {
    fn validate(&self) -> Result<(), validator::ValidationErrors> {
        use validator::ValidationErrors;
        
        let mut errors = ValidationErrors::new();
        
        if let Err(e) = self.pipeline.validate() {
            errors.merge(e);
        }
        
        if let Err(e) = self.index.validate() {
            errors.merge(e);
        }
        
        if let Err(e) = self.output.validate() {
            errors.merge(e);
        }
        
        if let Err(e) = self.system.validate() {
            errors.merge(e);
        }
        
        if let Err(e) = self.recovery.validate() {
            errors.merge(e);
        }
        
        if errors.is_empty() {
            Ok(())
        } else {
            Err(errors)
        }
    }
}

/// Configuration profiles for different use cases
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConfigProfile {
    pub name: String,
    pub description: String,
    pub config: RepriseConfig,
}

/// Configuration manager for REPrise
pub struct ConfigManager {
    config: RepriseConfig,
    profiles: HashMap<String, ConfigProfile>,
}

impl ConfigManager {
    /// Create a new configuration manager with default settings
    pub fn new() -> Self {
        let mut manager = Self {
            config: RepriseConfig::default(),
            profiles: HashMap::new(),
        };
        manager.load_builtin_profiles();
        manager
    }

    /// Load configuration from file
    pub fn load_from_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let content = std::fs::read_to_string(path)
            .map_err(|e| RepriseError::io_error(format!("Failed to read config file: {}", e)))?;

        let config: RepriseConfig = match path.extension().and_then(|ext| ext.to_str()) {
            Some("toml") => toml::from_str(&content)
                .map_err(|e| RepriseError::config(format!("TOML parse error: {}", e)))?,
            Some("yaml") | Some("yml") => serde_yaml::from_str(&content)
                .map_err(|e| RepriseError::config(format!("YAML parse error: {}", e)))?,
            Some("json") => serde_json::from_str(&content)
                .map_err(|e| RepriseError::config(format!("JSON parse error: {}", e)))?,
            _ => {
                return Err(RepriseError::config(
                    "Unsupported config file format. Use .toml, .yaml, .yml, or .json".to_string(),
                ))
            }
        };

        // Validate configuration
        config
            .validate()
            .map_err(|e| RepriseError::config(format!("Configuration validation failed: {}", e)))?;

        let mut manager = Self {
            config,
            profiles: HashMap::new(),
        };
        manager.load_builtin_profiles();
        Ok(manager)
    }

    /// Load configuration from environment variables and merge with current config
    pub fn load_from_env(&mut self) -> Result<()> {
        // Override configuration values from environment variables
        if let Ok(workers) = env::var("REPRISE_WORKERS") {
            self.config.pipeline.num_workers = workers.parse()
                .map_err(|e| RepriseError::config(format!("Invalid REPRISE_WORKERS: {}", e)))?;
        }

        if let Ok(memory) = env::var("REPRISE_MAX_MEMORY") {
            self.config.pipeline.max_memory_usage = memory.parse()
                .map_err(|e| RepriseError::config(format!("Invalid REPRISE_MAX_MEMORY: {}", e)))?;
        }

        if let Ok(log_level) = env::var("REPRISE_LOG_LEVEL") {
            self.config.logging.level = match log_level.to_lowercase().as_str() {
                "trace" => crate::logging::LogLevel::Trace,
                "debug" => crate::logging::LogLevel::Debug,
                "info" => crate::logging::LogLevel::Info,
                "warn" => crate::logging::LogLevel::Warn,
                "error" => crate::logging::LogLevel::Error,
                _ => return Err(RepriseError::config(format!("Invalid log level: {}", log_level))),
            };
        }

        if let Ok(json_logs) = env::var("REPRISE_JSON_LOGS") {
            self.config.logging.json_format = json_logs.parse()
                .map_err(|e| RepriseError::config(format!("Invalid REPRISE_JSON_LOGS: {}", e)))?;
        }

        if let Ok(output_dir) = env::var("REPRISE_OUTPUT_DIR") {
            self.config.output.output_dir = Some(PathBuf::from(output_dir));
        }

        // Validate after environment overrides
        self.config
            .validate()
            .map_err(|e| RepriseError::config(format!("Configuration validation failed after env override: {}", e)))?;

        Ok(())
    }

    /// Apply a configuration profile
    pub fn apply_profile(&mut self, profile_name: &str) -> Result<()> {
        let profile = self.profiles.get(profile_name).ok_or_else(|| {
            RepriseError::config(format!("Unknown profile: {}", profile_name))
        })?.clone();

        self.config = profile.config;
        self.config
            .validate()
            .map_err(|e| RepriseError::config(format!("Profile validation failed: {}", e)))?;

        Ok(())
    }

    /// Get current configuration
    pub fn config(&self) -> &RepriseConfig {
        &self.config
    }

    /// Get mutable configuration
    pub fn config_mut(&mut self) -> &mut RepriseConfig {
        &mut self.config
    }

    /// List available profiles
    pub fn list_profiles(&self) -> Vec<&str> {
        self.profiles.keys().map(|s| s.as_str()).collect()
    }

    /// Get profile description
    pub fn profile_description(&self, name: &str) -> Option<&str> {
        self.profiles.get(name).map(|p| p.description.as_str())
    }

    /// Save current configuration to file
    pub fn save_to_file<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let path = path.as_ref();
        let content = match path.extension().and_then(|ext| ext.to_str()) {
            Some("toml") => toml::to_string_pretty(&self.config)
                .map_err(|e| RepriseError::config(format!("TOML serialize error: {}", e)))?,
            Some("yaml") | Some("yml") => serde_yaml::to_string(&self.config)
                .map_err(|e| RepriseError::config(format!("YAML serialize error: {}", e)))?,
            Some("json") => serde_json::to_string_pretty(&self.config)
                .map_err(|e| RepriseError::config(format!("JSON serialize error: {}", e)))?,
            _ => {
                return Err(RepriseError::config(
                    "Unsupported config file format. Use .toml, .yaml, .yml, or .json".to_string(),
                ))
            }
        };

        std::fs::write(path, content)
            .map_err(|e| RepriseError::io_error(format!("Failed to write config file: {}", e)))?;

        Ok(())
    }

    /// Convert configuration to pipeline config
    pub fn to_pipeline_config(&self) -> PipelineConfig {
        let settings = &self.config.pipeline;
        PipelineConfig {
            num_workers: if settings.num_workers == 0 {
                rayon::current_num_threads()
            } else {
                settings.num_workers
            },
            channel_capacity: settings.channel_capacity,
            batch_size: settings.batch_size,
            max_memory_usage: settings.max_memory_usage,
            min_frequency: settings.min_frequency,
            max_frequency: settings.max_frequency,
            region_extension: settings.region_extension,
            max_region_size: settings.max_region_size,
            enable_backpressure: settings.enable_backpressure,
            backpressure_threshold: settings.backpressure_threshold,
            min_alignment_score: settings.min_alignment_score,
            min_identity: settings.min_identity,
            alignment_config: crate::alignment::RepeatAlignmentConfig::default(),
        }
    }

    /// Load built-in configuration profiles
    fn load_builtin_profiles(&mut self) {
        // Fast mode: prioritize speed over accuracy
        let fast_profile = ConfigProfile {
            name: "fast".to_string(),
            description: "Fast processing mode with reduced accuracy".to_string(),
            config: RepriseConfig {
                pipeline: PipelineSettings {
                    min_frequency: 5,
                    max_frequency: Some(5000),
                    region_extension: 50,
                    min_alignment_score: 5,
                    min_identity: 0.40,
                    ..Default::default()
                },
                index: IndexSettings {
                    parallel: true,
                    max_positions_per_kmer: 5000,
                    ..Default::default()
                },
                ..Default::default()
            },
        };

        // Accurate mode: prioritize accuracy over speed
        let accurate_profile = ConfigProfile {
            name: "accurate".to_string(),
            description: "High accuracy mode with slower processing".to_string(),
            config: RepriseConfig {
                pipeline: PipelineSettings {
                    min_frequency: 2,
                    max_frequency: Some(20000),
                    region_extension: 200,
                    min_alignment_score: 20,
                    min_identity: 0.65,
                    ..Default::default()
                },
                index: IndexSettings {
                    parallel: true,
                    max_positions_per_kmer: 20000,
                    ..Default::default()
                },
                ..Default::default()
            },
        };

        // Large genome mode: optimized for genomes >1GB
        let large_genome_profile = ConfigProfile {
            name: "large_genome".to_string(),
            description: "Optimized for large genomes (>1GB)".to_string(),
            config: RepriseConfig {
                pipeline: PipelineSettings {
                    channel_capacity: 50000,
                    batch_size: 500,
                    max_memory_usage: 2 * 1024 * 1024 * 1024, // 2GB
                    min_frequency: 10,
                    max_frequency: Some(50000),
                    region_extension: 100,
                    enable_backpressure: true,
                    backpressure_threshold: 0.7,
                    ..Default::default()
                },
                index: IndexSettings {
                    parallel: true,
                    max_positions_per_kmer: 10000,
                    expected_memory_usage: Some(1024 * 1024 * 1024), // 1GB
                    ..Default::default()
                },
                system: SystemSettings {
                    memory_warning_threshold: 70.0,
                    memory_critical_threshold: 85.0,
                    enable_gc_hints: true,
                    ..Default::default()
                },
                recovery: RecoverySettings {
                    enable_checkpoints: true,
                    checkpoint_interval: Duration::from_secs(600), // 10 minutes
                    ..Default::default()
                },
                ..Default::default()
            },
        };

        // Small genome mode: optimized for genomes <100MB
        let small_genome_profile = ConfigProfile {
            name: "small_genome".to_string(),
            description: "Optimized for small genomes (<100MB)".to_string(),
            config: RepriseConfig {
                pipeline: PipelineSettings {
                    num_workers: (rayon::current_num_threads() / 2).max(1),
                    channel_capacity: 5000,
                    batch_size: 50,
                    max_memory_usage: 128 * 1024 * 1024, // 128MB
                    min_frequency: 2,
                    max_frequency: Some(5000),
                    region_extension: 150,
                    ..Default::default()
                },
                index: IndexSettings {
                    parallel: false,
                    max_positions_per_kmer: 5000,
                    ..Default::default()
                },
                ..Default::default()
            },
        };

        self.profiles.insert(fast_profile.name.clone(), fast_profile);
        self.profiles.insert(accurate_profile.name.clone(), accurate_profile);
        self.profiles.insert(large_genome_profile.name.clone(), large_genome_profile);
        self.profiles.insert(small_genome_profile.name.clone(), small_genome_profile);
    }
}

impl Default for ConfigManager {
    fn default() -> Self {
        Self::new()
    }
}

/// Custom serde module for Duration serialization
mod duration_serde {
    use serde::{Deserialize, Deserializer, Serializer};
    use std::time::Duration;

    pub fn serialize<S>(duration: &Duration, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_u64(duration.as_secs())
    }

    pub fn deserialize<'de, D>(deserializer: D) -> Result<Duration, D::Error>
    where
        D: Deserializer<'de>,
    {
        let secs = u64::deserialize(deserializer)?;
        Ok(Duration::from_secs(secs))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_default_config() {
        let config = RepriseConfig::default();
        assert!(config.validate().is_ok());
    }

    #[test]
    fn test_config_serialization() {
        let config = RepriseConfig::default();
        
        // Test TOML
        let toml_str = toml::to_string(&config).unwrap();
        let _: RepriseConfig = toml::from_str(&toml_str).unwrap();

        // Test JSON
        let json_str = serde_json::to_string(&config).unwrap();
        let _: RepriseConfig = serde_json::from_str(&json_str).unwrap();

        // Test YAML
        let yaml_str = serde_yaml::to_string(&config).unwrap();
        let _: RepriseConfig = serde_yaml::from_str(&yaml_str).unwrap();
    }

    #[test]
    fn test_config_manager_profiles() {
        let manager = ConfigManager::new();
        let profiles = manager.list_profiles();
        
        assert!(profiles.contains(&"fast"));
        assert!(profiles.contains(&"accurate"));
        assert!(profiles.contains(&"large_genome"));
        assert!(profiles.contains(&"small_genome"));

        assert!(manager.profile_description("fast").is_some());
        assert!(manager.profile_description("nonexistent").is_none());
    }

    #[test]
    fn test_config_file_loading() {
        let config = RepriseConfig::default();
        
        // Test TOML file
        let mut toml_file = NamedTempFile::new().unwrap();
        let toml_content = toml::to_string_pretty(&config).unwrap();
        toml_file.write_all(toml_content.as_bytes()).unwrap();
        toml_file.flush().unwrap();

        let manager = ConfigManager::load_from_file(toml_file.path()).unwrap();
        assert!(manager.config().validate().is_ok());
    }

    #[test]
    fn test_profile_application() {
        let mut manager = ConfigManager::new();
        
        // Apply fast profile
        manager.apply_profile("fast").unwrap();
        assert_eq!(manager.config().pipeline.min_frequency, 5);
        assert_eq!(manager.config().pipeline.min_identity, 0.40);

        // Apply accurate profile
        manager.apply_profile("accurate").unwrap();
        assert_eq!(manager.config().pipeline.min_frequency, 2);
        assert_eq!(manager.config().pipeline.min_identity, 0.65);
    }

    #[test]
    fn test_environment_variable_override() {
        let mut manager = ConfigManager::new();
        
        // Set environment variables
        env::set_var("REPRISE_WORKERS", "8");
        env::set_var("REPRISE_MAX_MEMORY", "1073741824"); // 1GB
        env::set_var("REPRISE_LOG_LEVEL", "debug");
        env::set_var("REPRISE_JSON_LOGS", "true");
        
        manager.load_from_env().unwrap();
        
        assert_eq!(manager.config().pipeline.num_workers, 8);
        assert_eq!(manager.config().pipeline.max_memory_usage, 1073741824);
        assert_eq!(manager.config().logging.level, crate::logging::LogLevel::Debug);
        assert_eq!(manager.config().logging.json_format, true);
        
        // Clean up environment variables
        env::remove_var("REPRISE_WORKERS");
        env::remove_var("REPRISE_MAX_MEMORY");
        env::remove_var("REPRISE_LOG_LEVEL");
        env::remove_var("REPRISE_JSON_LOGS");
    }

    #[test]
    fn test_pipeline_config_conversion() {
        let manager = ConfigManager::new();
        let pipeline_config = manager.to_pipeline_config();
        
        assert_eq!(pipeline_config.min_frequency, manager.config().pipeline.min_frequency);
        assert_eq!(pipeline_config.max_frequency, manager.config().pipeline.max_frequency);
        assert_eq!(pipeline_config.region_extension, manager.config().pipeline.region_extension);
    }

    #[test]
    fn test_config_validation() {
        let mut config = RepriseConfig::default();
        
        // Valid configuration should pass
        assert!(config.validate().is_ok());
        
        // Invalid worker count should fail
        config.pipeline.num_workers = 1000;
        assert!(config.validate().is_err());
        
        // Fix and test memory settings
        config.pipeline.num_workers = 4;
        config.pipeline.max_memory_usage = 50; // Too small
        assert!(config.validate().is_err());
        
        // Fix and test identity range
        config.pipeline.max_memory_usage = 256 * 1024 * 1024;
        config.pipeline.min_identity = 1.5; // Invalid range
        assert!(config.validate().is_err());
    }
}