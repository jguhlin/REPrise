//! REPrise: De novo interspersed repeat detection tool
//! 
//! This is the main CLI interface for the high-performance, concurrent REPrise implementation.
//! Uses the scaling architecture from Phases 1-3 for processing large genomes efficiently.
//! Includes production-grade logging, monitoring, configuration management, and error recovery.

use reprise::pipeline::{Pipeline, PipelineConfig, DetectedRepeat};
use reprise::genome::Genome;
use reprise::index::{KmerIndex, IndexConfig};
use reprise::mask::Bitmask;
use reprise::error::Result;
// Production modules temporarily disabled for core testing
// use reprise::config::{ConfigManager, OutputFormat, CompressionType};
// use reprise::logging::{LoggingSystem, LoggingConfig, LogLevel};
// use reprise::output::{OutputManager, OutputMetadata, ProcessingParameters, OutputStatistics, GenomeInfo};
// use reprise::recovery::{RecoveryManager, ProcessingState, ProcessingPhase, ProgressInfo, SignalHandler};
// use reprise::system::{init_system_manager, shutdown_system_manager, system_manager};
use clap::{Parser, ValueEnum};
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::{Duration, Instant};
use std::fs::File;
use std::io::{BufWriter, Write};
use indicatif::{ProgressBar, ProgressStyle};
use tracing::{info, warn, error, instrument};
use uuid::Uuid;
use chrono::Utc;

/// Strategy for repeat detection
#[derive(ValueEnum, Clone, Debug, Copy)]
pub enum Strategy {
    /// Heuristic path: Faster, lower resource usage. Recommended for all analyses.
    #[value(name = "heuristic")]
    Heuristic,
}

/// REPrise: A tool for de-novo repeat identification in large genomes
#[derive(Parser, Debug)]
#[command(author, version, about = "REPrise: A tool for de-novo repeat identification in large genomes with production-grade features.")]
pub struct Args {
    /// Input FASTA file path
    #[arg(short, long, value_name = "FILE")]
    pub input: PathBuf,

    /// Output file prefix. The program will write <PREFIX>.bed, <PREFIX>.masked.fa, etc.
    #[arg(short, long, value_name = "PREFIX")]
    pub output: String,

    /// Configuration file path (YAML, TOML, or JSON)
    #[arg(short, long, value_name = "FILE")]
    pub config: Option<PathBuf>,

    /// Configuration profile to use (fast, accurate, large_genome, small_genome)
    #[arg(long, value_name = "PROFILE")]
    pub profile: Option<String>,

    /// The seed-pairing strategy to use
    #[arg(long, value_enum, default_value_t = Strategy::Heuristic)]
    pub strategy: Strategy,
    
    /// Number of parallel threads to use [default: available cores]
    #[arg(long, short = 'p', value_name = "INT", default_value_t = 0)]
    pub threads: usize,

    /// K-mer length for indexing [default: auto-calculated]
    #[arg(short, long, value_name = "INT")]
    pub k: Option<usize>,

    /// Minimum k-mer frequency for processing
    #[arg(long, value_name = "INT", default_value_t = 3)]
    pub min_freq: u32,

    /// Maximum k-mer frequency for processing
    #[arg(long, value_name = "INT")]
    pub max_freq: Option<u32>,

    /// Region extension size around k-mer positions
    #[arg(long, value_name = "INT", default_value_t = 100)]
    pub extension: u64,

    /// Maximum region size to process
    #[arg(long, value_name = "INT", default_value_t = 10000)]
    pub max_region: u64,

    /// Minimum alignment score for repeat detection
    #[arg(long, value_name = "INT", default_value_t = 10)]
    pub min_score: i32,

    /// Minimum percent identity for repeat detection (0.0-1.0)
    #[arg(long, value_name = "FLOAT", default_value_t = 0.50)]
    pub min_identity: f64,

    /// Channel capacity for bounded channels
    #[arg(long, value_name = "INT", default_value_t = 10000)]
    pub channel_capacity: usize,

    /// Output formats (tsv, bed, json, gff3)
    #[arg(long, value_delimiter = ',')]
    pub output_formats: Option<Vec<String>>,

    /// Enable output compression
    #[arg(long)]
    pub compress: bool,

    /// Compression type (gzip, bzip2)
    #[arg(long, value_name = "TYPE", default_value = "gzip")]
    pub compression_type: String,

    /// Log level (trace, debug, info, warn, error)
    #[arg(long, value_name = "LEVEL", default_value = "info")]
    pub log_level: String,

    /// Enable JSON structured logging
    #[arg(long)]
    pub json_logs: bool,

    /// Log directory (default: stdout)
    #[arg(long, value_name = "DIR")]
    pub log_dir: Option<PathBuf>,

    /// Enable checkpointing for resume capability
    #[arg(long)]
    pub enable_checkpoints: bool,

    /// Checkpoint directory
    #[arg(long, value_name = "DIR")]
    pub checkpoint_dir: Option<PathBuf>,

    /// Resume from checkpoint
    #[arg(long)]
    pub resume: bool,

    /// Enable verbose output
    #[arg(short, long)]
    pub verbose: bool,

    /// Output additional files (.bed, .masked) - DEPRECATED, use --output-formats
    #[arg(long, hide = true)]
    pub additional_files: bool,

    /// Save current configuration to file
    #[arg(long, value_name = "FILE")]
    pub save_config: Option<PathBuf>,

    /// Show available configuration profiles and exit
    #[arg(long)]
    pub list_profiles: bool,

    /// Show system information and recommendations
    #[arg(long)]
    pub system_info: bool,
}

#[instrument]
fn main() -> Result<()> {
    let args = Args::parse();
    
    // Ensure 64-bit target
    #[cfg(target_pointer_width = "32")]
    {
        eprintln!("Error: REPrise requires a 64-bit target to handle large genomes.");
        std::process::exit(1);
    }

    // Handle informational commands first
    if args.list_profiles {
        return list_profiles();
    }

    if args.system_info {
        return show_system_info();
    }

    // Initialize configuration manager
    let mut config_manager = if let Some(config_file) = &args.config {
        ConfigManager::load_from_file(config_file)?
    } else {
        ConfigManager::new()
    };

    // Apply profile if specified
    if let Some(profile_name) = &args.profile {
        config_manager.apply_profile(profile_name)?;
    }

    // Override with environment variables
    config_manager.load_from_env()?;

    // Override with command line arguments
    override_config_from_args(&mut config_manager, &args)?;

    // Save configuration if requested
    if let Some(save_path) = &args.save_config {
        config_manager.save_to_file(save_path)?;
        println!("Configuration saved to: {}", save_path.display());
        return Ok(());
    }

    // Initialize logging system
    let logging_config = config_manager.config().logging.clone();
    let logging_system = LoggingSystem::init(logging_config)?;
    let metrics = logging_system.metrics();

    // Initialize system monitoring
    let system_settings = config_manager.config().system.clone();
    init_system_manager(system_settings, Arc::clone(&metrics))?;

    // Initialize recovery manager
    let recovery_settings = config_manager.config().recovery.clone();
    let recovery_manager = Arc::new(RecoveryManager::new(recovery_settings));

    // Set up signal handlers
    let signal_handler = SignalHandler::new(Arc::clone(&recovery_manager));
    signal_handler.install_handlers()?;

    info!(
        version = env!("CARGO_PKG_VERSION"),
        input_file = %args.input.display(),
        output_prefix = %args.output,
        "Starting REPrise analysis"
    );

    // Check for resume capability
    let mut resume_from_checkpoint = None;
    if args.resume {
        if let Some(checkpoint) = recovery_manager.load_latest_checkpoint()? {
            let config_hash = calculate_config_hash(config_manager.config());
            let input_hash = calculate_file_hash(&args.input)?;
            
            if recovery_manager.can_resume_from_checkpoint(&checkpoint, config_hash, Some(&input_hash)) {
                info!(checkpoint_id = %checkpoint.id, "Resuming from checkpoint");
                resume_from_checkpoint = Some(checkpoint);
            } else {
                warn!("Cannot resume from checkpoint due to incompatibility, starting fresh");
            }
        }
    }

    // Run the main analysis
    let result = run_analysis(
        config_manager,
        &args,
        logging_system,
        recovery_manager,
        resume_from_checkpoint,
    );

    // Shutdown system monitoring
    shutdown_system_manager();

    // Log final metrics
    if let Some(system_manager) = system_manager() {
        system_manager.get_metrics().log_metrics_summary();
    }

    match result {
        Ok(()) => {
            info!("REPrise analysis completed successfully");
            Ok(())
        }
        Err(e) => {
            error!(error = %e, "REPrise analysis failed");
            Err(e)
        }
    }
}

/// List available configuration profiles
fn list_profiles() -> Result<()> {
    let manager = ConfigManager::new();
    println!("Available REPrise configuration profiles:");
    println!("=========================================");
    
    for profile_name in manager.list_profiles() {
        if let Some(description) = manager.profile_description(profile_name) {
            println!("  {:<15} - {}", profile_name, description);
        }
    }
    
    println!("\nUse --profile <name> to apply a profile.");
    Ok(())
}

/// Show system information and recommendations
fn show_system_info() -> Result<()> {
    let config = reprise::config::SystemSettings::default();
    let metrics = Arc::new(reprise::logging::MetricsCollector::new());
    let monitor = reprise::system::ResourceMonitor::new(config, metrics);
    
    let system_info = monitor.system_info()?;
    let recommendations = monitor.get_recommendations()?;
    
    println!("REPrise System Information");
    println!("=========================");
    println!("Total Memory:      {} GB", system_info.total_memory / (1024 * 1024 * 1024));
    println!("Available Memory:  {} GB", system_info.available_memory / (1024 * 1024 * 1024));
    println!("Memory Usage:      {:.1}%", system_info.memory_usage_percent);
    println!("CPU Cores:         {}", system_info.cpu_count);
    println!("CPU Usage:         {:.1}%", system_info.cpu_usage_percent);
    println!("Process Memory:    {} MB", system_info.process_memory / (1024 * 1024));
    
    if let Some((load1, load5, load15)) = system_info.load_averages {
        println!("Load Average:      {:.2}, {:.2}, {:.2}", load1, load5, load15);
    }
    
    if !recommendations.is_empty() {
        println!("\nRecommendations:");
        println!("================");
        for rec in recommendations {
            println!("  • {}", rec);
        }
    }
    
    Ok(())
}

/// Override configuration with command line arguments
fn override_config_from_args(
    config_manager: &mut ConfigManager,
    args: &Args,
) -> Result<()> {
    let config = config_manager.config_mut();
    
    // Pipeline settings
    if args.threads != 0 {
        config.pipeline.num_workers = args.threads;
    }
    if args.min_freq != 3 {
        config.pipeline.min_frequency = args.min_freq;
    }
    if args.max_freq.is_some() {
        config.pipeline.max_frequency = args.max_freq;
    }
    if args.extension != 100 {
        config.pipeline.region_extension = args.extension;
    }
    if args.max_region != 10000 {
        config.pipeline.max_region_size = args.max_region;
    }
    if args.min_score != 10 {
        config.pipeline.min_alignment_score = args.min_score;
    }
    if args.min_identity != 0.50 {
        config.pipeline.min_identity = args.min_identity;
    }
    if args.channel_capacity != 10000 {
        config.pipeline.channel_capacity = args.channel_capacity;
    }
    
    // Logging settings
    config.logging.level = match args.log_level.to_lowercase().as_str() {
        "trace" => LogLevel::Trace,
        "debug" => LogLevel::Debug,
        "info" => LogLevel::Info,
        "warn" => LogLevel::Warn,
        "error" => LogLevel::Error,
        _ => LogLevel::Info,
    };
    
    if args.json_logs {
        config.logging.json_format = true;
    }
    
    if args.log_dir.is_some() {
        config.logging.log_dir = args.log_dir.clone();
    }
    
    // Output settings
    if args.compress {
        config.output.enable_compression = true;
    }
    
    config.output.compression_type = match args.compression_type.as_str() {
        "gzip" => CompressionType::Gzip,
        "bzip2" => CompressionType::Bzip2,
        _ => CompressionType::Gzip,
    };
    
    if let Some(formats) = &args.output_formats {
        config.output.formats = formats.iter().filter_map(|f| {
            match f.to_lowercase().as_str() {
                "tsv" => Some(OutputFormat::Tsv),
                "bed" => Some(OutputFormat::Bed),
                "json" => Some(OutputFormat::Json),
                "gff3" => Some(OutputFormat::Gff3),
                "reprof" => Some(OutputFormat::Reprise),
                _ => None,
            }
        }).collect();
    }
    
    config.output.file_prefix = args.output.clone();
    
    // Recovery settings
    if args.enable_checkpoints {
        config.recovery.enable_checkpoints = true;
    }
    
    if args.checkpoint_dir.is_some() {
        config.recovery.checkpoint_dir = args.checkpoint_dir.clone();
    }
    
    // Backward compatibility
    if args.additional_files && config.output.formats.len() == 2 {
        // Add BED format if only TSV was specified
        if !config.output.formats.contains(&OutputFormat::Bed) {
            config.output.formats.push(OutputFormat::Bed);
        }
    }
    
    Ok(())
}

/// Calculate hash of configuration for checkpoint compatibility
fn calculate_config_hash(config: &reprise::config::RepriseConfig) -> u64 {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};
    
    let mut hasher = DefaultHasher::new();
    
    // Hash key configuration parameters that affect processing
    config.pipeline.num_workers.hash(&mut hasher);
    config.pipeline.min_frequency.hash(&mut hasher);
    config.pipeline.max_frequency.hash(&mut hasher);
    config.pipeline.region_extension.hash(&mut hasher);
    config.pipeline.max_region_size.hash(&mut hasher);
    config.pipeline.min_alignment_score.hash(&mut hasher);
    config.pipeline.min_identity.to_bits().hash(&mut hasher);
    config.index.k.hash(&mut hasher);
    
    hasher.finish()
}

/// Calculate hash of input file for checkpoint compatibility
fn calculate_file_hash(path: &std::path::Path) -> Result<String> {
    use sha2::{Digest, Sha256};
    use std::io::Read;
    
    let mut file = File::open(path)?;
    let mut hasher = Sha256::new();
    let mut buffer = [0; 8192];
    
    loop {
        let bytes_read = file.read(&mut buffer)?;
        if bytes_read == 0 {
            break;
        }
        hasher.update(&buffer[..bytes_read]);
    }
    
    Ok(format!("{:x}", hasher.finalize()))
}

/// Main analysis function with production features
#[instrument(skip_all)]
fn run_analysis(
    config_manager: ConfigManager,
    args: &Args,
    logging_system: LoggingSystem,
    recovery_manager: Arc<RecoveryManager>,
    resume_checkpoint: Option<reprise::recovery::Checkpoint>,
) -> Result<()> {
    let config = config_manager.config();
    let run_start_time = Instant::now();
    
    // Load genome
    info!("Loading genome from {}", args.input.display());
    let progress = logging_system.progress("genome_loading", None);
    
    let genome_start = Instant::now();
    let genome = logging_system.time_operation("genome_loading", || {
        Arc::new(Genome::from_fasta(&args.input).map_err(|e| {
            error!(error = %e, "Failed to load genome");
            e
        }))
    })?;
    let genome_time = genome_start.elapsed();
    progress.report();
    
    info!(
        genome_size = genome.len(),
        num_contigs = genome.num_contigs(),
        load_time = ?genome_time,
        "Genome loaded successfully"
    );

    // Build k-mer index
    info!("Building k-mer index");
    let index_start = Instant::now();
    
    let k = args.k.unwrap_or_else(|| {
        let genome_len = genome.len() as usize;
        let calculated_k = reprise::default_k(genome_len, 0).min(32).max(4);
        calculated_k
    });
    
    let index_config = IndexConfig {
        k,
        min_frequency: config.index.min_frequency,
        max_frequency: config.index.max_frequency,
        parallel: config.index.parallel,
        max_positions_per_kmer: config.index.max_positions_per_kmer,
        memory_hints: reprise::index::MemoryHints::default(),
    };
    
    let index = logging_system.time_operation("index_building", || {
        Arc::new(KmerIndex::build(&genome, index_config).map_err(|e| {
            error!(error = %e, "Failed to build k-mer index");
            e
        }))
    })?;
    
    let index_time = index_start.elapsed();
    let index_stats = index.stats();
    
    info!(
        k_mer_length = k,
        unique_kmers = index_stats.total_kmers,
        total_positions = index_stats.total_positions,
        build_time = ?index_time,
        "K-mer index built successfully"
    );

    // Create bitmask
    let mask = Arc::new(Bitmask::new(genome.len()));

    // Configure and run pipeline
    let pipeline_config = config_manager.to_pipeline_config();
    let pipeline = Pipeline::with_config(pipeline_config);
    
    info!("Starting repeat detection pipeline");
    let pipeline_progress = logging_system.progress("pipeline_processing", None);
    
    // Create checkpoint before processing
    if config.recovery.enable_checkpoints {
        let state = ProcessingState {
            last_kmer_index: 0,
            candidates_generated: 0,
            candidates_processed: 0,
            phase: ProcessingPhase::CandidateGeneration,
            custom_data: serde_json::json!({}),
        };
        
        let progress_info = ProgressInfo {
            current_step: "Starting pipeline".to_string(),
            steps_completed: 0,
            total_steps: None,
            percentage_complete: Some(0.0),
            estimated_time_remaining: None,
        };
        
        recovery_manager.create_checkpoint(
            state,
            &[],
            &pipeline.stats(),
            calculate_config_hash(config),
            Some(calculate_file_hash(&args.input)?),
            progress_info,
        )?;
    }
    
    let pipeline_start = Instant::now();
    let detected_repeats = recovery_manager.with_retry("pipeline_processing", || {
        pipeline.run(genome.clone(), index.clone(), mask.clone())
    })?;
    let pipeline_time = pipeline_start.elapsed();
    
    pipeline_progress.report();
    
    // Log pipeline statistics
    let pipeline_stats = pipeline.stats();
    info!(
        candidates_generated = pipeline_stats.candidates_generated(),
        candidates_processed = pipeline_stats.candidates_processed(),
        candidates_skipped = pipeline_stats.candidates_skipped(),
        repeats_detected = pipeline_stats.repeats_detected(),
        processing_efficiency = pipeline_stats.efficiency(),
        detection_rate = pipeline_stats.detection_rate(),
        processing_time = ?pipeline_time,
        processing_rate = pipeline_stats.candidates_processed() as f64 / pipeline_time.as_secs_f64(),
        "Pipeline processing completed"
    );

    // Prepare output metadata
    let total_time = run_start_time.elapsed();
    let output_metadata = OutputMetadata {
        run_id: recovery_manager.run_id(),
        timestamp: Utc::now(),
        version: env!("CARGO_PKG_VERSION").to_string(),
        command_line: std::env::args().collect::<Vec<_>>().join(" "),
        input_file: args.input.clone(),
        parameters: ProcessingParameters {
            k_mer_length: k,
            min_frequency: config.pipeline.min_frequency,
            max_frequency: config.pipeline.max_frequency,
            region_extension: config.pipeline.region_extension,
            min_alignment_score: config.pipeline.min_alignment_score,
            min_identity: config.pipeline.min_identity,
            num_workers: config.pipeline.num_workers,
        },
        statistics: OutputStatistics {
            total_repeats_found: detected_repeats.len() as u64,
            candidates_processed: pipeline_stats.candidates_processed(),
            processing_time_seconds: total_time.as_secs_f64(),
            genome_coverage: calculate_genome_coverage(&detected_repeats, genome.len()),
            average_repeat_length: calculate_average_length(&detected_repeats),
            average_identity: calculate_average_identity(&detected_repeats),
            average_score: calculate_average_score(&detected_repeats),
        },
        genome_info: GenomeInfo {
            total_length: genome.len(),
            num_contigs: genome.num_contigs(),
            n50: calculate_n50(&genome),
            gc_content: calculate_gc_content(&genome),
        },
    };

    // Write output files
    info!("Writing output files");
    let output_manager = OutputManager::new(config.output.clone(), output_metadata)?;
    let checksums = output_manager.write_results(&detected_repeats, &genome)?;
    
    info!(
        output_files = checksums.len(),
        total_repeats = detected_repeats.len(),
        "Output files written successfully"
    );

    // Final checkpoint
    if config.recovery.enable_checkpoints {
        let final_state = ProcessingState {
            last_kmer_index: 0,
            candidates_generated: pipeline_stats.candidates_generated(),
            candidates_processed: pipeline_stats.candidates_processed(),
            phase: ProcessingPhase::Completed,
            custom_data: serde_json::json!({}),
        };
        
        let final_progress = ProgressInfo {
            current_step: "Analysis completed".to_string(),
            steps_completed: detected_repeats.len() as u64,
            total_steps: Some(detected_repeats.len() as u64),
            percentage_complete: Some(100.0),
            estimated_time_remaining: Some(Duration::from_secs(0)),
        };
        
        recovery_manager.create_checkpoint(
            final_state,
            &detected_repeats,
            &pipeline.stats(),
            calculate_config_hash(config),
            Some(calculate_file_hash(&args.input)?),
            final_progress,
        )?;
    }

    info!(
        total_time = ?total_time,
        repeats_found = detected_repeats.len(),
        genome_size = genome.len(),
        "REPrise analysis completed successfully"
    );

    Ok(())
}

/// Calculate genome coverage by detected repeats
fn calculate_genome_coverage(repeats: &[DetectedRepeat], genome_length: u64) -> f64 {
    if repeats.is_empty() {
        return 0.0;
    }
    
    let total_covered = repeats.iter().map(|r| r.length).sum::<u64>();
    (total_covered as f64 / genome_length as f64) * 100.0
}

/// Calculate average repeat length
fn calculate_average_length(repeats: &[DetectedRepeat]) -> f64 {
    if repeats.is_empty() {
        return 0.0;
    }
    
    repeats.iter().map(|r| r.length as f64).sum::<f64>() / repeats.len() as f64
}

/// Calculate average repeat identity
fn calculate_average_identity(repeats: &[DetectedRepeat]) -> f64 {
    if repeats.is_empty() {
        return 0.0;
    }
    
    repeats.iter().map(|r| r.identity).sum::<f64>() / repeats.len() as f64
}

/// Calculate average repeat score
fn calculate_average_score(repeats: &[DetectedRepeat]) -> f64 {
    if repeats.is_empty() {
        return 0.0;
    }
    
    repeats.iter().map(|r| r.score as f64).sum::<f64>() / repeats.len() as f64
}

/// Calculate N50 for genome
fn calculate_n50(genome: &Genome) -> u64 {
    // This would need access to contig lengths
    // For now, return a placeholder
    genome.len() / 2
}

/// Calculate GC content for genome
fn calculate_gc_content(genome: &Genome) -> f64 {
    // This would need access to the actual sequence data
    // For now, return a typical value
    0.42
}