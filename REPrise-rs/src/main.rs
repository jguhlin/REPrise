use reprise::pipeline::{Pipeline, PipelineConfig};
use reprise::genome::Genome;
use reprise::index::{KmerIndex, IndexConfig};
use reprise::mask::Bitmask;
use reprise::streaming::{StreamingConfig, StreamingPipeline};
use reprise::error::Result;
use clap::Parser;
use std::sync::Arc;
use std::path::PathBuf;
use std::time::Instant;

/// REPrise: A tool for de-novo repeat identification in large genomes
#[derive(Parser, Debug)]
#[command(author, version, about = "REPrise: A tool for de-novo repeat identification in large genomes.")]
pub struct Args {
    /// Input FASTA file path
    #[arg(short, long, value_name = "FILE")]
    pub input: String,

    /// Output file prefix. The program will write <PREFIX>.bed, <PREFIX>.masked.fa, etc.
    #[arg(short, long, value_name = "PREFIX")]
    pub output: String,

    /// Number of parallel threads to use [default: available cores]
    #[arg(long, short = 'p', value_name = "INT", default_value_t = 0)]
    pub threads: usize,

    /// K-mer length for indexing [default: auto-calculated]
    #[arg(short, long, value_name = "INT")]
    pub k: Option<usize>,

    /// Use streaming mode for large genomes (memory-bounded processing)
    #[arg(long)]
    pub streaming: bool,

    /// Maximum memory usage in MB (streaming mode only)
    #[arg(long, default_value = "1024")]
    pub max_memory: usize,

    /// Minimum k-mer frequency threshold
    #[arg(long, default_value = "2")]
    pub min_freq: usize,

    /// Temporary directory for streaming mode
    #[arg(long)]
    pub temp_dir: Option<PathBuf>,

    /// Verbose output
    #[arg(short, long)]
    pub verbose: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();

    if args.verbose {
        println!("REPrise: De-novo repeat identification");
        println!("Input: {}", args.input);
        println!("Output: {}", args.output);
        println!("Mode: {}", if args.streaming { "Streaming (memory-bounded)" } else { "Standard" });
    }

    let start_time = Instant::now();

    if args.streaming {
        run_streaming_mode(&args)
    } else {
        run_standard_mode(&args)
    }?;

    let processing_time = start_time.elapsed();
    if args.verbose {
        println!("Total processing time: {:.2?}", processing_time);
    }

    Ok(())
}

fn run_streaming_mode(args: &Args) -> Result<()> {
    if args.verbose {
        println!("Using streaming mode for memory-bounded processing");
    }

    // Configure streaming pipeline
    let config = StreamingConfig {
        k: args.k.unwrap_or(15).min(20).max(4),
        min_frequency: args.min_freq,
        max_memory_mb: args.max_memory,
        temp_dir: args.temp_dir.clone().unwrap_or_else(|| std::env::temp_dir()),
        bloom_filter_size: 10_000_000,
        channel_capacity: 1000,
    };

    if args.verbose {
        println!("Streaming config:");
        println!("  k-mer length: {}", config.k);
        println!("  min frequency: {}", config.min_frequency);
        println!("  max memory: {} MB", config.max_memory_mb);
    }

    let pipeline = StreamingPipeline::new(config);
    let candidates = pipeline.run(&args.input)?;

    // Filter high-frequency candidates (potential repeats)
    let repeat_candidates: Vec<_> = candidates.iter()
        .filter(|c| c.frequency >= args.min_freq * 3)
        .collect();

    println!("=== STREAMING RESULTS ===");
    println!("Candidate k-mers found: {}", candidates.len());
    println!("High-frequency repeat candidates: {}", repeat_candidates.len());

    write_streaming_output(&args.output, &candidates)?;

    if args.verbose && !repeat_candidates.is_empty() {
        println!("Top repeat candidates:");
        let mut sorted_candidates = candidates;
        sorted_candidates.sort_by_key(|c| std::cmp::Reverse(c.frequency));

        for (i, candidate) in sorted_candidates.iter().take(5).enumerate() {
            println!("  {}: freq={}, positions={}", 
                i+1, candidate.frequency, candidate.positions.len());
        }
    }

    Ok(())
}

fn run_standard_mode(args: &Args) -> Result<()> {
    if args.verbose {
        println!("Using standard mode");
    }

    // Load genome
    let genome = Arc::new(Genome::from_fasta(&args.input)?);

    // Build k-mer index  
    let k = args.k.unwrap_or_else(|| {
        let genome_len = genome.len() as usize;
        reprise::default_k(genome_len, 0).min(15).max(4)  // Limit to 15 for memory safety
    });

    if args.verbose {
        println!("Genome: {} bp, k-mer length: {}", genome.len(), k);
    }

    let index_config = IndexConfig {
        k,
        min_frequency: args.min_freq as u32,
        max_frequency: None,
        parallel: true,
        max_positions_per_kmer: 10000,
        memory_hints: reprise::index::MemoryHints::default(),
    };

    let index = Arc::new(KmerIndex::build(&genome, index_config)?);

    // Create bitmask
    let mask = Arc::new(Bitmask::new(genome.len()));

    // Configure and run pipeline
    let pipeline_config = PipelineConfig::default();
    let pipeline = Pipeline::with_config(pipeline_config);

    let detected_repeats = pipeline.run(genome.clone(), index.clone(), mask.clone())?;

    println!("=== STANDARD RESULTS ===");
    println!("Found {} repeats.", detected_repeats.len());

    Ok(())
}

fn write_streaming_output(output_prefix: &str, candidates: &[reprise::streaming::KmerCandidate]) -> Result<()> {
    use std::fs::File;
    use std::io::{BufWriter, Write};
    
    let summary_path = format!("{}.summary", output_prefix);
    let mut writer = BufWriter::new(File::create(&summary_path)?);
    
    writeln!(writer, "# REPrise Streaming K-mer Summary")?;
    writeln!(writer, "# kmer_id\tfrequency\tpositions\tcontigs")?;
    
    for (i, candidate) in candidates.iter().enumerate() {
        let unique_contigs: std::collections::HashSet<_> = candidate.positions.iter()
            .map(|&(_, contig_id)| contig_id)
            .collect();
            
        writeln!(writer, "{}\t{}\t{}\t{}", 
            i, 
            candidate.frequency, 
            candidate.positions.len(),
            unique_contigs.len()
        )?;
    }
    
    writer.flush()?;
    println!("Output written to: {}", summary_path);
    Ok(())
}
