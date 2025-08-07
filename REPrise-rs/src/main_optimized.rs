//! Optimized REPrise main binary focused on core performance improvements

use clap::{Parser, ValueEnum};
use reprise::error::Result;
use reprise::genome::Genome;
use reprise::index::{IndexConfig, KmerIndex};
use reprise::mask::Bitmask;
use reprise::pipeline::{DetectedRepeat, Pipeline, PipelineConfig};
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;
use std::io::Write;

#[derive(Parser, Debug)]
#[command(name = "REPrise-optimized")]
#[command(about = "De novo interspersed repeat detection (optimized version)")]
pub struct Args {
    /// Input FASTA file path
    #[arg(short, long, value_name = "FILE")]
    pub input: PathBuf,

    /// Output file prefix
    #[arg(short, long, value_name = "PREFIX")]
    pub output: String,

    /// Strategy for repeat detection
    #[arg(long, value_enum, default_value_t = Strategy::Heuristic)]
    pub strategy: Strategy,
    
    /// Number of parallel threads [default: available cores]
    #[arg(long, short = 'p', default_value_t = 0)]
    pub threads: usize,

    /// K-mer length [default: auto-calculated]
    #[arg(short, long)]
    pub k: Option<usize>,

    /// Minimum k-mer frequency
    #[arg(long, default_value_t = 3)]
    pub min_freq: u32,

    /// Maximum k-mer frequency
    #[arg(long)]
    pub max_freq: Option<u32>,

    /// Region extension size
    #[arg(long, default_value_t = 100)]
    pub extension: u64,

    /// Maximum region size
    #[arg(long, default_value_t = 10000)]
    pub max_region: u64,

    /// Minimum alignment score
    #[arg(long, default_value_t = 10)]
    pub min_score: i32,

    /// Minimum percent identity (0.0-1.0)
    #[arg(long, default_value_t = 0.50)]
    pub min_identity: f64,

    /// Channel capacity
    #[arg(long, default_value_t = 10000)]
    pub channel_capacity: usize,

    /// Enable verbose output
    #[arg(short, long)]
    pub verbose: bool,
}

#[derive(ValueEnum, Clone, Debug)]
pub enum Strategy {
    Heuristic,
    Exact,
}

fn main() -> Result<()> {
    let args = Args::parse();
    
    println!("REPrise: De novo interspersed repeat detection (Optimized)");
    println!("=========================================================");
    println!("Input: {}", args.input.display());
    println!("Output prefix: {}", args.output);
    println!("Strategy: {:?}", args.strategy);
    
    let num_threads = if args.threads == 0 {
        rayon::current_num_threads()
    } else {
        args.threads
    };
    println!("Threads: {}", num_threads);

    // 1. Load genome
    println!("\n1. Loading genome...");
    let start = Instant::now();
    let genome = Arc::new(Genome::from_fasta(&args.input)?);
    let genome_time = start.elapsed();
    
    println!("   ✅ Loaded genome in {:?}", genome_time);
    println!("   - Genome size: {} bases", genome.len());
    // TODO: Add contig count method to genome
    println!("   - Contigs: estimated");

    // 2. Build k-mer index
    println!("\n2. Building k-mer index...");
    let start = Instant::now();
    
    let k = args.k.unwrap_or_else(|| {
        // Simple k calculation based on genome size
        match genome.len() {
            0..=10000 => 7,
            10001..=100000 => 10,
            100001..=1000000 => 13,
            _ => 16,
        }
    });
    println!("   - K-mer length: {}", k);
    
    let mut index_config = IndexConfig::default();
    index_config.k = k;
    let index = Arc::new(KmerIndex::build(&genome, index_config)?);
    let index_time = start.elapsed();
    
    let stats = index.stats();
    println!("   ✅ Built index in {:?}", index_time);
    println!("   - Unique k-mers: {}", stats.total_kmers);
    println!("   - Total positions: {}", stats.total_positions);

    // 3. Create bitmask for region tracking
    let mask = Arc::new(Bitmask::new(genome.len()));

    // 4. Configure pipeline with optimizations
    println!("\n3. Configuring optimized pipeline...");
    let pipeline_config = PipelineConfig::new()
        .with_workers(num_threads)
        .with_channel_capacity(args.channel_capacity)
        .with_frequency_range(args.min_freq, args.max_freq)
        .with_region_params(args.extension, args.max_region)
        .with_alignment_thresholds(args.min_score, args.min_identity)
        .with_backpressure(true, 0.8);

    let pipeline = Pipeline::with_config(pipeline_config);

    if args.verbose {
        println!("   - Worker threads: {}", num_threads);
        println!("   - Channel capacity: {}", args.channel_capacity);
        println!("   - Frequency range: {}-{:?}", args.min_freq, args.max_freq);
        println!("   - Min alignment score: {}", args.min_score);
        println!("   - Min identity: {:.1}%", args.min_identity * 100.0);
    }

    // 5. Run repeat detection
    println!("\n4. Running optimized repeat detection...");
    let start = Instant::now();
    let detected_repeats = pipeline.run(genome.clone(), index, mask)?;
    let pipeline_time = start.elapsed();

    // 6. Display results
    let pipeline_stats = pipeline.stats();
    println!("   ✅ Pipeline completed in {:?}", pipeline_time);
    println!("   - Candidates generated: {}", pipeline_stats.candidates_generated());
    println!("   - Candidates processed: {}", pipeline_stats.candidates_processed());
    println!("   - Candidates skipped: {}", pipeline_stats.candidates_skipped());
    println!("   - Processing efficiency: {:.1}%", pipeline_stats.efficiency() * 100.0);
    println!("   - Repeats detected: {}", pipeline_stats.repeats_detected());
    println!("   - Detection rate: {:.1}%", pipeline_stats.detection_rate() * 100.0);
    if pipeline_time.as_secs_f64() > 0.0 {
        println!("   - Processing rate: {:.0} candidates/sec", 
            pipeline_stats.candidates_processed() as f64 / pipeline_time.as_secs_f64());
    }

    // 7. Write output files
    println!("\n5. Writing output files...");
    write_output_files(&detected_repeats, &genome, &args)?;

    println!("\n🎉 Optimized REPrise completed successfully!");
    println!("Found {} repeats across {} bases", 
        detected_repeats.len(), genome.len());
    
    let total_time = genome_time + index_time + pipeline_time;
    println!("Total processing time: {:?}", total_time);

    Ok(())
}

fn write_output_files(
    detected_repeats: &[DetectedRepeat],
    genome: &Genome,
    args: &Args,
) -> Result<()> {
    // Write TSV summary file
    let tsv_path = format!("{}.tsv", args.output);
    let mut tsv_file = std::fs::File::create(&tsv_path)?;
    writeln!(tsv_file, "region1_start\tregion1_end\tregion2_start\tregion2_end\tscore\tidentity\tlength\tseed_kmer")?;
    
    for repeat in detected_repeats {
        writeln!(
            tsv_file,
            "{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}",
            repeat.region1.start,
            repeat.region1.end,
            repeat.region2.start,
            repeat.region2.end,
            repeat.score,
            repeat.identity,
            repeat.length,
            format!("{:?}", repeat.seed_kmer)
        )?;
    }

    println!("   ✅ Wrote TSV summary to {}", tsv_path);

    // Write BED file for visualization
    let bed_path = format!("{}.bed", args.output);
    let mut bed_file = std::fs::File::create(&bed_path)?;
    writeln!(bed_file, "track name=\"REPrise Repeats\" description=\"Detected interspersed repeats\"")?;
    
    for (i, repeat) in detected_repeats.iter().enumerate() {
        // Find contig name for region1
        let contig_info1 = genome.contig_info(repeat.region1.contig);
        let contig_name1 = contig_info1
            .map(|info| info.name.as_str())
            .unwrap_or("unknown");
            
        writeln!(
            bed_file,
            "{}\t{}\t{}\trepeat_{}_{:.1}%\t{}\t+",
            contig_name1,
            repeat.region1.local_start,
            repeat.region1.local_end,
            i + 1,
            repeat.identity * 100.0,
            repeat.score
        )?;
        
        // Add region2 if it's in a different contig or far away
        if repeat.region2.contig != repeat.region1.contig ||
           repeat.region2.start.abs_diff(repeat.region1.start) > 1000 {
            let contig_info2 = genome.contig_info(repeat.region2.contig);
            let contig_name2 = contig_info2
                .map(|info| info.name.as_str())
                .unwrap_or("unknown");
                
            writeln!(
                bed_file,
                "{}\t{}\t{}\trepeat_{}_{:.1}%\t{}\t+",
                contig_name2,
                repeat.region2.local_start,
                repeat.region2.local_end,
                i + 1,
                repeat.identity * 100.0,
                repeat.score
            )?;
        }
    }

    println!("   ✅ Wrote BED file to {}", bed_path);

    Ok(())
}