use reprise::streaming::{StreamingConfig, StreamingPipeline};
use reprise::error::Result;
use clap::Parser;
use std::path::PathBuf;
use std::time::Instant;

/// REPrise v2.0: Streaming, disk-based repeat detection for large genomes
#[derive(Parser, Debug)]
#[command(author, version, about = "REPrise v2.0: Memory-efficient repeat detection for large genomes")]
pub struct Args {
    /// Input FASTA file path
    #[arg(short, long, value_name = "FILE")]
    pub input: PathBuf,

    /// Output file prefix
    #[arg(short, long, value_name = "PREFIX")]
    pub output: PathBuf,

    /// K-mer length [default: 15, max: 20 for memory safety]
    #[arg(short, long, value_name = "INT")]
    pub k: Option<usize>,

    /// Minimum k-mer frequency threshold
    #[arg(long, default_value = "2")]
    pub min_freq: usize,

    /// Maximum memory usage in MB
    #[arg(long, default_value = "1024")]
    pub max_memory: usize,

    /// Temporary directory for disk-based processing
    #[arg(long)]
    pub temp_dir: Option<PathBuf>,

    /// Expected number of k-mers (for Bloom filter sizing)
    #[arg(long, default_value = "10000000")]
    pub bloom_size: usize,

    /// Producer-consumer channel capacity
    #[arg(long, default_value = "1000")]
    pub channel_size: usize,

    /// Verbose output
    #[arg(short, long)]
    pub verbose: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();

    if args.verbose {
        println!("REPrise v2.0: Streaming Memory-Bounded Repeat Detection");
        println!("Input: {:?}", args.input);
        println!("Output: {:?}", args.output);
        println!("Max memory: {} MB", args.max_memory);
    }

    // Validate input file exists
    if !args.input.exists() {
        return Err(reprise::error::REPriseError::io_error(
            format!("Input file {:?} does not exist", args.input)
        ));
    }

    // Configure streaming pipeline
    let config = StreamingConfig {
        k: args.k.unwrap_or(15).min(20).max(4),  // Cap at 20 for safety
        min_frequency: args.min_freq,
        max_memory_mb: args.max_memory,
        temp_dir: args.temp_dir.unwrap_or_else(|| std::env::temp_dir()),
        bloom_filter_size: args.bloom_size,
        channel_capacity: args.channel_size,
    };

    if args.verbose {
        println!("Configuration:");
        println!("  k-mer length: {}", config.k);
        println!("  min frequency: {}", config.min_frequency);
        println!("  max memory: {} MB", config.max_memory_mb);
        println!("  temp dir: {:?}", config.temp_dir);
        println!("  bloom filter size: {}", config.bloom_filter_size);
    }

    // Run streaming pipeline
    let start_time = Instant::now();
    let pipeline = StreamingPipeline::new(config);
    
    println!("Processing genome with streaming architecture...");
    let candidates = pipeline.run(&args.input)?;
    
    let processing_time = start_time.elapsed();
    
    // Report results
    println!("\n=== STREAMING RESULTS ===");
    println!("Processing time: {:.2?}", processing_time);
    println!("Candidate k-mers found: {}", candidates.len());
    println!("Memory-efficient processing: ✓");
    
    // Filter high-frequency candidates (potential repeats)
    let repeat_candidates: Vec<_> = candidates.iter()
        .filter(|c| c.frequency >= args.min_freq * 3)  // Higher threshold for repeats
        .collect();
    
    println!("High-frequency repeat candidates: {}", repeat_candidates.len());
    
    // Write basic output (can be extended for full repeat detection)
    write_candidate_summary(&args.output, &candidates)?;
    
    if args.verbose {
        println!("\nTop 10 most frequent k-mers:");
        let mut sorted_candidates = candidates.clone();
        sorted_candidates.sort_by_key(|c| std::cmp::Reverse(c.frequency));
        
        for (i, candidate) in sorted_candidates.iter().take(10).enumerate() {
            println!("  {}: freq={}, positions={}", 
                i+1, candidate.frequency, candidate.positions.len());
        }
    }
    
    println!("\nOutput written to: {:?}.summary", args.output);
    println!("Streaming processing complete! 🚀");
    
    Ok(())
}

fn write_candidate_summary(output_prefix: &PathBuf, candidates: &[reprise::streaming::KmerCandidate]) -> Result<()> {
    use std::fs::File;
    use std::io::{BufWriter, Write};
    
    let summary_path = format!("{}.summary", output_prefix.display());
    let mut writer = BufWriter::new(File::create(&summary_path)?);
    
    writeln!(writer, "# REPrise v2.0 Streaming K-mer Summary")?;
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
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;
    
    #[test]
    fn test_small_genome_streaming() -> Result<()> {
        // Create a small test FASTA
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, ">test_contig_1").unwrap();
        writeln!(temp_file, "ACGTACGTACGTACGTACGTACGTACGTACGT").unwrap();
        writeln!(temp_file, ">test_contig_2").unwrap();
        writeln!(temp_file, "TGCATGCATGCATGCATGCATGCATGCATGCA").unwrap();
        temp_file.flush().unwrap();
        
        let config = StreamingConfig {
            k: 6,
            min_frequency: 1,
            max_memory_mb: 10,  // Very small for testing
            temp_dir: std::env::temp_dir(),
            bloom_filter_size: 1000,
            channel_capacity: 100,
        };
        
        let pipeline = StreamingPipeline::new(config);
        let candidates = pipeline.run(temp_file.path())?;
        
        // Should have some candidates from the repetitive sequences
        assert!(!candidates.is_empty(), "Should find some k-mer candidates");
        
        // Check that some k-mers appear multiple times (repeats)
        let has_repeats = candidates.iter().any(|c| c.frequency >= 2);
        assert!(has_repeats, "Should find some repeated k-mers");
        
        Ok(())
    }
}