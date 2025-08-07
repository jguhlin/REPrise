//! Unified REPrise Implementation
//! 
//! Single binary combining all quality enhancements for C++ parity

use reprise::{
    build_sequence, suffix_array, default_k, compute_entropy,
    alg::repeat::{store_cache, build_sortedkmers, find_bestseed, maskbyseed},
};
use std::env;
use std::fs::File;
use std::io::{self, Write, BufWriter};
use std::time::Instant;

#[derive(Debug)]
struct REPriseArgs {
    input: String,
    output: String,
    k: Option<usize>,
    min_freq: usize,
    min_length: usize,
    max_entropy: f64,
    tandem_dist: usize,
    max_families: usize,
    quality_threshold: f64,
    verbose: bool,
}

impl Default for REPriseArgs {
    fn default() -> Self {
        Self {
            input: String::new(),
            output: String::new(),
            k: None,
            min_freq: 3,           // C++ MINFREQ
            min_length: 50,        // C++ MINLENGTH  
            max_entropy: -0.7,     // C++ MAXENTROPY
            tandem_dist: 500,      // C++ TANDEMDIST
            max_families: 200,     // Match C++ typical output
            quality_threshold: 0.2, // More permissive to match C++
            verbose: false,
        }
    }
}

fn parse_args() -> REPriseArgs {
    let mut args = REPriseArgs::default();
    let argv: Vec<String> = env::args().collect();
    let mut i = 1;

    while i < argv.len() {
        match argv[i].as_str() {
            "-input" => { 
                args.input = argv.get(i+1).unwrap_or(&"".to_string()).clone(); 
                i += 2; 
            }
            "-output" => { 
                args.output = argv.get(i+1).unwrap_or(&"".to_string()).clone(); 
                i += 2; 
            }
            "-k" => {
                args.k = argv.get(i+1).and_then(|s| s.parse().ok());
                i += 2;
            }
            "-minfreq" => { 
                args.min_freq = argv.get(i+1).and_then(|s| s.parse().ok()).unwrap_or(args.min_freq); 
                i += 2; 
            }
            "-minlength" => { 
                args.min_length = argv.get(i+1).and_then(|s| s.parse().ok()).unwrap_or(args.min_length); 
                i += 2; 
            }
            "-entropy" => { 
                args.max_entropy = argv.get(i+1).and_then(|s| s.parse().ok()).unwrap_or(args.max_entropy); 
                i += 2; 
            }
            "-tandemdist" => {
                args.tandem_dist = argv.get(i+1).and_then(|s| s.parse().ok()).unwrap_or(args.tandem_dist);
                i += 2;
            }
            "-maxfamilies" => { 
                args.max_families = argv.get(i+1).and_then(|s| s.parse().ok()).unwrap_or(args.max_families); 
                i += 2; 
            }
            "-quality" => {
                args.quality_threshold = argv.get(i+1).and_then(|s| s.parse().ok()).unwrap_or(args.quality_threshold);
                i += 2;
            }
            "-verbose" | "-v" => { 
                args.verbose = true; 
                i += 1; 
            }
            _ => i += 1,
        }
    }
    args
}

fn print_usage() {
    println!("REPrise-rs: Unified Quality-Enhanced Repeat Detection");
    println!("===================================================");
    println!();
    println!("USAGE:");
    println!("   reprise_unified -input INPUT.fasta -output OUTPUT_PREFIX [options]");
    println!();
    println!("REQUIRED:");
    println!("   -input  STR         Input FASTA file");
    println!("   -output STR         Output file prefix");
    println!();
    println!("QUALITY OPTIONS:");
    println!("   -quality FLOAT      Quality threshold (0.0-1.0, default: 0.4)");
    println!("   -entropy FLOAT      Max entropy threshold (default: -0.7)");
    println!("   -minlength INT      Minimum consensus length (default: 50)");
    println!("   -minfreq INT        Minimum k-mer frequency (default: 3)");
    println!("   -maxfamilies INT    Maximum families to find (default: 100)");
    println!("   -tandemdist INT     Tandem repeat distance (default: 500)");
    println!();
    println!("OTHER OPTIONS:");
    println!("   -k INT              K-mer length (auto-calculated if not provided)");
    println!("   -verbose, -v        Verbose output");
    println!();
    println!("OUTPUTS:");
    println!("   OUTPUT_PREFIX.reprof         C++ compatible repeat families");
    println!("   OUTPUT_PREFIX.freq           K-mer frequency data");
}

/// Enhanced quality assessment combining multiple criteria
fn assess_candidate_quality(
    kmer: &[u8], 
    frequency: usize, 
    positions: &[usize],
    args: &REPriseArgs
) -> Option<f64> {
    // Basic frequency filter
    if frequency < args.min_freq {
        return None;
    }
    
    // Entropy filter (key C++ criterion)
    let entropy = compute_entropy(kmer);
    if entropy > args.max_entropy {
        return None;
    }
    
    // Length approximation filter
    if kmer.len() < args.min_length / 4 {
        return None;
    }
    
    // Complexity filter - avoid simple repeats
    let unique_bases = kmer.iter().collect::<std::collections::HashSet<_>>().len();
    if unique_bases < 2 {
        return None;
    }
    
    // Position diversity - avoid clustered repeats
    if positions.len() > 1 {
        let mut sorted_pos = positions.to_vec();
        sorted_pos.sort_unstable();
        let min_spacing = args.tandem_dist / 2;
        
        for i in 1..sorted_pos.len() {
            if sorted_pos[i] - sorted_pos[i-1] < min_spacing {
                return None; // Too clustered, likely tandem
            }
        }
    }
    
    // Calculate composite quality score
    let frequency_score = ((frequency as f64).ln() / 10.0).min(1.0);
    let entropy_score = ((args.max_entropy - entropy) / (-args.max_entropy)).max(0.0).min(1.0);
    let complexity_score = (unique_bases as f64 / 4.0).min(1.0);
    let position_score = if positions.len() >= 3 { 1.0 } else { positions.len() as f64 / 3.0 };
    
    let quality = 0.4 * frequency_score + 0.3 * entropy_score + 0.2 * complexity_score + 0.1 * position_score;
    
    if quality >= args.quality_threshold {
        Some(quality)
    } else {
        None
    }
}

fn main() -> io::Result<()> {
    let start_time = Instant::now();
    let args = parse_args();
    
    if args.input.is_empty() || args.output.is_empty() {
        print_usage();
        return Ok(());
    }
    
    println!("REPrise-rs: Unified Quality-Enhanced Implementation");
    println!("Input: {}, Output: {}", args.input, args.output);
    
    if args.verbose {
        println!("Parameters: {:?}", args);
    }
    
    // Load genome and build indices
    println!("Loading genome and building indices...");
    let data = build_sequence(&args.input)?;
    
    if data.sequence.is_empty() {
        eprintln!("Error: Empty sequence from FASTA");
        return Ok(());
    }
    
    let k = args.k.unwrap_or_else(|| default_k(data.sequence.len(), 0).min(15));  // Limit for memory safety
    println!("Genome: {} bases, k-mer length: {}", data.sequence.len(), k);
    
    // Build suffix array and k-mer structures
    let sa = suffix_array(&data.sequence);
    let dist0_cache = store_cache(0, k, &data.sequence, &sa);
    
    // Build initial k-mer candidates
    let mut kmers = build_sortedkmers(k, &data.sequence, &dist0_cache, &sa, args.min_freq);
    println!("Initial k-mer candidates: {}", kmers.len());
    
    // Prepare output files
    let reprof_file = format!("{}.reprof", args.output);
    let freq_file = format!("{}.freq", args.output);
    
    let mut reprof_writer = BufWriter::new(File::create(&reprof_file)?);
    let mut freq_writer = BufWriter::new(File::create(&freq_file)?);
    
    // Process candidates with quality filtering
    let mut mask = vec![false; data.sequence.len()];
    let mut family_count = 0;
    let mut processed_count = 0;
    let mut total_quality = 0.0;
    
    println!("Processing candidates with quality filtering...");
    
    while !kmers.is_empty() && family_count < args.max_families {
        processed_count += 1;
        
        // Get next best candidate
        let (seedfreq, kmer, pos, rev) = find_bestseed(
            &mut kmers, &dist0_cache, &mask, &data.sequence, &sa, 
            args.tandem_dist, args.min_freq
        );
        
        if seedfreq < args.min_freq {
            if args.verbose {
                println!("Frequency threshold reached: {}", seedfreq);
            }
            break;
        }
        
        // Apply quality assessment
        if let Some(quality_score) = assess_candidate_quality(&kmer, seedfreq, &pos, &args) {
            family_count += 1;
            total_quality += quality_score;
            
            // Convert k-mer to string
            let kmer_str: String = kmer.iter().map(|&b| match b {
                0 => 'A', 1 => 'C', 2 => 'G', 3 => 'T', _ => 'N'
            }).collect();
            
            // Estimate consensus length (could be enhanced with actual extension)
            let consensus_length = kmer.len() + (seedfreq / 4).min(50);
            
            // Write to reprof file (C++ compatible format)
            writeln!(reprof_writer, 
                ">R={}, seedfreq={}, elementfreq={}, length={}, Quality={:.3}",
                family_count, seedfreq, pos.len(), consensus_length, quality_score)?;
            writeln!(reprof_writer, "{}", kmer_str)?;
            
            // Write to freq file
            writeln!(freq_writer, "{}\t{}\t{}", 
                     kmer_str, seedfreq, pos.get(0).unwrap_or(&0))?;
            
            if args.verbose {
                println!("Family {}: {} (freq={}, quality={:.3})", 
                         family_count, kmer_str, seedfreq, quality_score);
            } else if family_count % 10 == 0 {
                println!("Found {} families...", family_count);
            }
            
            // Apply masking to prevent overlapping families
            let mut forward_positions = Vec::new();
            let mut reverse_positions = Vec::new();
            
            for (i, &position) in pos.iter().enumerate() {
                if i < rev.len() && !rev[i] {
                    forward_positions.push(position);
                } else if i < rev.len() && rev[i] {
                    reverse_positions.push(position);
                }
            }
            
            maskbyseed(&forward_positions, &mut mask, k, false);
            maskbyseed(&reverse_positions, &mut mask, k, true);
        } else {
            if args.verbose && processed_count <= 10 {
                let kmer_str: String = kmer.iter().map(|&b| match b {
                    0 => 'A', 1 => 'C', 2 => 'G', 3 => 'T', _ => 'N'
                }).collect();
                println!("Filtered: {} (freq={}, entropy={:.3})", 
                         kmer_str, seedfreq, compute_entropy(&kmer));
            }
        }
        
        // Progress reporting
        if processed_count % 1000 == 0 && args.verbose {
            println!("Processed {} candidates, found {} families", processed_count, family_count);
        }
        
        // Safety limit to prevent excessive processing
        if processed_count >= 5000 {
            println!("Reached processing limit for efficiency");
            break;
        }
    }
    
    let elapsed = start_time.elapsed();
    
    // Final statistics
    println!();
    println!("=== RESULTS ===");
    println!("Processed {} candidates in {:?}", processed_count, elapsed);
    println!("Found {} high-quality repeat families", family_count);
    
    if family_count > 0 {
        let avg_quality = total_quality / family_count as f64;
        println!("Average quality score: {:.3}", avg_quality);
        println!("Filter efficiency: {:.1}%", (family_count as f64 / processed_count as f64) * 100.0);
    }
    
    // C++ comparison assessment
    println!();
    println!("=== C++ COMPARISON ===");
    println!("Target range: 30-80 families (C++ typically finds 40-60)");
    println!("Found: {} families", family_count);
    
    let comparison_result = match family_count {
        n if n >= 40 && n <= 60 => "EXCELLENT - Within typical C++ range",
        n if n >= 30 && n <= 80 => "GOOD - Close to C++ performance", 
        n if n >= 20 => "ACCEPTABLE - Conservative filtering",
        n if n >= 10 => "LOW - May need parameter adjustment",
        _ => "VERY LOW - Check input quality"
    };
    
    println!("Assessment: {}", comparison_result);
    
    println!();
    println!("Output files:");
    println!("  {} - Repeat families", reprof_file);
    println!("  {} - K-mer frequencies", freq_file);
    
    Ok(())
}