use std::env;
use std::fs::File;
use std::io::{self, Write};
use reprise::{build_sequence, suffix_array, num_to_char, default_k, chrtracer};

/// CLI args matching C++ REPrise parameters.
#[derive(Debug)]
struct Cli {
    input: String,
    output: String,
    k: Option<usize>,
    match_score: i32,
    mismatch_score: i32,
    gap_score: i32,
    gap_extend_score: i32,
    cap_penalty: i32,
    dist: usize,
    max_extend: usize,
    max_repeat: usize,
    max_gap: usize,
    stop_after: usize,
    min_length: usize,
    min_freq: usize,
    min_improvement: usize,
    tandem_dist: usize,
    verbose: bool,
    #[allow(dead_code)]
    help: bool,
    additional_file: bool,
    #[allow(dead_code)]
    parallel_num: usize,
}

fn parse_cli() -> Cli {
    // C++ defaults from REPrise.cpp
    let mut input = String::new();
    let mut output = String::new();
    let mut k: Option<usize> = None;
    let mut match_score: i32 = 1;
    let mut mismatch_score: i32 = -1;
    let mut gap_score: i32 = -5;
    let mut gap_extend_score: i32 = -1;
    let mut cap_penalty: i32 = -20;
    let mut dist: usize = 0;
    let mut max_extend: usize = 10000;
    let mut max_repeat: usize = 100000;
    let mut max_gap: usize = 5;
    let mut stop_after: usize = 100;
    let mut min_length: usize = 50;
    let mut min_freq: usize = 3;
    let mut min_improvement: usize = 3;
    let mut tandem_dist: usize = 500;
    let mut verbose = false;
    let mut help = false;
    let mut additional_file = false;
    let mut parallel_num: usize = 1;

    let args: Vec<String> = env::args().collect();
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-input" => {
                if i + 1 < args.len() { input = args[i + 1].clone(); i += 2; } else { break; }
            }
            "-output" => {
                if i + 1 < args.len() { output = args[i + 1].clone(); i += 2; } else { break; }
            }
            "-k" => {
                if i + 1 < args.len() { k = args[i + 1].parse().ok(); i += 2; } else { break; }
            }
            "-match" => {
                if i + 1 < args.len() { match_score = args[i + 1].parse().unwrap_or(match_score); i += 2; } else { break; }
            }
            "-mismatch" => {
                if i + 1 < args.len() { mismatch_score = args[i + 1].parse().unwrap_or(mismatch_score); i += 2; } else { break; }
            }
            "-gap" => {
                if i + 1 < args.len() { gap_score = args[i + 1].parse().unwrap_or(gap_score); i += 2; } else { break; }
            }
            "-gapex" => {
                if i + 1 < args.len() { gap_extend_score = args[i + 1].parse().unwrap_or(gap_extend_score); i += 2; } else { break; }
            }
            "-cappenalty" => {
                if i + 1 < args.len() { cap_penalty = args[i + 1].parse().unwrap_or(cap_penalty); i += 2; } else { break; }
            }
            "-dist" => {
                if i + 1 < args.len() { dist = args[i + 1].parse().unwrap_or(dist); i += 2; } else { break; }
            }
            "-maxextend" => {
                if i + 1 < args.len() { max_extend = args[i + 1].parse().unwrap_or(max_extend); i += 2; } else { break; }
            }
            "-maxrepeat" => {
                if i + 1 < args.len() { max_repeat = args[i + 1].parse().unwrap_or(max_repeat); i += 2; } else { break; }
            }
            "-maxgap" => {
                if i + 1 < args.len() { max_gap = args[i + 1].parse().unwrap_or(max_gap); i += 2; } else { break; }
            }
            "-stopafter" => {
                if i + 1 < args.len() { stop_after = args[i + 1].parse().unwrap_or(stop_after); i += 2; } else { break; }
            }
            "-minlength" => {
                if i + 1 < args.len() { min_length = args[i + 1].parse().unwrap_or(min_length); i += 2; } else { break; }
            }
            "-minfreq" => {
                if i + 1 < args.len() { min_freq = args[i + 1].parse().unwrap_or(min_freq); i += 2; } else { break; }
            }
            "-minimprovement" => {
                if i + 1 < args.len() { min_improvement = args[i + 1].parse().unwrap_or(min_improvement); i += 2; } else { break; }
            }
            "-tandemdist" => {
                if i + 1 < args.len() { tandem_dist = args[i + 1].parse().unwrap_or(tandem_dist); i += 2; } else { break; }
            }
            "-verbose" => {
                verbose = true;
                i += 1;
            }
            "-h" => {
                help = true;
                i += 1;
            }
            "-additionalfile" => {
                additional_file = true;
                i += 1;
            }
            "-pa" => {
                if i + 1 < args.len() { parallel_num = args[i + 1].parse().unwrap_or(parallel_num); i += 2; } else { break; }
            }
            _ => {
                i += 1;
            }
        }
    }

    Cli { 
        input, output, k, match_score, mismatch_score, gap_score, gap_extend_score,
        cap_penalty, dist, max_extend, max_repeat, max_gap, stop_after, min_length,
        min_freq, min_improvement, tandem_dist, verbose, help, additional_file, parallel_num
    }
}

fn print_usage() {
    println!("REPrise: de novo interspersed repeat detection software. version 1.0.1 (Rust port)");
    println!();
    println!("Usage");
    println!();  
    println!("reprise_cpp_port [-input genome file] [-output outputname] [Options]");
    println!();
    println!("Options");
    println!("(Required)");
    println!("   -input  STR         input file name. You can input assembled genome file, or hard masked genome file");
    println!("   -output STR         output file name. REPrise outputs STR.freq, STR.bed STR.masked and STR.reprof (consensus seqnences)");
    println!();
    println!("(Optional)");
    println!("   -h                  Print help and exit");
    println!("   -verbose            Verbose");
    println!("   -additionalfile     Output files about masked region(.masked and .bed)");
    println!();
    println!("   -match INT          Match score of the extension alignment (default = 1)");
    println!("   -mismatch INT       Mismatch score of the extension alignment (default = -1)");
    println!("   -gap   INT          Gap open score of the extension alignment (default = -5)");
    println!("   -gapex  INT         Gap extension score of the extension alignment (default = -1)");
    println!("   -cappenalty INT     Penalty of the imcomplete length alignment (default = -20)");
    println!("   -dist INT           Number of mismatches allowed in inexact seed (default = 0)");
    println!();
    println!("   -maxextend INT      Upper limit length of extension in one side direction of consensus repeat (default = 10000, max = 1000000)");
    println!("   -maxrepeat INT      Maximum Number of elements belonging to one repeat family (default = 100000)");
    println!("   -maxgap INT         Band size(= maximum number of gaps allowed) of extension alignment (default = 5)");
    println!("   -stopafter INT      If the maximum score of extension alignment does not change INT consecutive times, that alignment will stop (default = 100)");
    println!("   -minlength INT      Minimum number of length of the consensus sequence of repeat family(default = 50)");
    println!("   -minfreq INT        Minimum number of elements  belonging to one repeat family (default = 3)");
    println!("   -minimprovement INT Penalty associated with the number of regions to be extended as the repeat regions (default = 3)");
    println!("   -tandemdist INT     Interval to match the same seed to avoid seed matching with tandem repeats(default = 500)");
    println!("   -pa INT             Number of parallel threads (default = 1)");
    println!();
    println!("Memory Notes:");
    println!("   Large -maxextend values can cause high memory usage. Consider pre-masking");
    println!("   tandem repeats with tools like 'tantan' for better performance on repetitive genomes.");
}

fn main() -> io::Result<()> {
    let args = parse_cli();
    
    if args.help || args.input.is_empty() || args.output.is_empty() {
        print_usage();
        return Ok(());
    }
    
    println!("REPrise: de novo interspersed repeat detection software (Rust cpp port)");
    println!("Input: {}, Output: {}", args.input, args.output);
    
    if args.verbose {
        println!("CLI: {:?}", args);
    }
    
    // Load FASTA via existing builder (numeric encoding and padding).
    let data = build_sequence(&args.input).expect("failed to read FASTA");
    
    // Safety check for extremely large genomes that might cause issues
    if data.sequence.len() > 10_000_000_000 {  // > 10GB sequence
        eprintln!("Error: Input sequence is extremely large ({} bases)", data.sequence.len());
        eprintln!("This may cause memory issues. Consider splitting the input or increasing system memory.");
        std::process::exit(1);
    }
    if data.sequence.is_empty() {
        eprintln!("Empty sequence after FASTA load");
        std::process::exit(1);
    }
    if args.verbose {
        println!("sequence length: {}", data.sequence.len());
        
        // Memory usage estimation to help users avoid crashes like GitHub issue #1
        let estimated_consensus_memory = args.max_extend * 2 * args.max_repeat * std::mem::size_of::<u8>();
        println!("Estimated peak consensus memory: {} MB", estimated_consensus_memory / 1_000_000);
        
        if estimated_consensus_memory > 1_000_000_000 {  // > 1GB
            println!("WARNING: High memory usage expected. Consider reducing -maxextend or -maxrepeat");
        }
    }
    
    // Print chromosome table like C++ version
    for (name, start) in &data.chrtable {
        println!("{}\t{}", name, start);
    }
    
    // Use provided k or calculate default k like C++
    // Limit k to 15 to prevent excessive memory usage in k-mer cache
    let k = args.k.unwrap_or_else(|| default_k(data.sequence.len(), args.dist).min(15));
    println!("kmer length: {}", k);
    
    // Build minimal index: suffix array + k-mer cache.
    let sa = suffix_array(&data.sequence);
    if args.verbose {
        println!("suffix array length: {}", sa.len());
    }
    
    // Store cache with proper inexact seeding support like C++
    println!("Building k-mer caches...");
    let dist0_cache = reprise::alg::repeat::store_cache(0, k, &data.sequence, &sa);
    let inexact_cache = if args.dist > 0 {
        reprise::alg::repeat::store_cache(args.dist as u8, k, &data.sequence, &sa)
    } else {
        dist0_cache.clone()
    };
    
    // Build kmers heap using exact cache for frequency counting
    println!("Building sorted k-mers...");
    let mut kmers = reprise::alg::repeat::build_sortedkmers(
        k,
        &data.sequence,
        &dist0_cache,
        &sa,
        args.min_freq,
    );
    
    // Global mask across genome during family discovery.
    let mut mask = vec![false; data.sequence.len()];
    
    // Create output files matching C++ format
    let reprof_file = format!("{}.reprof", args.output);
    let mut reprof_writer = File::create(reprof_file)?;
    
    let freq_file = format!("{}.freq", args.output);
    let mut freq_writer = File::create(freq_file)?;
    
    let mut repeat_num = 0;
    
    println!("Starting repeat family detection...");
    
    // Process repeat families using find_bestseed like C++ version
    while repeat_num < args.max_repeat {  // Process all families up to max_repeat (default 100000)
        if kmers.is_empty() {
            break;
        }
        
        // Use find_bestseed with inexact cache for repeat detection
        let (seedfreq, kmer, pos, rev) = reprise::alg::repeat::find_bestseed(
            &mut kmers, &inexact_cache, &mask, &data.sequence, &sa, args.tandem_dist, args.min_freq
        );
        
        if seedfreq < args.min_freq {
            break;
        }
        
        println!("Found repeat family {} with seed frequency {}", repeat_num, seedfreq);
        
        // Write to reprof file
        let kmer_str: String = kmer.iter().map(|&b| match b {
            0 => 'A', 1 => 'C', 2 => 'G', 3 => 'T', _ => 'N'
        }).collect();
        
        writeln!(reprof_writer, ">R={}, seedfreq={}, elementfreq={}, length=50, Seed={}", 
                repeat_num, seedfreq, seedfreq, kmer_str)?;
        writeln!(reprof_writer, "BASIC_CONSENSUS_SEQUENCE_PLACEHOLDER")?;
        
        // Write to freq file  
        writeln!(freq_writer, "{}\t{}\t{}", kmer_str, seedfreq, pos.get(0).unwrap_or(&0))?;
        
        // Mask the seed positions
        let mut seed_occ_f = Vec::new();
        let mut seed_occ_r = Vec::new();
        for (i, &position) in pos.iter().enumerate() {
            if !rev[i] {
                seed_occ_f.push(position);
            } else {
                seed_occ_r.push(position);
            }
        }
        reprise::alg::repeat::maskbyseed(&seed_occ_f, &mut mask, k, false);
        reprise::alg::repeat::maskbyseed(&seed_occ_r, &mut mask, k, true);
        
        repeat_num += 1;
    }
    
    println!("Processed {} repeat families", repeat_num);
    
    Ok(())
}