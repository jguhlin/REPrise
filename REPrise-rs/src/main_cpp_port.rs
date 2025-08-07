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
            "-verbose" => { verbose = true; i += 1; }
            "-h" => { help = true; i += 1; }
            "-additionalfile" => { additional_file = true; i += 1; }
            "-pa" => {
                if i + 1 < args.len() { parallel_num = args[i + 1].parse().unwrap_or(parallel_num); i += 2; } else { break; }
            }
            _ => { i += 1; }
        }
    }

    if help {
        print_usage();
        std::process::exit(0);
    }

    if input.is_empty() || output.is_empty() {
        eprintln!("Error: -input and -output are required");
        print_usage();
        std::process::exit(1);
    }

    // Memory safety validation - prevent heap-buffer-overflow like GitHub issue #1
    if max_extend > 1_000_000 {
        eprintln!("Error: -maxextend {} is too large (max: 1,000,000)", max_extend);
        eprintln!("Large values can cause memory issues. Consider using tandem repeat masking first.");
        std::process::exit(1);
    }

    if max_repeat > 10_000_000 {
        eprintln!("Error: -maxrepeat {} is too large (max: 10,000,000)", max_repeat);
        std::process::exit(1);
    }

    // Warn about potentially problematic parameter combinations
    if max_extend > 50_000 {
        eprintln!("Warning: -maxextend {} is very large and may consume significant memory", max_extend);
        eprintln!("Consider pre-masking tandem repeats with tools like tantan for better performance");
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
    println!("REPrise [-input genome file] [-output outputname] [Options]");
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

#[derive(Debug)]
#[allow(dead_code)]
struct RepeatHit {
    contig_id: String,
    start: usize,
    end: usize,
    length: usize,
    orientation: char, // '+' or '-'
    score: i32,
}

fn main() -> io::Result<()> {
    let args = parse_cli();
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
    let k = args.k.unwrap_or_else(|| default_k(data.sequence.len(), args.dist).min(15));  // Limit for memory safety
    println!("kmer length: {}", k);

    // Build minimal index: suffix array + k-mer cache.
    let sa = suffix_array(&data.sequence);
    if args.verbose {
        println!("suffix array length: {}", sa.len());
    }
    
    // Store cache with proper inexact seeding support like C++
    let dist0_cache = reprise::alg::repeat::store_cache(0, k, &data.sequence, &sa);
    let inexact_cache = if args.dist > 0 {
        reprise::alg::repeat::store_cache(args.dist as u8, k, &data.sequence, &sa)
    } else {
        dist0_cache.clone()
    };

    // Build kmers heap using exact cache for frequency counting
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
    
    // Create BED and masked files if requested
    let mut bed_writer = if args.additional_file {
        let bed_file = format!("{}.bed", args.output);
        Some(File::create(bed_file)?)
    } else {
        None
    };
    
    let mut masked_writer = if args.additional_file {
        let masked_file = format!("{}.masked", args.output);
        Some(File::create(masked_file)?)
    } else {
        None
    };
    
    let mut repeat_num = 0;
    
    // Process repeat families using find_bestseed like C++ version with full extension alignment
    while repeat_num < args.max_repeat {
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

        // Immediately mask the seed positions to prevent overlapping k-mers
        // This is critical for matching C++ behavior
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

        // Create consensus sequence with extension alignment - safe allocation
        let consensus_size = 2 * args.max_extend + k;
        
        // Memory safety check - prevent excessive allocations that caused C++ crashes
        if consensus_size > 10_000_000 {
            eprintln!("Error: Consensus sequence size {} too large for repeat family {}", 
                     consensus_size, repeat_num);
            eprintln!("This can happen with very large -maxextend values on repetitive genomes");
            break; // Skip this family and continue processing
        }
        
        // Safe allocation - Rust prevents buffer overflows but we check size limits
        let mut consensus = vec![0u8; consensus_size];
        
        // Initialize consensus with the seed sequence at MAXEXTEND position
        for (i, &base) in kmer.iter().enumerate() {
            if args.max_extend + i < consensus.len() {
                consensus[args.max_extend + i] = base;
            }
        }
        
        // Perform extension alignment in both directions
        let mut repeatstart = vec![0i32; seedfreq];
        let mut repeatend = vec![0i32; seedfreq];
        let mut seed_ext = vec![-1i32; seedfreq];
        
        // Right extension
        let right_ext = reprise::alg::repeat::extend(
            true, // is_right
            seedfreq,
            &pos,
            &rev,
            &data.sequence,
            &mut consensus,
            &mut seed_ext,
            args.match_score,
            args.mismatch_score,
            args.gap_score,
            args.gap_extend_score,
            args.cap_penalty,
            args.max_extend,
            args.max_gap,
            args.stop_after,
            k,
            args.min_improvement as i32,
        );
        
        // Update seed_ext for right extension
        for i in 0..seedfreq {
            repeatend[i] = seed_ext[i];
        }
        
        // Reset seed_ext for left extension
        seed_ext.fill(-1);
        
        // Left extension  
        let left_ext = reprise::alg::repeat::extend(
            false, // is_right
            seedfreq,
            &pos,
            &rev,
            &data.sequence,
            &mut consensus,
            &mut seed_ext,
            args.match_score,
            args.mismatch_score,
            args.gap_score,
            args.gap_extend_score,
            args.cap_penalty,
            args.max_extend,
            args.max_gap,
            args.stop_after,
            k,
            args.min_improvement as i32,
        );
        
        // Update seed_ext for left extension
        for i in 0..seedfreq {
            repeatstart[i] = -seed_ext[i] - 1;
        }
        
        // Calculate consensus boundaries  
        let consensusstart = args.max_extend.saturating_sub(left_ext.max(0) as usize);
        let consensusend = args.max_extend + k - 1 + right_ext.max(0) as usize;
        let consensus_length = consensusend.saturating_sub(consensusstart) + 1;
        
        // Only create repeat family if it meets minimum length requirement
        if consensus_length >= args.min_length {
            let mut element_count = 0;
            
            // Perform masking alignment for each occurrence
            for i in 0..seedfreq {
                let (elementstart, elementend) = reprise::alg::repeat::masking_align(
                    i,
                    consensusstart,
                    consensusend,
                    &consensus,
                    &data.sequence,
                    pos[i],
                    k,
                    rev[i],
                    args.match_score,
                    args.mismatch_score,
                    args.gap_score,
                    args.gap_extend_score,
                    args.max_extend,
                    args.stop_after,
                );
                
                // Check if element meets minimum length
                if (elementend - elementstart + 1).abs() >= args.min_length as isize {
                    element_count += 1;
                    
                    // Write to BED file if requested
                    if let Some(ref mut bed) = bed_writer {
                        let element_pos = pos[i];
                        let (chr_name, chr_offset) = chrtracer(element_pos, &data.chrtable);
                        let bed_start = element_pos + elementstart.max(0) as usize - chr_offset;
                        let bed_end = element_pos + elementend.max(0) as usize - chr_offset;
                        let element_length = (elementend - elementstart + 1).abs();
                        let strand = if rev[i] { "-" } else { "+" };
                        writeln!(bed, "{}\t{}\t{}\tR={}\t{}\t{}", 
                            chr_name, bed_start, bed_end, repeat_num, element_length, strand)?;
                    }
                    
                    // Mask the element
                    reprise::alg::repeat::maskbyrepeat_element(i, elementstart.max(0) as usize, elementend.max(0) as usize, &mut mask, &pos);
                }
            }
            
            // Only output if we have valid elements
            if element_count > 0 {
                // Create consensus sequence output
                write!(reprof_writer, ">R={}, seedfreq={}, elementfreq={}, length={}, Seed=", 
                       repeat_num, seedfreq, element_count, consensus_length)?;
                for &base in &kmer {
                    write!(reprof_writer, "{}", num_to_char(base))?;
                }
                writeln!(reprof_writer)?;
                
                // Write the extended consensus sequence
                for i in consensusstart..=consensusend {
                    write!(reprof_writer, "{}", num_to_char(consensus[i]))?;
                    if (i - consensusstart + 1) % 80 == 0 {
                        writeln!(reprof_writer)?;
                    }
                }
                if (consensusend - consensusstart + 1) % 80 != 0 {
                    writeln!(reprof_writer)?;
                }
                
                // Write to freq file (format: repeat_id \t seedfreq \t element_count)
                writeln!(freq_writer, "R={}\t{}\t{}", repeat_num, seedfreq, element_count)?;
                
                repeat_num += 1;
            }
            // Note: Seeds are already masked immediately after find_bestseed
        }
        // Note: Seeds are already masked immediately after find_bestseed
    }

    // Write masked sequence file if requested
    if let Some(ref mut masked) = masked_writer {
        for (chr_name, chr_start) in &data.chrtable {
            if chr_name == "unknown" || chr_name == "padding" {
                continue; // Skip special chromosome entries
            }
            
            writeln!(masked, ">{}", chr_name)?;
            
            // Find the end of this chromosome
            let chr_end = if let Some((_name, next_start)) = data.chrtable.iter()
                .find(|(n, _)| n != chr_name && n != "unknown" && n != "padding") {
                *next_start
            } else {
                data.sequence.len()
            };
            
            // Write sequence with masked regions in lowercase
            for i in *chr_start..chr_end.min(data.sequence.len()) {
                if data.sequence[i] > 3 {
                    continue; // Skip invalid bases
                }
                let base_char = num_to_char(data.sequence[i]);
                let output_char = if mask[i] {
                    base_char.to_ascii_lowercase()
                } else {
                    base_char
                };
                write!(masked, "{}", output_char)?;
                
                // Line wrap at 80 characters
                if (i - *chr_start + 1) % 80 == 0 {
                    writeln!(masked)?;
                }
            }
            
            // Add final newline if needed
            let chr_len = chr_end.min(data.sequence.len()).saturating_sub(*chr_start);
            if chr_len % 80 != 0 {
                writeln!(masked)?;
            }
        }
    }

    println!("Processed {} repeat families", repeat_num);
    Ok(())
}

