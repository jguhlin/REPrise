use std::env;
use std::fs::File;
use std::io::{self, Write};
use reprise::{build_sequence, suffix_array, reverse_complement, num_to_char, default_k};

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
    help: bool,
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

    Cli { 
        input, output, k, match_score, mismatch_score, gap_score, gap_extend_score,
        cap_penalty, dist, max_extend, max_repeat, max_gap, stop_after, min_length,
        min_freq, min_improvement, tandem_dist, verbose, help
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
    println!();
    println!("   -match INT          Match score of the extension alignment (default = 1)");
    println!("   -mismatch INT       Mismatch score of the extension alignment (default = -1)");
    println!("   -gap   INT          Gap open score of the extension alignment (default = -5)");
    println!("   -gapex  INT         Gap extension score of the extension alignment (default = -1)");
    println!("   -cappenalty INT     Penalty of the imcomplete length alignment (default = -20)");
    println!("   -dist INT           Number of mismatches allowed in inexact seed (default = 0)");
    println!();
    println!("   -maxextend INT      Upper limit length of extension in one side direction of consensus repeat (default = 10000)");
    println!("   -maxrepeat INT      Maximum Number of elements belonging to one repeat family (default = 100000)");
    println!("   -maxgap INT         Band size(= maximum number of gaps allowed) of extension alignment (default = 5)");
    println!("   -stopafter INT      If the maximum score of extension alignment does not change INT consecutive times, that alignment will stop (default = 100)");
    println!("   -minlength INT      Minimum number of length of the consensus sequence of repeat family(default = 50)");
    println!("   -minfreq INT        Minimum number of elements  belonging to one repeat family (default = 3)");
    println!("   -minimprovement INT Penalty associated with the number of regions to be extended as the repeat regions (default = 3)");
    println!("   -tandemdist INT     Interval to match the same seed to avoid seed matching with tandem repeats(default = 500)");
    println!();
}

#[derive(Debug)]
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
    if data.sequence.is_empty() {
        eprintln!("Empty sequence after FASTA load");
        std::process::exit(1);
    }
    if args.verbose {
        println!("sequence length: {}", data.sequence.len());
    }

    // Print chromosome table like C++ version
    for (name, start) in &data.chrtable {
        println!("{}\t{}", name, start);
    }

    // Use provided k or calculate default k like C++
    let k = args.k.unwrap_or_else(|| default_k(data.sequence.len(), args.dist));
    println!("kmer length: {}", k);

    // Build minimal index: suffix array + k-mer cache.
    let sa = suffix_array(&data.sequence);
    if args.verbose {
        println!("suffix array length: {}", sa.len());
    }
    
    // Store cache with edit distance 0 for exact matching
    let cache = reprise::alg::repeat::store_cache(0, k, &data.sequence, &sa);

    // Build kmers heap and iterate families as a proxy for repeats.
    let mut kmers = reprise::alg::repeat::build_sortedkmers(
        k,
        &data.sequence,
        &cache,
        &sa,
        args.min_freq,
    );

    // Global mask across genome during family discovery.
    let mut mask = vec![false; data.sequence.len()];

    // For now, create a simple output file matching C++ .reprof format
    let reprof_file = format!("{}.reprof", args.output);
    let mut reprof_writer = File::create(reprof_file)?;
    
    let mut repeat_num = 0;
    
    // Process repeat families (simplified version for now)
    while let Some((freq, kmer)) = kmers.pop() {
        if freq < args.min_freq { break; }
        if repeat_num >= args.max_repeat { break; }

        // Forward occurrences
        let mut occ_f = reprise::alg::repeat::findkmer(&kmer, &cache, &data.sequence, &sa);
        reprise::alg::repeat::removetandem(&mut occ_f, args.tandem_dist);
        reprise::alg::repeat::removemasked(&mut occ_f, &mask, k, false);

        // Reverse complement occurrences
        let rc_kmer = reverse_complement(&kmer);
        let mut occ_r = reprise::alg::repeat::findkmer(&rc_kmer, &cache, &data.sequence, &sa);
        reprise::alg::repeat::removetandem(&mut occ_r, args.tandem_dist);
        reprise::alg::repeat::removemasked(&mut occ_r, &mask, k, true);

        let total_freq = occ_f.len() + occ_r.len();
        if total_freq < args.min_freq { continue; }

        // Mask the used positions
        reprise::alg::repeat::maskbyseed(&occ_f, &mut mask, k, false);
        reprise::alg::repeat::maskbyseed(&occ_r, &mut mask, k, true);

        // Create a simple consensus sequence (just the seed for now)
        write!(reprof_writer, ">R={}, seedfreq={}, elementfreq={}, length={}, Seed=", 
               repeat_num, freq, total_freq, k)?;
        for &base in &kmer {
            write!(reprof_writer, "{}", num_to_char(base))?;
        }
        writeln!(reprof_writer)?;
        
        // Write the consensus (just the seed sequence for now)
        for &base in &kmer {
            write!(reprof_writer, "{}", num_to_char(base))?;
        }
        writeln!(reprof_writer)?;

        repeat_num += 1;
    }

    println!("Processed {} repeat families", repeat_num);
    Ok(())
}

