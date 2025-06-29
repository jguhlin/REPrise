use clap::Parser;
use reprise::build_sequence;
use reprise::suffix_array;

/// Normalize single-hyphen long options (e.g. `-input`) to clap style `--input`.
fn normalize_args() -> Vec<String> {
    let mut args = std::env::args().collect::<Vec<_>>();
    for arg in &mut args[1..] {
        if arg.starts_with('-') && !arg.starts_with("--") && arg.len() > 2 {
            *arg = format!("--{}", &arg[1..]);
        }
    }
    args
}

/// Command line arguments for REPrise-rs.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Input genome file
    #[arg(long, value_name = "FILE")]
    input: String,

    /// Output file name prefix
    #[arg(long, value_name = "PREFIX")]
    output: String,

    /// Match score of the extension alignment
    #[arg(long = "match", value_name = "INT")]
    match_score: Option<i32>,

    /// Mismatch score of the extension alignment
    #[arg(long, value_name = "INT")]
    mismatch: Option<i32>,

    /// Gap open score
    #[arg(long, value_name = "INT")]
    gap: Option<i32>,

    /// Gap extension score
    #[arg(long, value_name = "INT")]
    gapex: Option<i32>,

    /// Penalty of the incomplete length alignment
    #[arg(long = "capplenalty", value_name = "INT")]
    capplenalty: Option<i32>,

    /// Number of mismatches allowed in inexact seed
    #[arg(long, value_name = "INT")]
    dist: Option<i32>,

    /// Upper limit length of extension in one side direction
    #[arg(long, value_name = "INT")]
    maxextend: Option<i32>,

    /// Maximum number of elements in one repeat family
    #[arg(long, value_name = "INT")]
    maxrepeat: Option<i32>,

    /// Band size of extension alignment
    #[arg(long, value_name = "INT")]
    maxgap: Option<i32>,

    /// Stop after INT consecutive non improvements
    #[arg(long = "stopafter", value_name = "INT")]
    stopafter: Option<i32>,

    /// Minimum length of consensus sequence
    #[arg(long, value_name = "INT")]
    minlength: Option<i32>,

    /// Minimum number of elements belonging to one repeat family
    #[arg(long, value_name = "INT")]
    minfreq: Option<i32>,

    /// Penalty associated with the number of regions to be extended
    #[arg(long, value_name = "INT")]
    minimprovement: Option<i32>,

    /// Interval to match the same seed to avoid tandem repeat
    #[arg(long, value_name = "INT")]
    tandemdist: Option<i32>,

    /// Number of parallel cores
    #[arg(long = "pa", value_name = "INT")]
    pa: Option<i32>,

    /// Verbose output
    #[arg(short = 'v', long)]
    v: bool,

    /// Output files about masked region
    #[arg(long)]
    additonalfile: bool,

    /// k-mer length
    #[arg(long)]
    k: Option<i32>,
}

fn main() {
    let args = Cli::parse_from(normalize_args());
    println!("{:#?}", args);

    // Build sequence from the input FASTA file.
    let data = build_sequence(&args.input).expect("failed to read FASTA");
    println!("sequence length: {}", data.sequence.len());

    // Construct suffix array using the Rust implementation.
    let sa = suffix_array(&data.sequence);
    println!("suffix array length: {}", sa.len());

    // Demo of the newly ported routines
    let k = args.k.unwrap_or(4) as usize;
    let cache = reprise::alg::repeat::store_cache(0, k, &data.sequence, &sa);
    let mut kmers = reprise::alg::repeat::build_sortedkmers(
        k,
        &data.sequence,
        &cache,
        &sa,
        args.minfreq.unwrap_or(1) as usize,
    );
    println!("top kmer freq: {:?}", kmers.peek());
    let families = reprise::alg::repeat::build_repeat_families(
        kmers,
        &data.sequence,
        &cache,
        &sa,
        args.minfreq.unwrap_or(1) as usize,
    );
    println!("repeat families detected: {}", families.len());
}
