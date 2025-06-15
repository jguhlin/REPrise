use clap::Parser;

/// Command line arguments for REPrise-rs.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input genome file
    #[arg(long)]
    input: String,

    /// Output file name prefix
    #[arg(long)]
    output: String,

    /// Match score of the extension alignment
    #[arg(long = "match")]
    match_score: Option<i32>,

    /// Mismatch score of the extension alignment
    #[arg(long)]
    mismatch: Option<i32>,

    /// Gap open score
    #[arg(long)]
    gap: Option<i32>,

    /// Gap extension score
    #[arg(long)]
    gapex: Option<i32>,

    /// Penalty of the incomplete length alignment
    #[arg(long = "cappenalty")]
    cap_penalty: Option<i32>,

    /// Number of mismatches allowed in inexact seed
    #[arg(long)]
    dist: Option<i32>,

    /// Upper limit length of extension in one side direction
    #[arg(long)]
    maxextend: Option<i32>,

    /// Maximum number of elements in one repeat family
    #[arg(long)]
    maxrepeat: Option<i32>,

    /// Band size of extension alignment
    #[arg(long)]
    maxgap: Option<i32>,

    /// Stop after INT consecutive non improvements
    #[arg(long = "stopafter")]
    stop_after: Option<i32>,

    /// Minimum length of consensus sequence
    #[arg(long)]
    minlength: Option<i32>,

    /// Minimum number of elements belonging to one repeat family
    #[arg(long)]
    minfreq: Option<i32>,

    /// Penalty associated with the number of regions to be extended
    #[arg(long)]
    minimprovement: Option<i32>,

    /// Interval to match the same seed to avoid tandem repeat
    #[arg(long)]
    tandemdist: Option<i32>,

    /// Number of parallel cores
    #[arg(long = "pa")]
    parallel: Option<i32>,

    /// Verbose output
    #[arg(long)]
    verbose: bool,

    /// Output files about masked region
    #[arg(long)]
    additonalfile: bool,

    /// k-mer length
    #[arg(long)]
    k: Option<i32>,
}

fn main() {
    let args = Args::parse();
    println!("{:#?}", args);
}
