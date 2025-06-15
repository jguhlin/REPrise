use clap::Parser;

#[derive(Parser, Debug)]
#[command(author, version, about)]
pub struct Config {
    /// input file name
    #[arg(long)]
    pub input: String,

    /// output file name prefix
    #[arg(long)]
    pub output: String,

    /// match score
    #[arg(long, default_value_t = 1)]
    pub match_score: i32,

    /// mismatch score
    #[arg(long, default_value_t = -1)]
    pub mismatch: i32,

    /// gap open score
    #[arg(long, default_value_t = -5)]
    pub gap: i32,

    /// gap extension score
    #[arg(long, default_value_t = -1)]
    pub gapex: i32,

    /// penalty of incomplete length alignment
    #[arg(long, default_value_t = -20)]
    pub cappenalty: i32,

    /// number of mismatches allowed in seed
    #[arg(long, default_value_t = 0)]
    pub dist: i32,

    /// maximum extend length
    #[arg(long, default_value_t = 10000)]
    pub maxextend: i32,

    /// maximum elements of a repeat family
    #[arg(long, default_value_t = 100000)]
    pub maxrepeat: i32,

    /// band size of extension alignment
    #[arg(long, default_value_t = 5)]
    pub maxgap: i32,

    /// stop extension after consecutive unchanged scores
    #[arg(long, default_value_t = 100)]
    pub stopafter: i32,

    /// minimum consensus length
    #[arg(long, default_value_t = 50)]
    pub minlength: i32,

    /// minimum frequency
    #[arg(long, default_value_t = 3)]
    pub minfreq: i32,

    /// minimum improvement threshold
    #[arg(long, default_value_t = 3)]
    pub minimprovement: i32,

    /// interval to avoid tandem repeat seeds
    #[arg(long, default_value_t = 500)]
    pub tandemdist: i32,

    /// number of parallel threads
    #[arg(long, default_value_t = 1)]
    pub pa: i32,

    /// enable verbose output
    #[arg(long, default_value_t = false)]
    pub verbose: bool,

    /// output masked region files (.masked and .bed)
    #[arg(long, default_value_t = false)]
    pub additonalfile: bool,
}

