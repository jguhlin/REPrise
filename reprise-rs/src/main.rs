use clap::Parser;
use reprise_rs::{cli, reprise};

fn main() {
    let config = cli::Config::parse();
    println!("{:?}", config);

    match reprise::read_fasta(&config.input) {
        Ok(genome) => {
            println!("Loaded sequence length: {}", genome.sequence.len());
            println!("Chromosomes: {}", genome.chrtable.len());

            let k = reprise::default_k(genome.sequence.len(), config.dist as usize);
            println!("Default k: {}", k);

            if let Some((kmer, count)) = reprise::best_kmer(&genome.sequence, k) {
                let kmer_str: String = kmer.iter().map(|&b| reprise_rs::utils::num_to_char(b)).collect();
                println!("Most frequent kmer {} occurs {} times", kmer_str, count);
            }
        }
        Err(e) => eprintln!("Failed to read input: {e}")
    }
}
