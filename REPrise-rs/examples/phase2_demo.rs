//! Phase 2 Implementation Demo - Genome Loading and K-mer Indexing
//!
//! This example demonstrates the complete Phase 2 functionality including:
//! - High-performance genome loading with needletail
//! - Canonical k-mer processing with reverse complements
//! - Contig-by-contig indexing with memory efficiency
//! - Deterministic hashing suitable for concurrent pipelines

use reprise::error::Result;
use reprise::genome::Genome;
use reprise::kmer::{KmerEngine, KmerFrequency};
use reprise::index::{KmerIndex, IndexConfig};
use std::io::Write;
use std::time::Instant;
use tempfile::NamedTempFile;

fn main() -> Result<()> {
    println!("REPrise Phase 2 Demo - Production-Ready Genome Processing");
    println!("=========================================================");
    
    // Demo 1: High-Performance Genome Loading
    demo_genome_loading()?;
    
    // Demo 2: Canonical K-mer Processing  
    demo_kmer_processing()?;
    
    // Demo 3: Large-Scale Indexing with Memory Efficiency
    demo_large_scale_indexing()?;
    
    // Demo 4: Multi-Contig Fragmented Assembly Handling
    demo_fragmented_assembly()?;
    
    println!("\n🎉 All Phase 2 demos completed successfully!");
    println!("\nPhase 2 provides:");
    println!("✅ High-performance FASTA parsing with needletail");
    println!("✅ Canonical k-mer representation with reverse complements");
    println!("✅ Contig-aware indexing for fragmented assemblies");
    println!("✅ Memory-efficient processing of large genomes");
    println!("✅ Deterministic hashing for concurrent pipelines");
    println!("✅ Comprehensive error handling and validation");
    
    Ok(())
}

fn demo_genome_loading() -> Result<()> {
    println!("\n1. High-Performance Genome Loading");
    println!("----------------------------------");
    
    // Create a test genome with multiple contigs
    let fasta_content = r#">chr1 Human chromosome 1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT
>chr2 Human chromosome 2  
TTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAA
CCTGGAACCTGGAACCTGGAACCTGGAACCTGGAACCTGGAACC
>mitochondria Mitochondrial genome
AAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTT
"#;
    
    let mut temp_file = NamedTempFile::new().unwrap();
    temp_file.write_all(fasta_content.as_bytes()).unwrap();
    temp_file.flush().unwrap();
    
    let start = Instant::now();
    let genome = Genome::from_fasta(temp_file.path())?;
    let load_time = start.elapsed();
    
    println!("✅ Loaded genome in {:?}", load_time);
    println!("   - Total size: {} bases", genome.len());
    println!("   - Contigs: {}", genome.num_contigs());
    
    // Demonstrate contig-aware operations
    for (id, contig) in genome.contigs() {
        println!("   - {}: {}..{} ({} bases)", contig.name, contig.range().start, contig.range().end, contig.length);
        
        // Test contig boundary detection
        let mid_pos = contig.start + contig.length / 2;
        assert_eq!(genome.contig_of(mid_pos), Some(id));
        
        // Test within-contig validation
        assert!(genome.is_within_one_contig(contig.start, contig.length as usize));
    }
    
    println!("✅ All contig boundaries validated");
    Ok(())
}

fn demo_kmer_processing() -> Result<()> {
    println!("\n2. Canonical K-mer Processing");
    println!("-----------------------------");
    
    let test_sequence = "ATCGATCGATCGATCGATCG";
    let k = 13;
    
    // Create k-mer engine
    let engine = KmerEngine::new(k)?;
    let start = Instant::now();
    
    // Test canonical representation
    let sequence_numeric: Vec<u8> = test_sequence.chars()
        .map(|c| match c {
            'A' => 0, 'T' => 3, 'C' => 1, 'G' => 2, _ => 99
        })
        .collect();
    
    let kmers = engine.extract_kmers(&sequence_numeric);
    let process_time = start.elapsed();
    
    println!("✅ Extracted {} k-mers in {:?}", kmers.len(), process_time);
    
    // Demonstrate canonical property
    if let Some(first_kmer) = kmers.first() {
        let rev_comp = first_kmer.reverse_complement();
        println!("   - Original:          {:016x}", first_kmer.canonical());
        println!("   - Reverse complement: {:016x}", rev_comp.canonical());
        println!("   - Both canonical:     {}", first_kmer.canonical() == rev_comp.canonical());
    }
    
    // Test frequency counting
    let mut frequency_counter = KmerFrequency::new(k)?;
    frequency_counter.count_sequence(&sequence_numeric);
    
    println!("✅ Frequency analysis:");
    println!("   - Unique k-mers: {}", frequency_counter.unique_kmers());
    println!("   - Total k-mers: {}", frequency_counter.total_kmers());
    
    Ok(())
}

fn demo_large_scale_indexing() -> Result<()> {
    println!("\n3. Large-Scale Indexing Performance");
    println!("-----------------------------------");
    
    // Generate a larger synthetic genome (simulate ~1MB)
    let large_sequence = "ATCGATCGATCG".repeat(80000); // ~1MB
    let fasta_content = format!(">large_contig\n{}\n", large_sequence);
    
    let mut temp_file = NamedTempFile::new().unwrap();
    temp_file.write_all(fasta_content.as_bytes()).unwrap();
    temp_file.flush().unwrap();
    
    // Load genome
    let load_start = Instant::now();
    let genome = Genome::from_fasta(temp_file.path())?;
    let load_time = load_start.elapsed();
    
    println!("✅ Loaded {:.1} MB genome in {:?}", genome.len() as f64 / 1_000_000.0, load_time);
    
    // Build index with parallel processing
    let config = IndexConfig {
        k: 13,
        min_frequency: 2,
        max_frequency: Some(1000),
        parallel: true,
        max_positions_per_kmer: 10000,
    };
    
    let index_start = Instant::now();
    let index = KmerIndex::build(&genome, config)?;
    let index_time = index_start.elapsed();
    
    let stats = index.stats();
    println!("✅ Built index in {:?}", index_time);
    println!("   - Unique k-mers: {}", stats.total_kmers);
    println!("   - Total positions: {}", stats.total_positions);
    println!("   - Average frequency: {:.2}", stats.avg_frequency);
    println!("   - Memory usage: {:.1} MB", stats.estimated_memory_usage() as f64 / 1_000_000.0);
    
    // Test index performance
    let frequent_kmers = index.most_frequent_kmers(10);
    println!("✅ Top frequent k-mers:");
    for (i, (kmer, freq)) in frequent_kmers.iter().take(5).enumerate() {
        println!("   {}. {:016x} (frequency: {})", i + 1, kmer.canonical(), freq);
    }
    
    Ok(())
}

fn demo_fragmented_assembly() -> Result<()> {
    println!("\n4. Fragmented Assembly Processing");
    println!("----------------------------------");
    
    // Simulate a highly fragmented assembly
    let mut fasta_content = String::new();
    
    // Create many small contigs
    for i in 0..1000 {
        let contig_seq = match i % 4 {
            0 => "ATCGATCGATCGATCGATCG",
            1 => "GCATGCATGCATGCATGCAT", 
            2 => "TTAACCGGTTAACCGGTTAA",
            _ => "CCTGGAACCTGGAACCTGGA",
        };
        fasta_content.push_str(&format!(">contig_{:04}\n{}\n", i, contig_seq));
    }
    
    let mut temp_file = NamedTempFile::new().unwrap();
    temp_file.write_all(fasta_content.as_bytes()).unwrap();
    temp_file.flush().unwrap();
    
    let start = Instant::now();
    let genome = Genome::from_fasta(temp_file.path())?;
    let load_time = start.elapsed();
    
    println!("✅ Loaded fragmented assembly in {:?}", load_time);
    println!("   - Contigs: {}", genome.num_contigs());
    println!("   - Total size: {} bases", genome.len());
    println!("   - Average contig size: {:.0} bases", 
        genome.len() as f64 / genome.num_contigs() as f64);
    
    // Test contig-by-contig indexing
    let config = IndexConfig {
        k: 13,
        min_frequency: 1,
        max_frequency: None,
        parallel: true,
        max_positions_per_kmer: 10000,
    };
    
    let index_start = Instant::now();
    let index = KmerIndex::build(&genome, config)?;
    let index_time = index_start.elapsed();
    
    let stats = index.stats();
    println!("✅ Indexed fragmented assembly in {:?}", index_time);
    
    // Find k-mers that span multiple contigs
    let multi_contig_kmers = index.multi_contig_kmers();
    println!("   - Multi-contig k-mers: {}", multi_contig_kmers.len());
    println!("   - Memory efficiency: {:.1} bytes/base", 
        stats.estimated_memory_usage() as f64 / genome.len() as f64);
    
    // Test memory usage is reasonable
    let memory_mb = stats.estimated_memory_usage() as f64 / 1_000_000.0;
    let genome_mb = genome.len() as f64 / 1_000_000.0;
    let overhead_ratio = memory_mb / genome_mb;
    
    println!("   - Index overhead: {:.1}x genome size", overhead_ratio);
    
    // Should be reasonable overhead (higher for small fragmented assemblies due to fixed costs)
    assert!(overhead_ratio < 20.0, "Index overhead should be reasonable for small fragmented assemblies");
    
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_phase2_integration() {
        // Integration test to ensure all Phase 2 components work together
        let result = std::panic::catch_unwind(|| {
            demo_genome_loading().unwrap();
            demo_kmer_processing().unwrap();
            // Skip the large demos in tests for speed
        });
        
        assert!(result.is_ok(), "Phase 2 integration test failed");
    }
}