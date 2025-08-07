// Scaling benchmark for suffix array implementations
use std::time::Instant;
use suffix::SuffixTable;

fn naive_suffix_array(seq: &[u8]) -> Vec<i64> {
    let mut sa: Vec<i64> = (0..seq.len() as i64).collect();
    sa.sort_by(|&a, &b| seq[a as usize..].cmp(&seq[b as usize..]));
    sa
}

fn suffix_crate_array(seq: &[u8]) -> Vec<i64> {
    let text = String::from_utf8_lossy(seq);
    let table = SuffixTable::new(text);
    table.table().iter().map(|&x| x as i64).collect()
}

fn create_random_dna(length: usize) -> Vec<u8> {
    (0..length).map(|i| (i % 4) as u8).collect()
}

fn benchmark_size(size: usize) {
    println!("Testing sequence size: {}", size);
    let seq = create_random_dna(size);
    
    // Test naive implementation (only for smaller sizes)
    let naive_time = if size <= 50000 {
        let start = Instant::now();
        let _naive_sa = naive_suffix_array(&seq);
        Some(start.elapsed())
    } else {
        None
    };
    
    // Test suffix crate implementation
    let start = Instant::now();
    let _suffix_sa = suffix_crate_array(&seq);
    let suffix_time = start.elapsed();
    
    if let Some(naive_time) = naive_time {
        let speedup = naive_time.as_secs_f64() / suffix_time.as_secs_f64();
        println!("  Naive:     {:?}", naive_time);
        println!("  Suffix:    {:?}", suffix_time);
        println!("  Speedup:   {:.2}x", speedup);
    } else {
        println!("  Suffix:    {:?} (naive too slow)", suffix_time);
    }
    println!();
}

fn main() {
    println!("Suffix Array Scaling Benchmark");
    println!("===============================");
    println!();

    let sizes = vec![1000, 5000, 10000, 25000, 50000, 100000, 250000];
    
    for size in sizes {
        benchmark_size(size);
    }
}