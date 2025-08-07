#!/usr/bin/env python3
"""Compare REPrise C++ and Rust outputs"""

import sys
import os

def parse_freq_file(filename):
    """Parse a .freq file and return statistics"""
    if not os.path.exists(filename):
        return None
    
    families = []
    total_occurrences = 0
    
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                kmer = parts[0]
                freq = int(parts[1])
                families.append((kmer, freq))
                total_occurrences += freq
    
    return {
        'num_families': len(families),
        'total_occurrences': total_occurrences,
        'top_families': families[:10] if families else [],
        'freq_distribution': sorted([f[1] for f in families], reverse=True)[:20]
    }

def main():
    print("=" * 60)
    print("REPrise Output Comparison")
    print("=" * 60)
    
    # Analyze C++ output on full E. coli
    cpp_stats = parse_freq_file('REPrise-rs/comparison_cpp.freq')
    if cpp_stats:
        print("\nC++ REPrise (Full E. coli):")
        print(f"  - Repeat families detected: {cpp_stats['num_families']:,}")
        print(f"  - Total occurrences: {cpp_stats['total_occurrences']:,}")
        print(f"  - Top frequencies: {cpp_stats['freq_distribution'][:10]}")
        print(f"  - Average frequency: {cpp_stats['total_occurrences']/cpp_stats['num_families']:.1f}")
    
    # Check for Rust outputs
    rust_files = [
        'REPrise-rs/test_rust_ecoli.freq',
        'REPrise-rs/ecoli_test_inexact.freq',
        'subset_rust.freq'
    ]
    
    for rust_file in rust_files:
        if os.path.exists(rust_file):
            rust_stats = parse_freq_file(rust_file)
            if rust_stats:
                print(f"\nRust implementation ({rust_file}):")
                print(f"  - Repeat families detected: {rust_stats['num_families']:,}")
                print(f"  - Total occurrences: {rust_stats['total_occurrences']:,}")
                if rust_stats['num_families'] > 0:
                    print(f"  - Average frequency: {rust_stats['total_occurrences']/rust_stats['num_families']:.1f}")
    
    # Compare subset runs
    print("\n" + "=" * 60)
    print("Subset Comparison (70KB of E. coli):")
    print("=" * 60)
    
    cpp_subset = parse_freq_file('subset_cpp.freq')
    if cpp_subset:
        print(f"\nC++ on subset:")
        print(f"  - Families: {cpp_subset['num_families']}")
        print(f"  - Total occurrences: {cpp_subset['total_occurrences']}")
        print(f"  - Top 5 frequencies: {cpp_subset['freq_distribution'][:5]}")

if __name__ == '__main__':
    main()