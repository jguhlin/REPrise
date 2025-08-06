# TODO: Complete Feature Parity with C++ REPrise

## Priority 1: Fix Algorithmic Discrepancies

### 1.1 Fix Seed Masking Coordination
**Issue**: Rust finds more repeat families than C++ (6 vs 1 on simple test, ~1500 vs 45 on E. coli)
**Root Cause**: Overlapping k-mers from the same repeat region are being selected as separate seeds

**Implementation**:
- **File**: `src/alg/repeat.rs` 
- **Function**: `find_bestseed()`
- **Changes**:
  - After selecting a seed, immediately mask ALL occurrences before heap is re-examined
  - Ensure k-mers overlapping with already-masked regions are skipped
  - Add check: when popping from heap, verify k-mer isn't already masked

```rust
// In find_bestseed, after selecting best seed:
// 1. Mask all forward occurrences immediately
// 2. Mask all reverse occurrences immediately  
// 3. Re-validate heap entries against updated mask
```

### 1.2 Optimize Inexact Seeding Performance
**Issue**: Inexact seeding with edit distance is extremely slow on large genomes
**Root Cause**: Pre-computing all k-mer variants is computationally expensive

**Implementation**:
- **File**: `src/alg/repeat.rs`
- **Function**: `store_cache()`
- **Options**:
  1. Lazy evaluation: Only compute inexact matches when needed
  2. Parallel processing: Use rayon for cache generation
  3. Limit inexact matching to high-frequency k-mers only

## Priority 2: Add Missing Output Files

### 2.1 Implement .freq File Output
**Purpose**: Frequency and statistics for each repeat family
**Format**: Tab-delimited with columns: repeat_id, frequency, occurrences

**Implementation**:
- **File**: `src/main.rs`
- **Location**: After processing each repeat family (around line 374)
- **Code Structure**:
```rust
// Open freq file at start
let mut freq_writer = File::create(format!("{}.freq", args.output))?;

// For each repeat family:
writeln!(freq_writer, "R={}\t{}\t{}", repeat_num, seedfreq, element_count)?;
```

### 2.2 Implement .bed File Output  
**Purpose**: BED format genomic coordinates for visualization in genome browsers
**Format**: chromosome, start, end, name, score, strand

**Implementation**:
- **File**: `src/main.rs` 
- **Dependencies**: Need chromosome tracing functionality first (see 2.4)
- **Code Structure**:
```rust
// Open BED file at start
let mut bed_writer = File::create(format!("{}.bed", args.output))?;

// For each element in repeat family:
let (chr_name, chr_offset) = get_chromosome_info(element_start, &data.chrtable);
writeln!(bed_writer, "{}\t{}\t{}\tR={}\t{}\t{}", 
    chr_name, 
    element_start - chr_offset,
    element_end - chr_offset,
    repeat_num,
    element_length,
    if is_reverse { "-" } else { "+" }
)?;
```

### 2.3 Implement .masked File Output
**Purpose**: Output genome sequence with repeat regions masked (lowercase or 'N')
**Format**: FASTA with repeats in lowercase

**Implementation**:
- **File**: `src/main.rs`
- **Location**: After all repeat families processed
- **Code Structure**:
```rust
// After all families processed:
let mut masked_writer = File::create(format!("{}.masked", args.output))?;

// Write masked sequence
for (chr_name, chr_start) in &data.chrtable {
    writeln!(masked_writer, ">{}", chr_name)?;
    
    let chr_end = get_next_chr_start(chr_name, &data.chrtable);
    for i in chr_start..chr_end {
        let base = if mask[i] {
            num_to_char(data.sequence[i]).to_ascii_lowercase()
        } else {
            num_to_char(data.sequence[i])
        };
        write!(masked_writer, "{}", base)?;
        if (i - chr_start + 1) % 80 == 0 {
            writeln!(masked_writer)?;
        }
    }
}
```

### 2.4 Add Chromosome Position Tracking
**Purpose**: Convert absolute positions to chromosome-relative coordinates

**Implementation**:
- **File**: `src/lib.rs`
- **New Function**: `chrtracer()`
```rust
pub fn chrtracer(pos: usize, chrtable: &[(String, usize)]) -> (String, usize) {
    // Binary search or iterate to find chromosome
    for i in (0..chrtable.len()).rev() {
        if pos >= chrtable[i].1 {
            return (chrtable[i].0.clone(), chrtable[i].1);
        }
    }
    ("unknown".to_string(), 0)
}
```

### 2.5 Add -additionalfile Flag
**Purpose**: Control whether .bed and .masked files are generated

**Implementation**:
- **File**: `src/main.rs`
- **Changes**:
  1. Add field to `Cli` struct: `additional_file: bool`
  2. Add CLI parsing: `"-additionalfile" => { additional_file = true; }`
  3. Conditionally create .bed and .masked files based on flag

## Priority 3: Code Quality Improvements

### 3.1 Remove Unused Constants Warning
**File**: `src/alg/repeat.rs`
- Remove unused constants (lines 822-828) or mark with `#[allow(dead_code)]`

### 3.2 Clean Up Unused Imports
**File**: `src/main.rs`
- Remove `reverse_complement` from imports (line 4)
- Remove unused `RepeatHit` struct (line 166)

### 3.3 Add Integration Tests
**File**: `tests/full_equivalence.rs`
- Test all output files (.reprof, .freq, .bed, .masked)
- Compare with C++ outputs character by character
- Test with various parameter combinations

## Testing Strategy

1. **Simple Dataset** (`test/tst.fa`):
   - Verify 1 repeat family found (matching C++)
   - Check all output files generated correctly

2. **E. coli Dataset** (`data/ecoli.fasta`):
   - Target: 45 repeat families (matching C++)
   - Verify performance is acceptable (<60 seconds)
   - Check output file sizes match C++ version

3. **Parameter Testing**:
   - Test with `-dist 1` for inexact matching
   - Test with different `-minfreq` values
   - Test with `-additionalfile` flag

## Implementation Order

1. **Week 1**: Fix seed masking coordination (Priority 1.1)
   - This is the critical bug causing wrong repeat counts
   
2. **Week 1**: Add .freq file output (Priority 2.1)
   - Simple addition, helps with debugging
   
3. **Week 2**: Add chromosome tracking and .bed output (Priority 2.4, 2.2)
   - Required for proper genomic coordinates
   
4. **Week 2**: Add .masked file output (Priority 2.3)
   - Useful for downstream analysis
   
5. **Week 3**: Optimize inexact seeding if needed (Priority 1.2)
   - Only if performance is unacceptable
   
6. **Week 3**: Clean up code warnings (Priority 3)
   - Final polish

## Success Metrics

- [ ] Rust and C++ find same number of repeat families on all test datasets
- [ ] All output files (.reprof, .freq, .bed, .masked) match C++ format
- [ ] Performance within 2x of C++ version
- [ ] No compiler warnings
- [ ] Integration tests pass
- [ ] Can process E. coli genome in under 60 seconds