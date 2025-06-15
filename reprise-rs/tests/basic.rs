use reprise_rs::reprise::{
    read_fasta, compute_entropy, default_k, build_suffix_array, find_kmer,
    kmer_frequencies, best_kmer, remove_tandem, remove_masked, mask_by_seed,
    chrtracer,
};

#[test]
fn load_test_fasta() {
    let genome = read_fasta("../test/tst.fa").expect("read fasta");
    assert!(!genome.sequence.is_empty());
}

#[test]
fn entropy_and_default_k() {
    let kmer = vec![0u8, 1, 2, 3];
    let e = compute_entropy(&kmer);
    let expected = (0.25f64).ln();
    assert!((e - expected).abs() < 1e-6);

    let k = default_k(1000, 0);
    assert!(k > 0);
}

#[test]
fn suffix_array_search() {
    let seq = b"ACGACGTT".to_vec();
    let sa = build_suffix_array(&seq);
    let occ = find_kmer(&seq, &sa, b"ACG");
    assert_eq!(occ, vec![0, 3]);
}

#[test]
fn kmer_count_and_best() {
    let seq = b"ACGACGTT".to_vec();
    let counts = kmer_frequencies(&seq, 3);
    assert_eq!(counts.get(&b"ACG".to_vec()), Some(&2usize));

    let best = best_kmer(&seq, 3).unwrap();
    assert_eq!(best.0, b"ACG".to_vec());
    assert_eq!(best.1, 2);
}

#[test]
fn tandem_and_masking() {
    let mut occs = vec![1usize, 5, 6, 15];
    remove_tandem(&mut occs, 5);
    assert_eq!(occs, vec![1, 6, 15]);

    let mask = vec![false, false, true, false, false, false];
    let mut occs2 = vec![0usize, 1, 4];
    remove_masked(&mut occs2, &mask, false, 2);
    assert_eq!(occs2, vec![0, 4]);

    let mut mask_vec = vec![false; 6];
    mask_by_seed(&occs2, &mut mask_vec, false, 2);
    assert!(mask_vec[0] && mask_vec[1]);
    assert!(!mask_vec[2]);
}

#[test]
fn chromosome_tracking() {
    let genome = read_fasta("../test/tst.fa").expect("read fasta");
    let start = genome.chrtable[1].1;
    let (chr, pos) = chrtracer(&genome.chrtable, start + 10);
    assert_eq!(chr, genome.chrtable[1].0);
    assert_eq!(pos, genome.chrtable[1].1);
}
