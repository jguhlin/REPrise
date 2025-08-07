use crate::kmer::Kmer;
use crate::genome::Genome;
use crate::index::KmerIndex;
use crate::mask::{Bitmask, ClaimGuard};
use crate::error::Result;
use crossbeam::channel::{self, Sender, Receiver};
use rayon::prelude::*;
use std::sync::Arc;
use std::ops::Range;

/// Represents a candidate pair of genomic regions for repeat analysis.
pub struct CandidatePair {
    pub region1: GenomicRegion,
    pub region2: GenomicRegion,
    pub seed_kmer: Kmer,
    pub seed_frequency: u32,
    pub priority: u32,
}

/// Represents a genomic region with contig awareness.
pub struct GenomicRegion {
    pub position: u64,
    pub contig: u32,
    pub local_position: u64,
    pub length: u64,
}

impl GenomicRegion {
    pub fn range(&self) -> Range<u64> {
        self.position..self.position + self.length
    }
}

/// Configuration for the repeat detection pipeline.
pub struct PipelineConfig {
    pub channel_capacity: usize,
    pub num_workers: usize,
    pub extension_length: u64,
    pub min_seed_frequency: u32,
    pub max_seed_frequency: Option<u32>,
    pub max_memory_usage: u64,
    pub backpressure_threshold: f64,
    pub enable_inter_contig: bool,
    pub max_candidates: Option<u64>,
}

impl Default for PipelineConfig {
    fn default() -> Self {
        Self {
            channel_capacity: 1_000_000,
            num_workers: num_cpus::get(),
            extension_length: 1000,
            min_seed_frequency: 5,
            max_seed_frequency: Some(10000),
            max_memory_usage: u64::MAX,
            backpressure_threshold: 0.8,
            enable_inter_contig: true,
            max_candidates: None,
        }
    }
}

/// Builder for creating a Pipeline instance.
pub struct PipelineBuilder {
    config: PipelineConfig,
}

impl PipelineBuilder {
    pub fn new() -> Self {
        Self { config: PipelineConfig::default() }
    }

    pub fn channel_capacity(mut self, capacity: usize) -> Self {
        self.config.channel_capacity = capacity;
        self
    }

    pub fn num_workers(mut self, num_workers: usize) -> Self {
        self.config.num_workers = num_workers;
        self
    }

    pub fn extension_length(mut self, length: u64) -> Self {
        self.config.extension_length = length;
        self
    }

    pub fn seed_frequency_range(mut self, min: u32, max: Option<u32>) -> Self {
        self.config.min_seed_frequency = min;
        self.config.max_seed_frequency = max;
        self
    }

    pub fn enable_inter_contig(mut self, enable: bool) -> Self {
        self.config.enable_inter_contig = enable;
        self
    }

    pub fn build(self) -> Pipeline {
        Pipeline::with_config(self.config)
    }
}

pub struct DetectedRepeat {
    pub family_id: u32,
    pub start: u64,
    pub end: u64,
    pub score: f64,
    pub identity: f64,
    pub length: u64,
}

pub struct Pipeline {
    config: PipelineConfig,
}

impl Pipeline {
    pub fn with_config(config: PipelineConfig) -> Self {
        Self { config }
    }

    pub fn run(
        self,
        genome: Arc<Genome>,
        index: Arc<KmerIndex>,
        mask: Arc<Bitmask>,
    ) -> Result<Vec<DetectedRepeat>> {
        let (tx, rx): (Sender<CandidatePair>, Receiver<CandidatePair>) = 
            channel::bounded(self.config.channel_capacity);

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.config.num_workers)
            .build()?;

        let results = pool.install(|| {
            let (producer_res, consumer_res) = rayon::join(
                || self.producer(tx, index, &mask),
                || self.consumer(rx, genome, &mask),
            );
            producer_res?; // Propagate producer errors
            consumer_res
        })?;

        Ok(results)
    }

    fn producer(&self, tx: Sender<CandidatePair>, index: Arc<KmerIndex>, mask: &Bitmask) -> Result<()> {
        let candidates = index.get_candidate_pairs(self.config.min_seed_frequency, self.config.max_seed_frequency);
        for candidate in candidates {
            if !mask.is_range_free(&candidate.region1.range()) || !mask.is_range_free(&candidate.region2.range()) {
                continue;
            }
            tx.send(candidate)?;
        }
        Ok(())
    }

    fn consumer(&self, rx: Receiver<CandidatePair>, genome: Arc<Genome>, mask: &Bitmask) -> Result<Vec<DetectedRepeat>> {
        let mut repeats = Vec::new();
        while let Ok(candidate) = rx.recv() {
            let guard1 = ClaimGuard::new(mask, candidate.region1.range());
            if guard1.is_some() {
                let guard2 = ClaimGuard::new(mask, candidate.region2.range());
                if guard2.is_some() {
                    // Successfully claimed both regions
                    // In a real implementation, perform alignment here
                    repeats.push(DetectedRepeat {
                        family_id: 0,
                        start: candidate.region1.position,
                        end: candidate.region1.position + candidate.region1.length,
                        score: 0.0,
                        identity: 0.0,
                        length: candidate.region1.length,
                    });
                }
            }
        }
        Ok(repeats)
    }
}
