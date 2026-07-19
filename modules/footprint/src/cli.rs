// Copyright (c) 2026 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Core module for collapsing BED files with deduplication and indexing
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main functions for efficiently collapsing BED files
//! by identifying and deduplicating identical genomic intervals. The module provides
//! flexible output modes including collapsed BED files with queue annotations or
//! separate binary index files for space-efficient storage, extending the original
//! BED file with additional columns for read name and queue information and fast
//! lookups fo read names.
//!
//! In short, every BED entry is fingerprinted using byte-level hashing accounting
//! for specific columns from the original BED file to ensure that only truly identical
//! rows are grouped together. The deduplicated entries are held in memory alongside
//! their corresponding read identifier queues [maintaining original order for reconstruction].

use clap::{self, Parser, ValueEnum};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[clap(name="collapse", author = env!("CARGO_PKG_AUTHORS"), version = env!("CARGO_PKG_VERSION"), about = "shrink your .beds", long_about = None)]
pub struct Args {
    #[arg(
        short = 't',
        long = "threads",
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,

    #[arg(
        short = 'L',
        long = "level",
        help = "Logging level",
        value_name = "LEVEL",
        default_value_t = log::Level::Info,
    )]
    pub level: log::Level,

    #[arg(
        short = 'b',
        long = "bed",
        required = true,
        value_name = "PATHS",
        value_delimiter = ',',
        num_args = 1..,
        help = "Paths to BED12 files delimited by comma"
    )]
    pub bed: Vec<PathBuf>,

    #[arg(
        short = 'c',
        long = "collapse-mode",
        help = "Mode to collapse reads. Affects the way the binnary key is processed and thus any equivalence",
        value_name = "COLLAPSE_MODE",
        default_value = "transcript"
    )]
    pub collapse_mode: CollapseMode,

    #[arg(
        short = 'o',
        long = "outdir",
        help = "Path to output directory (/collapse)",
        value_name = "PATH",
        required = false,
        default_value = "."
    )]
    pub outdir: PathBuf,

    #[arg(
        short = 'n',
        long = "name",
        help = "Name of the output file",
        value_name = "NAME",
        required = false,
        default_value = "collapsed.bed"
    )]
    pub name: String,

    #[arg(
        short = 'U',
        long = "max-five-utr-length",
        help = "Maximum 5'UTR distance between queues to allow collapsing reads in bp",
        value_name = "INT",
        default_value = "0"
    )]
    pub max_five_utr_length: u32,

    #[arg(
        short = 'u',
        long = "max-three-utr-length",
        help = "Maximum 3'UTR distance between queues to allow collapsing reads in bp",
        value_name = "INT",
        default_value = "0"
    )]
    pub max_three_utr_length: u32,
}

#[derive(Parser, Debug, Clone, ValueEnum, Copy)]
pub enum CollapseMode {
    #[clap(name = "exon")]
    Exon, // INFO: cds_start, cds_end are 0
    #[clap(name = "cds")]
    Coding, // INFO: tx_start, tx_end are 0
    #[clap(name = "transcript")]
    Transcript, // INFO: cannonical
    #[clap(name = "gapped-cds")]
    GappedCoding, // INFO: tx_start, tx_end are 0, all gaps preserved
}

impl std::str::FromStr for CollapseMode {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "exon" => Ok(CollapseMode::Exon),
            "cds" => Ok(CollapseMode::Coding),
            "transcript" => Ok(CollapseMode::Transcript),
            "gapped-cds" => Ok(CollapseMode::GappedCoding),
            _ => Err(format!("Invalid collapse mode: {s}")),
        }
    }
}
