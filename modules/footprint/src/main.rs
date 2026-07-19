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

use collapse::{cli::Args, utils::__write_collapsed, utils::unpack};

use clap::Parser;
use simple_logger::init_with_level;

fn main() {
    let start = std::time::Instant::now();

    let args = Args::parse();
    init_with_level(args.level).unwrap_or_else(|e| log::error!("ERROR: Logger init error -> {e}"));

    let target = args.outdir.join("collapsed");
    std::fs::create_dir_all(&target)
        .unwrap_or_else(|e| panic!("ERROR: Could not create directory {:?} -> {e}", &target));

    let tracks = unpack(
        args.bed.clone(),
        args.collapse_mode,
        args.max_five_utr_length,
        args.max_three_utr_length,
    )
    .unwrap_or_else(|e| {
        panic!(
            "ERROR: There was an error while unpacking {:?} -> {e}",
            &args.bed
        )
    });

    let output = target.join(args.name);
    log::info!("INFO: Writing to collapsed file to {:?}", &output);

    __write_collapsed(tracks, &output); // WARN: extra column appended!

    let elapsed = start.elapsed();
    log::info!("INFO: Elapsed time: {elapsed:?}");
}
