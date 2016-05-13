#![feature(question_mark)]
extern crate transcriptome_translation;
extern crate memmap;
extern crate num_cpus;

use std::fs::File;
use std::process::exit;
use memmap::{Mmap, Protection};

fn main() {
    if let Err(e) = transcribe() {
        println!("Error: {}", e);
        exit(1);
    }
}

fn transcribe() -> transcriptome_translation::Result {
    let mmap = Mmap::open_path("tests/test-nucleo.FASTA", Protection::Read)?;
    let input = unsafe { mmap.as_slice() };
    let mut output_file = File::create("results.txt")?;
    let n_threads = num_cpus::get() as u32;

    transcriptome_translation::start_parse(input, &mut output_file, n_threads)?;

    output_file.sync_all()?;

    Ok(())
}