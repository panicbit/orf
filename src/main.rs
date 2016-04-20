extern crate transcriptome_translation;
extern crate memmap;
extern crate num_cpus;

use std::fs::File;
use memmap::{Mmap, Protection};

fn main() {
    let mmap = Mmap::open_path("test/test-nucleo.FASTA", Protection::Read).unwrap();
    let input = unsafe { mmap.as_slice() };
    let mut output_file = File::create("results.txt").unwrap();
    let n_threads = num_cpus::get() as u32;

    transcriptome_translation::start_parse(input, &mut output_file, n_threads);

    output_file.sync_all();
}
