extern crate transcriptome_translation;

use std::fs::OpenOptions;

fn main() {
    let input = "test/test-nucleo.FASTA";
    let mut output_file = OpenOptions::new()
        .create(true)
        .read(false)
        .write(true)
        .open("results.txt")
        .unwrap();
    let n_threads = 4;

    transcriptome_translation::start_parse(input, &mut output_file, n_threads);

    output_file.sync_all();
}
