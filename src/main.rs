extern crate transcriptome_translation;

use std::fs::File;

fn main() {
    let input = "test/test-nucleo.FASTA";
    let mut output_file = File::create("results.txt").unwrap();
    let n_threads = 4;

    transcriptome_translation::start_parse(input, &mut output_file, n_threads);

    output_file.sync_all();
}
