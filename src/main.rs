extern crate transcriptome_translation;

fn main() {
    let input = "test/test-nucleo.FASTA";
    let output = "results.txt";
    let n_threads = 4;
    transcriptome_translation::start_parse(input, output, n_threads);
}
