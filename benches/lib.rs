#![feature(test)]
extern crate transcriptome_translation;
extern crate test;
extern crate num_cpus;

use test::Bencher;
use std::io::sink;

#[bench]
fn bench_100_simple_trans(b: &mut Bencher) {
    let input = include_bytes!("../tests/test-nucleo.FASTA");
    let mut output = sink();
    let n_threads = num_cpus::get() as u32;

    b.iter(|| {
        for _ in 0 .. 100 {
            transcriptome_translation::start_parse(input, &mut output, n_threads);
        }
    });
}
