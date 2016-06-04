#![feature(test)]
extern crate orf;
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
            orf::start_parse(input, &mut output, n_threads).expect("parse");
        }
    });
}

#[bench]
fn bench_many_small_trans(b: &mut Bencher) {
    let input = include_bytes!("../benches/test-nucleo_many_small.FASTA");
    bench_trans(b, input);
}

#[bench]
fn bench_many_big_trans(b: &mut Bencher) {
    let input = include_bytes!("../benches/test-nucleo_many_big.FASTA");
    bench_trans(b, input);
}

#[bench]
fn bench_single_huge(b: &mut Bencher) {
    let input = include_bytes!("../benches/test-nucleo_single_big.FASTA");
    bench_trans(b, input);
}

fn bench_trans(b: &mut Bencher, input: &[u8]) {
    let mut output = sink();
    let n_threads = num_cpus::get() as u32;
    b.iter(|| {
        orf::start_parse(input, &mut output, n_threads).expect("parse");
    });
}
