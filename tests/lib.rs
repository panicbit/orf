extern crate transcriptome_translation;

use std::io::prelude::*;
use std::fs::File;

#[test]
fn simple() {
    let input = read_file("tests/test-nucleo.FASTA");
    let expected_result = read_file("tests/results.txt");
    let mut result = Vec::new();

    transcriptome_translation::start_parse(&input, &mut result, 1).expect("parse");

    assert!(result == expected_result, "Result differs from expected result");
}

fn read_file(path: &str) -> Vec<u8> {
    let mut data = Vec::new();
    let mut file = File::open(path).unwrap();
    file.read_to_end(&mut data).unwrap();
    data
}
