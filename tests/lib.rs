extern crate transcriptome_translation;

use std::io::prelude::*;
use std::fs::File;

#[test]
fn simple() {
    let mut test_data = Vec::new();
    let mut test_file = File::open("tests/results.txt").unwrap();
    test_file.read_to_end(&mut test_data).unwrap();

    let mut result = Vec::new();
    transcriptome_translation::start_parse("tests/test-nucleo.FASTA", &mut result, 1);

    assert!(result == test_data, "Result differs from test data");

}
