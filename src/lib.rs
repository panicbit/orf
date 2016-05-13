#![feature(plugin)]
#![feature(type_macros)]
#![plugin(phf_macros)]
extern crate phf;

#[macro_use]
extern crate nom;
extern crate scoped_threadpool;
extern crate monster;

mod codons;
mod translator;

use std::str;
use std::io::prelude::*;
use nom::{not_line_ending,line_ending};
use nom::IResult;
use std::sync::Arc;
use scoped_threadpool::Pool as ThreadPool;
use std::sync::mpsc::channel;

#[derive(Debug)]
struct FASTA<'a> {
    pub id: &'a str,
    pub sequence: Vec<&'a str>
}

pub fn start_parse<Output>(input: &[u8], mut output: Output, n_threads: u32) where
    Output: Write
{
    let mut threadpool = ThreadPool::new(n_threads);
    let (tx, rx) = channel();
    if let IResult::Done(_,o) = fasta_deserialize(input) {
        threadpool.scoped(|threadpool| {
            for fasta in o {
                let tx = tx.clone();
                let amino_seq: Arc<Vec<u8>> = Arc::new(fasta.sequence
                    .iter()
                    .fold(Vec::new(), |mut acc, item| {
                            acc.extend(item.as_bytes());
                            acc
                    }));
                let fasta_id = fasta.id;
                let dispatch_decoding = |window, decoder: fn(&[u8]) -> String| {
                    let amino_seq = amino_seq.clone();
                    let tx = tx.clone();
                    threadpool.execute(move || {
                        tx.send(FASTA_Complete::new(window, fasta_id, decoder(&amino_seq)));
                    });
                };
                dispatch_decoding("> No Move|", translator::no_move);
                dispatch_decoding("> Shift Left One|", translator::nucleotide_shift_left_one);
                dispatch_decoding("> Shift Left Two|", translator::nucleotide_shift_left_two);
                dispatch_decoding("> Rev. No Move|", translator::rev_no_move);
                dispatch_decoding("> Rev. Shift Left One|", translator::rev_nucleotide_shift_left_one);
                dispatch_decoding("> Rev. Shift Left Two|", translator::rev_nucleotide_shift_left_two);
            }
        });

        drop(tx);

        for result in rx {
            output.write(result.window.as_bytes());
            output.write(result.id.as_bytes());
            output.write(result.sequence.as_bytes());
        }

    }
}

#[derive(Debug)]
struct FASTA_Complete<'a> {
    window: &'a str,
    id: &'a str,
    sequence: String,
}

impl<'a> FASTA_Complete<'a> {
    fn new(window: &'a str, id: &'a str, sequence: String) -> FASTA_Complete<'a> {
        FASTA_Complete {
            window: window,
            id: id,
            sequence: sequence
        }
    }
}

//FASTA_Complete.window are hardcoded to include labeling the id with '>
//and their reading frame.
//There is probably substantial room for improvement and the code could be deduplicated
//After memmapping the file.  I personally prefer laying it all out even if it does get
//a bit lengthy.
fn fasta_deserialize(input:&[u8]) -> IResult<&[u8], Vec<FASTA>>  {
    many0!(input,
      chain!(
        tag!(">") ~
        id: map_res!(not_line_ending, str::from_utf8) ~ line_ending ~
        sequence: many0!(terminated!(map_res!( is_not!(">\n"), str::from_utf8), tag!("\n"))),
        ||{
            FASTA {
                id: id,
                sequence: sequence
            }
        }
      )
   )
}
