#![feature(plugin)]
#![feature(type_macros)]
#![plugin(phf_macros)]
extern crate phf;

#[macro_use]
extern crate nom;
extern crate scoped_threadpool;
extern crate monster;

mod codons;
use codons::{CODONS, REV_CODONS};

use std::str;
use std::io::prelude::*;
use nom::{not_line_ending,line_ending};
use nom::IResult;
use std::vec::*;
use std::path::Path;
use std::sync::Arc;
use std::sync::mpsc::{Sender, Receiver};
use std::sync::mpsc;
use std::thread;
use std::mem;
use scoped_threadpool::Pool as ThreadPool;
use std::sync::mpsc::channel;
use std::rc::Rc;
use monster::incubation::{SliceDropFirst, SliceDropLast};

#[derive(Debug)]
pub struct FASTA<'a> {
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
                let fasta = Arc::new(fasta);
                let dispatch_decoding = |window, decoder: fn(&[u8]) -> String| {
                    let fasta = fasta.clone();
                    let amino_seq = amino_seq.clone();
                    let tx = tx.clone();
                    threadpool.execute(move || {
                        tx.send(FASTA_Complete::new(window, fasta.id, decoder(&amino_seq)));
                    });
                };
                dispatch_decoding("> No Move|", no_move);
                dispatch_decoding("> Shift Left One|", nucleotide_shift_left_one);
                dispatch_decoding("> Shift Left Two|", nucleotide_shift_left_two);
                dispatch_decoding("> Rev. No Move|", rev_no_move);
                dispatch_decoding("> Rev. Shift Left One|", rev_nucleotide_shift_left_one);
                dispatch_decoding("> Rev. Shift Left Two|", rev_nucleotide_shift_left_two);
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
pub struct FASTA_Complete<'a> {
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
pub fn rc_handler(read: FASTA) -> Rc<FASTA> {
    let amino_seq = Rc::new(read);
    amino_seq

}

pub fn rev_nucleotide_shift_left_two(mut amino_seq: &[u8]) -> String {
    // fn rev_nucleotide_shift_left_two does the following:
    // Reverses all elements in the Array.
    // Removes element at position '0'
    // Removes elemtne at position '0'
    // Then we check to see if the vector is a multiple of three
    // IF the vector is not a multiple of three we remove from the end
    // of the vector a single element then check until the vector is a multiple of three.
    // Then we convert our stream of nucleotides into groups of three
    // We then convert our groups of three into utf8 encoded String
    // We then use a phf to convert our utf8 encoded Strings into corrosponding
    // u8 Amino Acid encoding,
    // We then push the results of the Amino Acid encoding to a vector.
    // We push() a newline to the end of the String to assist with file encoding.
    let mut done = String::with_capacity(amino_seq.len() / 3 + 2);;
    done.push('\n');

    // Shift sequence two elements to the right
    amino_seq = amino_seq.drop_last(2);

    rev_trim_and_map(&amino_seq, &mut done);

    done.push('\n');
    done
}

pub fn rev_nucleotide_shift_left_one(mut amino_seq: &[u8]) -> String {
    // fn rev_nucleotide_shift_left_one does the following:
    // Reverses all elements in the Array.
    // Removes element at position '0'
    // Then we check to see if the vector is a multiple of three
    // IF the vector is not a multiple of three we remove from the end
    // of the vector a single element then check until the vector is a multiple of three.
    // Then we convert our stream of nucleotides into groups of three
    // We then convert our groups of three into utf8 encoded String
    // We then use a phf to convert our utf8 encoded Strings into corrosponding
    // u8 Amino Acid encoding,
    // We then push the results of the Amino Acid encoding to a vector.
    // We push() a newline to the end of the String to assist with file encoding.
    let mut done = String::new();
    done.push('\n');

    // Shift elements to the right once
    amino_seq = amino_seq.drop_last(1);

    rev_trim_and_map(&amino_seq, &mut done);

    done.push('\n');
    done
}

pub fn rev_no_move(amino_seq: &[u8]) -> String {
    // fn rev_no_move does the following:
    // Reverses all elements in the Array.
    // Then we check to see if the vector is a multiple of three
    // IF the vector is not a multiple of three we remove from the end
    // of the vector a single element then check until the vector is a multiple of three.
    // Then we convert our stream of nucleotides into groups of three
    // We then convert our groups of three into utf8 encoded String
    // We then use a phf to convert our utf8 encoded Strings into corrosponding
    // u8 Amino Acid encoding,
    // We then push the results of the Amino Acid encoding to a vector.
    // We push() a newline to the end of the String to assist with file encoding.
    let mut done = String::with_capacity(amino_seq.len() / 3 + 2);;
    done.push('\n');

    rev_trim_and_map(&amino_seq, &mut done);

    done.push('\n');
    done
}

pub fn nucleotide_shift_left_two(mut amino_seq: &[u8]) -> String {
    // fn nucleotide_shift_left_two does the following:
    // Removes element at position '0'
    // Removes elemtne at position '0'
    // Then we check to see if the vector is a multiple of three
    // IF the vector is not a multiple of three we remove from the end
    // of the vector a single element then check until the vector is a multiple of three.
    // Then we convert our stream of nucleotides into groups of three
    // We then convert our groups of three into utf8 encoded String
    // We then use a phf to convert our utf8 encoded Strings into corrosponding
    // u8 Amino Acid encoding,
    // We then push the results of the Amino Acid encoding to a vector.
    // We push() a newline to the end of the String to assist with file encoding.
    let mut done = String::with_capacity(amino_seq.len() / 3 + 2);;

    // Shift elements to the left twice
    amino_seq = amino_seq.drop_first(2);

    trim_and_map(&amino_seq, &mut done);

    done.push('\n');
    done
}

pub fn nucleotide_shift_left_one(mut amino_seq: &[u8]) -> String {
    // fn nucleotide_shift_left_one does the following:
    // Removes elemtne at position '0'
    // Then we check to see if the vector is a multiple of three
    // IF the vector is not a multiple of three we remove from the end
    // of the vector a single element then check until the vector is a multiple of three.
    // Then we convert our stream of nucleotides into groups of three
    // We then convert our groups of three into utf8 encoded String
    // We then use a phf to convert our utf8 encoded Strings into corrosponding
    // u8 Amino Acid encoding,
    // We then push the results of the Amino Acid encoding to a vector.
    // We push() a newline to the end of the String to assist with file encoding.
    let mut done = String::with_capacity(amino_seq.len() / 3 + 2);;
    done.push('\n');

    // Shift elements to the left once
    amino_seq = amino_seq.drop_first(1);

    trim_and_map(amino_seq, &mut done);

    done.push('\n');
    done
}

pub fn no_move<'a>(amino_seq: &[u8]) -> String {
    // fn no_move does the following:
    // Then we check to see if the vector is a multiple of three
    // IF the vector is not a multiple of three we remove from the end
    // of the vector a single element then check until the vector is a multiple of three.
    // Then we convert our stream of nucleotides into groups of three
    // We then convert our groups of three into utf8 encoded String
    // We then use a phf to convert our utf8 encoded Strings into corrosponding
    // u8 Amino Acid encoding,
    // We then push the results of the Amino Acid encoding to a vector.
    // We push() a newline to the end of the String to assist with file encoding.
    let mut done = String::with_capacity(amino_seq.len() / 3 + 2);;
    done.push('\n');

    trim_and_map(amino_seq, &mut done);

    done.push('\n');
    done
}

fn trim_and_map(mut amino_seq: &[u8], done: &mut String) {
    // Trim elements from the end until the length is a multiple of 3
    amino_seq = amino_seq.drop_last(amino_seq.len() % 3);
    debug_assert!(amino_seq.len() % 3 == 0);

    for aminos in amino_seq.chunks(3) {
        debug_assert!(aminos.len() == 3);
        match CODONS.get(aminos) {
            Some(&p) => done.push(p),
            None => println!("Done!"),
        }
    }
}

fn rev_trim_and_map(mut amino_seq: &[u8], done: &mut String) {
    // Trim elements from the beginning until the length is a multiple of 3
    amino_seq = amino_seq.drop_first(amino_seq.len() % 3);
    debug_assert!(amino_seq.len() % 3 == 0);

    for aminos in amino_seq.chunks(3).rev() {
        debug_assert!(aminos.len() == 3);
        match REV_CODONS.get(aminos) {
            Some(&p) => done.push(p),
            None => println!("Done!"),
        }
    }
}

pub fn fasta_deserialize(input:&[u8]) -> IResult<&[u8], Vec<FASTA>>  {
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
