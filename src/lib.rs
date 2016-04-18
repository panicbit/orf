#![feature(plugin)]
#![feature(type_macros)]
#![plugin(phf_macros)]
extern crate phf;

#[macro_use]
extern crate nom;
extern crate memmap;
extern crate scoped_threadpool;
extern crate monster;

use std::str;
use std::io::prelude::*;
use nom::{not_line_ending,line_ending};
use nom::IResult;
use memmap::{Mmap, Protection};
use std::vec::*;
use std::path::Path;
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

static CODONS: phf::Map<&'static [u8], char> = phf_map! {
    //Alanine
    b"GCA" => 'A',
    b"GCG" => 'A',
    b"GCC" => 'A',
    b"GCT" => 'A',
    //Aspartic_Acid (D)
    //Asparagine (N)
    //Cysteine
    b"TGT" => 'C',
    b"TGC" => 'C',
    //Aspartic_Acid
    b"GAC" => 'D',
    b"GAT" => 'D',
    //Glutamic_Acid
    b"GAA" => 'E',
    b"GAG" => 'E',
    //Phenylalanine
    b"TTT" => 'F',
    b"TTC" => 'F',
    //Glycine
    b"GGA" => 'G',
    b"GGG" => 'G',
    b"GGC" => 'G',
    b"GGT" => 'G',
    //Histidine
    b"CAC" => 'H',
    b"CAT" => 'H',
    //Isoleucine
    b"ATT" => 'I',
    b"ATC" => 'I',
    b"ATA" => 'I',
    //Leucine (L)
    b"TTG" => 'L',
    b"TTA" => 'L',
    b"CTA" => 'L',
    b"CTC" => 'L',
    b"CTG" => 'L',
    b"CTT" => 'L',
    //Lysine (K)
    b"AAA" => 'K',
    b"AAG" => 'K',
    //Methionine (M)
    b"ATG" => 'M',
    //Asparagine (N)
    b"AAT" => 'N',
    b"AAC" => 'N',
    //Pyrrolysine (O) Special Stop Codon
    b"UAG" => 'O',
    //Proline (P)
    b"CCA" => 'P',
    b"CCG" => 'P',
    b"CCC" => 'P',
    b"CCT" => 'P',
    //Glutamine (Q)
    b"CAA" => 'Q',
    b"CAG" => 'Q',
    //Arginine (R)
    b"AGA" => 'R',
    b"AGG" => 'R',
    b"CGT" => 'R',
    b"CGC" => 'R',
    b"CGA" => 'R',
    b"CGG" => 'R',
    //Serine (S)
    b"AGT" => 'S',
    b"AGC" => 'S',
    b"TCT" => 'S',
    b"TCC" => 'S',
    b"TCA" => 'S',
    b"TCG" => 'S',
    //Threonine (T)
    b"ACA" => 'T',
    b"ACG" => 'T',
    b"ACC" => 'T',
    b"ACT" => 'T',
    //Selenocysteine (U)
    b"UGA" => 'U',
    //Valine (V)
    b"GTA" => 'V',
    b"GTG" => 'V',
    b"GTC" => 'V',
    b"GTT" => 'V',
    //Tryptophan (W)
    b"TGG" => 'W',
    //Tyrosine (Y)
    b"TAT" => 'Y',
    b"TAC" => 'Y',
    //Stop Codons
    b"TGA" => '*',
    b"TAA" => '*',
    b"TAG" => '*',
    //Glutamic Acid (E) or glutamine (Q) (Z)
    //X = any of the 13
    //translation stop (*)
    //gap of indeterminate length (-)
};

// Reversed codon map
static REV_CODONS: phf::Map<&'static [u8], char> = phf_map! {
    //Alanine
    b"ACG" => 'A',
    b"GCG" => 'A',
    b"CCG" => 'A',
    b"TCG" => 'A',
    //Aspartic_Acid (D)
    //Asparagine (N)
    //Cysteine
    b"TGT" => 'C',
    b"CGT" => 'C',
    //Aspartic_Acid
    b"CAG" => 'D',
    b"TAG" => 'D',
    //Glutamic_Acid
    b"AAG" => 'E',
    b"GAG" => 'E',
    //Phenylalanine
    b"TTT" => 'F',
    b"CTT" => 'F',
    //Glycine
    b"AGG" => 'G',
    b"GGG" => 'G',
    b"CGG" => 'G',
    b"TGG" => 'G',
    //Histidine
    b"CAC" => 'H',
    b"TAC" => 'H',
    //Isoleucine
    b"TTA" => 'I',
    b"CTA" => 'I',
    b"ATA" => 'I',
    //Leucine (L)
    b"GTT" => 'L',
    b"ATT" => 'L',
    b"ATC" => 'L',
    b"CTC" => 'L',
    b"GTC" => 'L',
    b"TTC" => 'L',
    //Lysine (K)
    b"AAA" => 'K',
    b"GAA" => 'K',
    //Methionine (M)
    b"GTA" => 'M',
    //Asparagine (N)
    b"TAA" => 'N',
    b"CAA" => 'N',
    //Pyrrolysine (O) Special Stop Codon
    b"GAU" => 'O',
    //Proline (P)
    b"ACC" => 'P',
    b"GCC" => 'P',
    b"CCC" => 'P',
    b"TCC" => 'P',
    //Glutamine (Q)
    b"AAC" => 'Q',
    b"GAC" => 'Q',
    //Arginine (R)
    b"AGA" => 'R',
    b"GGA" => 'R',
    b"TGC" => 'R',
    b"CGC" => 'R',
    b"AGC" => 'R',
    b"GGC" => 'R',
    //Serine (S)
    b"TGA" => 'S',
    b"CGA" => 'S',
    b"TCT" => 'S',
    b"CCT" => 'S',
    b"ACT" => 'S',
    b"GCT" => 'S',
    //Threonine (T)
    b"ACA" => 'T',
    b"GCA" => 'T',
    b"CCA" => 'T',
    b"TCA" => 'T',
    //Selenocysteine (U)
    b"AGU" => 'U',
    //Valine (V)
    b"ATG" => 'V',
    b"GTG" => 'V',
    b"CTG" => 'V',
    b"TTG" => 'V',
    //Tryptophan (W)
    b"GGT" => 'W',
    //Tyrosine (Y)
    b"TAT" => 'Y',
    b"CAT" => 'Y',
    //Stop Codons
    b"AGT" => '*',
    b"AAT" => '*',
    b"GAT" => '*',
    //Glutamic Acid (E) or glutamine (Q) (Z)
    //X = any of the 13
    //translation stop (*)
    //gap of indeterminate length (-)
};

pub fn start_parse<Input, Output>(input_path: Input, mut output: Output, n_threads: u32) where
    Input: AsRef<Path>,
    Output: Write
{
    let file_mmap = Mmap::open_path(input_path, Protection::Read).unwrap();
    let bytes: &[u8] = unsafe {
        file_mmap.as_slice() };
//This mmap technique is extremely fast and extremely efficient on large datasets. +1 for memmap
    let mut threadpool = ThreadPool::new(n_threads);
    let (tx, rx) = channel();
    if let IResult::Done(_,o) = fasta_deserialize(bytes) {
        threadpool.scoped(|threadpool| {
            for fasta in o {
                let tx = tx.clone();
                threadpool.execute(move || {
                    let amino_seq: Vec<u8> = fasta.sequence
                        .into_iter()
                        .fold(Vec::new(), |mut acc, item| {
                                acc.extend(item.as_bytes());
                                acc
                        });
                    let mut vec: Vec<FASTA_Complete> = Vec::new();
                    let nomove = FASTA_Complete {
                        window: "> No Move|",
                        id: fasta.id,
                        sequence: no_move(&amino_seq)
                    };
                    let sl1 = FASTA_Complete {
                        window: "> Shift Left One|",
                        id: fasta.id,
                        sequence: nucleotide_shift_left_one(&amino_seq)
                    };
                    let sl2 = FASTA_Complete {
                        window: "> Shift Left Two|",
                        id: fasta.id,
                        sequence: nucleotide_shift_left_two(&amino_seq)
                    };
                    let rnm = FASTA_Complete {
                        window: "> Rev. No Move|",
                        id: fasta.id,
                        sequence: rev_no_move(&amino_seq)
                    };
                    let rsl1 = FASTA_Complete {
                        window: "> Rev. Shift Left One|",
                        id: fasta.id,
                        sequence: rev_nucleotide_shift_left_one(&amino_seq)
                    };
                    let rsl2 = FASTA_Complete {
                        window: "> Rev. Shift Left Two|",
                        id: fasta.id,
                        sequence: rev_nucleotide_shift_left_two(&amino_seq)
                    };
                    vec.push(nomove);
                    vec.push(sl1);
                    vec.push(sl2);
                    vec.push(rnm);
                    vec.push(rsl1);
                    vec.push(rsl2);
                    tx.send(vec).unwrap();
                });
            }
        });

        drop(tx);

        for results in rx {
            for results in results {
                output.write(results.window.as_bytes());
                output.write(results.id.as_bytes());
                output.write(results.sequence.as_bytes());
            }
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
    fn new() -> FASTA_Complete<'a> {
        let sequence = String::new();
        let window = "None";
        let id = "None";
        let FASTA_Complete = FASTA_Complete {
            window: window,
            id: id,
            sequence: sequence
        };
        FASTA_Complete
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
    let mut done = String::new();
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
    let mut done = String::new();
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
    let mut done = String::new();

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
    let mut done = String::new();
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
    let mut done = String::new();
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
