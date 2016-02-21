#![feature(plugin)]
#![feature(type_macros)]
#![plugin(phf_macros)]
extern crate phf;

#[macro_use]
extern crate nom;
extern crate memmap;
extern crate radix_trie;

use radix_trie::Trie;
use std::str;
use std::fs::File;
use std::prelude::*;
use std::path::Path;
use nom::{alphanumeric, not_line_ending,line_ending,eof,multispace, GetOutput, newline, eol};
use nom::IResult;
use std::iter::Map;
use std::str::Chars;
use std::fmt::Display;
use memmap::{Mmap, Protection};
use std::iter::IntoIterator;
use std::vec::*;
use std::boxed::Box;

#[derive(Debug)]
pub struct FASTA_AA<'a> {
    pub id: &'a str,
    pub sequence: Vec<&'a str>
}

#[derive(Debug)]
pub struct FASTA_Nucleotide<'a> {
    pub id: &'a str,
    pub sequence: Vec<Nucleotide_Base>
}

#[derive(Debug, Clone)]
pub struct AA_Base {
    pub key: AA_Key
}
#[derive(Debug, Clone)]
pub struct Nucleotide_Base {
    pub key: FASTA_Nucleotide_Key
}
#[derive(Debug,PartialEq,Eq,Clone)]
#[repr(u8)]
pub enum FASTA_Nucleotide_Key {
    A,
    C,
    G,
    T,
    U,
    R,
    Y,
    K,
    M,
    S,
    W,
    B,
    D,
    H,
    V,
    N,
    Hyphen = b'-',
    Asterisk = b'*'
}
#[derive(Debug, Clone)]
#[repr(u8)]
pub enum AA_Key {
    A, //Alanine
    B, //Aspartic Acid (D) or Asparagine (N)
    C, //Cysteine
    D, //Aspartic Acid
    E, //Glutamic Acid
    F, //Phenylalanine
    G, //Glycine
    H, //Histidine
    I, //Isoleucine
    J, //Leucine (L) or Isoleucine (I)
    K, //Lysine
    L, //Leucine
    M, //Methionine
    N, //Asparagine
    O, //Pyrrolysine
    P, //Proline
    Q, //Glutamine
    R, //Arginine
    S, //Serine
    T, //Threonine
    U, //Selenocysteine
    V, //Valine
    W, //Tryptophan
    Y, //Tyrosine
    Z, //Glutamic Acid (E) or Glutamine (Q)
    X,
    Asterisk = b'*',
    Hyphen = b'-',
    Error
}
/*
const Amino_Acid = {
    //Alanine
    'A' => "GCA",
    'A' => "GCG",
    'A' => "GCC",
    'A' => "GCT",
    //Cysteine
    'C' => "TGT",
    'C' => "TGC",
    //Apsartic Acid
    'D' => "GAC",
    'D' => "GAT",
    //Glutamic Acid
    'E' => "GAA",
    'E' => "GAG",
    //Phenylalanine
    'F' => "TTT",
    'F' => "TTC",
    'G' => "GGA",
    'G' => "GGG",
    'G' => "GGC",
    'G' => "GGT",
    //Histidine
    'H' => "CAC",
    'H' => "CAT",
    //Isoleucine
    'I' => "ATT",
    'I' => "ATC",
    'I' => "ATA",
    //Leucine
    'L' => "TTG",
    'L' => "TTA",
    'L' => "CTA",
    'L' => "CTC",
    'L' => "CTG",
    'L' => "CTT",
    //Lysine
    'K' => "AAA",
    'K' => "AAG",
    //Methionine
    'M' => "ATG",
    //Asparagine
    'N' => "AAT",
    'N' => "AAC",
    //Pyrrolysine
    'O' => "UAG",
    //Proline
    'P' => "CCA",
    'P' => "CCG",
    'P' => "CCC",
    'P' => "CCT",
    //Glutamine
    'Q' => "CAA",
    'Q' => "CAG",
    //Arginine
    'R' => "AGA",
    'R' => "AGG",
    //Serine
    'S' => "AGT",
    'S' => "AGC",
    //Threonine
    'T' => "ACA",
    'T' => "ACG",
    'T' => "ACC",
    'T' => "ACT",
    //Selenocysteine
    'U' => "UGA",
    //Valine
    'V' => "GTA",
    'V' => "GTG",
    'V' => "GTC",
    'V' => "GTT",
    //Tryptophan
    'W' => "TGG",
    'Y' => "TAT",
    'Y' => "TAC"
};
*/
static CODONS: phf::Map<&'static str, &'static str> = phf_map! {
    //Alanine
    "GCA" => "A",
    "GCG" => "A",
    "GCC" => "A",
    "GCT" => "A",
    //Aspartic_Acid (D)
    //Asparagine (N)
    //Cysteine
    "TGT" => "C",
    "TGC" => "C",
    //Aspartic_Acid
    "GAC" => "D",
    "GAT" => "D",
    //Glutamic_Acid
    "GAA" => "E",
    "GAG" => "E",
    //Phenylalanine
    "TTT" => "F",
    "TTC" => "F",
    //Glycine
    "GGA" => "G",
    "GGG" => "G",
    "GGC" => "G",
    "GGT" => "G",
    //Histidine
    "CAC" => "H",
    "CAT" => "H",
    //Isoleucine
    "ATT" => "I",
    "ATC" => "I",
    "ATA" => "I",
    //Leucine (L)
    "TTG" => "L",
    "TTA" => "L",
    "CTA" => "L",
    "CTC" => "L",
    "CTG" => "L",
    "CTT" => "L",
    //Lysine (K)
    "AAA" => "K",
    "AAG" => "K",
    //Methionine (M)
    "ATG" => "M",
    //Asparagine (N)
    "AAT" => "N",
    "AAC" => "N",
    //Pyrrolysine (O) Special Stop Codon
    "UAG" => "O",
    //Proline (P)
    "CCA" => "P",
    "CCG" => "P",
    "CCC" => "P",
    "CCT" => "P",
    //Glutamine (Q)
    "CAA" => "Q",
    "CAG" => "Q",
    //Arginine (R)
    "AGA" => "R",
    "AGG" => "R",
    "CGT" => "R",
    "CGC" => "R",
    "CGA" => "R",
    "CGG" => "R",
    //Serine (S)
    "AGT" => "S",
    "AGC" => "S",
    "TCT" => "S",
    "TCC" => "S",
    "TCA" => "S",
    "TCG" => "S",
    //Threonine (T)
    "ACA" => "T",
    "ACG" => "T",
    "ACC" => "T",
    "ACT" => "T",
    //Selenocysteine (U)
    "UGA" => "U",
    //Valine (V)
    "GTA" => "V",
    "GTG" => "V",
    "GTC" => "V",
    "GTT" => "V",
    //Tryptophan (W)
    "TGG" => "W",
    //Tyrosine (Y)
    "TAT" => "Y",
    "TAC" => "Y",
    //Stop Codons
    "TGA" => "*",
    "TAA" => "*",
    "TAG" => "*",
    //Glutamic Acid (E) or glutamine (Q) (Z)
    //X = any of the 13
    //translation stop (*)
    //gap of indeterminate length (-)
};

pub enum Human_Readable_AA {
    Adenine,
    Cytosine,
    Guanine,
    Thymine,
    Uracil,
    Purine,
    Pyrimidines,
}


pub fn start_parse() {
    let file_mmap = Mmap::open_path("/Path/To/File", Protection::Read).unwrap();
    let bytes: &[u8] = unsafe {
        file_mmap.as_slice() };
    if let IResult::Done( o, parsed) = FASTA_Read_AA(bytes) {
        //This comes from Nom. FASTA_Read_AA is taking the mmap'd file and parsing it into a
        //Vec<FASTA_AA>
        nucleo_to_amino(parsed);


        // Build_Constructs(parsed);
        // Is for debugging
    //    print!(" o is, {:?}", str::from_utf8(o) );
    }
}

pub fn nucleo_to_amino(read: Vec<FASTA_AA>) {
      for s in &read {
          let mut seq = s.sequence.clone();
          let mut id = s.id;
          let mut amino_seq: Vec<&str> = seq.drain(..).collect::<Vec<&str>>();
          let mut amino_seq = amino_seq.join("").into_bytes();
          println!("\n>{:?} | No Shift", id);
          no_move(amino_seq);
      }
      for s in &read {
          let mut seq = s.sequence.clone();
          let mut amino_seq: Vec<&str> = seq.drain(..).collect::<Vec<&str>>();
          let mut amino_seq = amino_seq.join("").into_bytes();
          println!("\n>{:?} | Shift Left One", s.id);
          nucleotide_shift_left_one(amino_seq);
      }
      for s in &read {
          let mut seq = s.sequence.clone();
          let mut amino_seq: Vec<&str> = seq.drain(..).collect::<Vec<&str>>();
          let mut amino_seq = amino_seq.join("").into_bytes();
          println!("\n>{:?} | Shift Left Two", s.id);
          nucleotide_shift_left_two(amino_seq);
      }
      for s in &read {
          let mut seq = s.sequence.clone();
          let mut amino_seq: Vec<&str> = seq.drain(..).collect::<Vec<&str>>();
          let mut amino_seq = amino_seq.join("");
          let mut amino_seq = amino_seq.into_bytes();
          println!("\n>{:?} | Rev No Shift", s.id);
          rev_no_move(amino_seq);
      }
      for s in &read {
          let mut seq = s.sequence.clone();
          let mut amino_seq: Vec<&str> = seq.drain(..).collect::<Vec<&str>>();
          let mut amino_seq = amino_seq.join("").into_bytes();

          println!("\n>{:?} | Rev Shift Left One", s.id);
          rev_nucleotide_shift_left_one(amino_seq);
      }
      for s in &read {
          let mut seq = s.sequence.clone();
          let mut amino_seq: Vec<&str> = seq.drain(..).collect::<Vec<&str>>();
          let mut amino_seq = amino_seq.join("").into_bytes();
          println!("\n>{:?} | Rev Shift Left Two", s.id);
          rev_nucleotide_shift_left_two(amino_seq);
      }
}

pub fn rev_nucleotide_shift_left_two(mut amino_seq: Vec<u8>) {
    amino_seq.reverse();
    amino_seq.remove(0);
    amino_seq.remove(0);
    if amino_seq.len() % 3 == 0 {
    //Pop other two nucleotides off the back to keep product
    //a remainder of three.
        while amino_seq.is_empty() == false {
            let mapped = amino_seq.drain(..3).collect::<Vec<u8>>();
            let mapped = String::from_utf8(mapped);
            for map in mapped {
                let mapped = CODONS.get(&*map);
                match mapped {
                    Some(ref p) => print!("{}", p),
                    None => println!("Done!"),
                }
            }
        }
    } else {
        amino_seq.pop();
        if amino_seq.len() % 3 == 0 {
            while amino_seq.is_empty() == false {
                let mapped = amino_seq.drain(..3).collect::<Vec<u8>>();
                let mapped = String::from_utf8(mapped);
                for map in mapped {
                    let mapped = CODONS.get(&*map);
                    match mapped {
                        Some(ref p) => print!("{}", p),
                        None => println!("Done!"),
                    }
                }
            }
        } else {
            amino_seq.pop();
            if amino_seq.len() % 3 == 0 {
                while amino_seq.is_empty() == false {
                    let mapped = amino_seq.drain(..3).collect::<Vec<u8>>();
                    let mapped = String::from_utf8(mapped);
                    for map in mapped {
                        let mapped = CODONS.get(&*map);
                        match mapped {
                            Some(ref p) => print!("{}", p),
                            None => println!("Done!"),
                        }
                    }
                }
            }
        }
    }
}

pub fn rev_nucleotide_shift_left_one(mut amino_clone: Vec<u8>) {
    amino_clone.reverse();
    amino_clone.remove(0);
    if amino_clone.len() %3 == 0 {
        //Pop other two nucleotides off the back to keep product
        //a remainder of three.
        while amino_clone.is_empty() == false {
            let mapped = amino_clone.drain(..3).collect::<Vec<u8>>();
            let mapped = String::from_utf8(mapped);
            for map in mapped {
                let mapped = CODONS.get(&*map);
                match mapped {
                    Some(ref p) => print!("{}", p),
                    None => println!("Done!"),
                }
            }
        }
    } else {
        amino_clone.pop();
        if amino_clone.len() % 3 == 0 {
            while amino_clone.is_empty() == false {
                let mapped = amino_clone.drain(..3).collect::<Vec<u8>>();
                let mapped = String::from_utf8(mapped);
                for map in mapped {
                    let mapped = CODONS.get(&*map);
                    match mapped {
                        Some(ref p) => print!("{}", p),
                        None => println!("Done!"),
                    }
                }
            }
        } else {
            amino_clone.pop();
            if amino_clone.len() % 3 == 0 {
                while amino_clone.is_empty() == false {
                    let mapped = amino_clone.drain(..3).collect::<Vec<u8>>();
                    let mapped = String::from_utf8(mapped);
                    for map in mapped {
                        let mapped = CODONS.get(&*map);
                        match mapped {
                            Some(ref p) => print!("{}", p),
                            None => println!("Done!"),
                        }
                    }
                }
            }
        }
    }
}

pub fn rev_no_move(mut amino_clone: Vec<u8>) {
    amino_clone.reverse();
    if amino_clone.len() % 3 == 0 {
    //Is it possible to do this without the loop.
        while amino_clone.is_empty() == false {
            let mapped = amino_clone.drain(..3).collect::<Vec<u8>>();
            let mapped = String::from_utf8(mapped);
            for map in mapped {
                let mapped = CODONS.get(&*map);
                match mapped {
                    Some(ref p) => print!("{}", p),
                    None => println!("Done!"),
                }
            }
        }
    } else {
        amino_clone.pop();
        if amino_clone.len() % 3 == 0 {
            while amino_clone.is_empty() == false {
                let mapped = amino_clone.drain(..3).collect::<Vec<u8>>();
                let mapped = String::from_utf8(mapped);
                for map in mapped {
                    let mapped = CODONS.get(&*map);
                    match mapped {
                        Some(ref p) => print!("{}", p),
                        None => println!("Done!"),
                    }
                }
            }
        } else {
            amino_clone.pop();
            if amino_clone.len() % 3 == 0 {
                while amino_clone.is_empty() == false {
                    let mapped = amino_clone.drain(..3).collect::<Vec<u8>>();
                    let mapped = String::from_utf8(mapped);
                    for map in mapped {
                        let mapped = CODONS.get(&*map);
                        match mapped {
                            Some(ref p) => print!("{}", p),
                            None => println!("Done!"),
                        }
                    }
                }
            }
        }
    }
}


pub fn nucleotide_shift_left_two(mut amino_seq: Vec<u8>) {
    amino_seq.remove(0);
    amino_seq.remove(0);
    if amino_seq.len() % 3 == 0 {
    //Pop other two nucleotides off the back to keep product
    //a remainder of three.
        while amino_seq.is_empty() == false {
            let mapped = amino_seq.drain(..3).collect::<Vec<u8>>();
            let mapped = String::from_utf8(mapped);
            for map in mapped {
                let mapped = CODONS.get(&*map);
                match mapped {
                    Some(ref p) => print!("{}", p),
                    None => println!("Done!"),
                }
            }
        }
    } else {
        amino_seq.pop();
        if amino_seq.len() % 3 == 0 {
            while amino_seq.is_empty() == false {
                let mapped = amino_seq.drain(..3).collect::<Vec<u8>>();
                let mapped = String::from_utf8(mapped);
                for map in mapped {
                    let mapped = CODONS.get(&*map);
                    match mapped {
                        Some(ref p) => print!("{}", p),
                        None => println!("Done!"),
                    }
                }
            }
        } else {
            amino_seq.pop();
            if amino_seq.len() % 3 == 0 {
                while amino_seq.is_empty() == false {
                    let mapped = amino_seq.drain(..3).collect::<Vec<u8>>();
                    let mapped = String::from_utf8(mapped);
                    for map in mapped {
                        let mapped = CODONS.get(&*map);
                        match mapped {
                            Some(ref p) => print!("{}", p),
                            None => println!("Done!"),
                        }
                    }
                }
            }
        }
    }
}

pub fn nucleotide_shift_left_one(mut amino_clone: Vec<u8>) {
    amino_clone.remove(0);
    if amino_clone.len() %3 == 0 {
        //Pop other two nucleotides off the back to keep product
        //a remainder of three.
        while amino_clone.is_empty() == false {
            let mapped = amino_clone.drain(..3).collect::<Vec<u8>>();
            let mapped = String::from_utf8(mapped);
            for map in mapped {
                let mapped = CODONS.get(&*map);
                match mapped {
                    Some(ref p) => print!("{}", p),
                    None => println!("Done!"),
                }
            }
        }
    } else {
        amino_clone.pop();
        if amino_clone.len() % 3 == 0 {
            while amino_clone.is_empty() == false {
                let mapped = amino_clone.drain(..3).collect::<Vec<u8>>();
                let mapped = String::from_utf8(mapped);
                for map in mapped {
                    let mapped = CODONS.get(&*map);
                    match mapped {
                        Some(ref p) => print!("{}", p),
                        None => println!("Done!"),
                    }
                }
            }
        } else {
            amino_clone.pop();
            if amino_clone.len() % 3 == 0 {
                while amino_clone.is_empty() == false {
                    let mapped = amino_clone.drain(..3).collect::<Vec<u8>>();
                    let mapped = String::from_utf8(mapped);
                    for map in mapped {
                        let mapped = CODONS.get(&*map);
                        match mapped {
                            Some(ref p) => print!("{}", p),
                            None => println!("Done!"),
                        }
                    }
                }
            }
        }
    }
}

pub fn no_move(mut amino_clone: Vec<u8>) {
    if amino_clone.len() % 3 == 0 {
    //Is it possible to do this without the loop.
        while amino_clone.is_empty() == false {
            let mapped = amino_clone.drain(..3).collect::<Vec<u8>>();
            let mapped = String::from_utf8(mapped);
            for map in mapped {
                let mapped = CODONS.get(&*map);
                match mapped {
                    Some(ref p) => print!("{}", p),
                    None => println!("Done!"),
                }
            }
        }
    } else {
        amino_clone.pop();
        if amino_clone.len() % 3 == 0 {
            while amino_clone.is_empty() == false {
                let mapped = amino_clone.drain(..3).collect::<Vec<u8>>();
                let mapped = String::from_utf8(mapped);
                for map in mapped {
                    let mapped = CODONS.get(&*map);
                    match mapped {
                        Some(ref p) => print!("{}", p),
                        None => println!("Done!"),
                    }
                }
            }
        } else {
            amino_clone.pop();
            if amino_clone.len() % 3 == 0 {
                while amino_clone.is_empty() == false {
                    let mapped = amino_clone.drain(..3).collect::<Vec<u8>>();
                    let mapped = String::from_utf8(mapped);
                    for map in mapped {
                        let mapped = CODONS.get(&*map);
                        match mapped {
                            Some(ref p) => print!("{}", p),
                            None => println!("Done!"),
                        }
                    }
                }
            }
        }
    }
}

pub fn FASTA_Read_AA(input:&[u8]) -> IResult<&[u8], Vec<FASTA_AA>>  {
    many0!(input,
      chain!(
        tag!(">") ~
        id: map_res!(not_line_ending, str::from_utf8) ~ line_ending ~
        sequence: many0!(terminated!(map_res!( is_not!(">\n"), str::from_utf8), tag!("\n"))),
        ||{
            FASTA_AA {
                id: id,
                sequence: sequence
            }
        }
      )
   )
}
fn main() {
    start_parse();
}
           
