#![feature(plugin)]
#![feature(type_macros)]
#![plugin(phf_macros)]
extern crate phf;

#[macro_use]
extern crate nom;
extern crate memmap;

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
static CODONS: phf::Map<&'static str, u8> = phf_map! {
    //Alanine
    "GCA" => 'A' as u8,
    "GCG" => 'A' as u8,
    "GCC" => 'A' as u8,
    "GCT" => 'A' as u8,
    //Aspartic_Acid (D)
    //Asparagine (N)
    //Cysteine
    "TGT" => 'C' as u8,
    "TGC" => 'C' as u8,
    //Aspartic_Acid
    "GAC" => 'D' as u8,
    "GAT" => 'D' as u8,
    //Glutamic_Acid
    "GAA" => 'E' as u8,
    "GAG" => 'E' as u8,
    //Phenylalanine
    "TTT" => 'F' as u8,
    "TTC" => 'F' as u8,
    //Glycine
    "GGA" => 'G' as u8,
    "GGG" => 'G' as u8,
    "GGC" => 'G' as u8,
    "GGT" => 'G' as u8,
    //Histidine
    "CAC" => 'H' as u8,
    "CAT" => 'H' as u8,
    //Isoleucine
    "ATT" => 'I' as u8,
    "ATC" => 'I' as u8,
    "ATA" => 'I' as u8,
    //Leucine (L)
    "TTG" => 'L' as u8,
    "TTA" => 'L' as u8,
    "CTA" => 'L' as u8,
    "CTC" => 'L' as u8,
    "CTG" => 'L' as u8,
    "CTT" => 'L' as u8,
    //Lysine (K)
    "AAA" => 'K' as u8,
    "AAG" => 'K' as u8,
    //Methionine (M)
    "ATG" => 'M' as u8,
    //Asparagine (N)
    "AAT" => 'N' as u8,
    "AAC" => 'N' as u8,
    //Pyrrolysine (O) Special Stop Codon
    "UAG" => 'O' as u8,
    //Proline (P)
    "CCA" => 'P' as u8,
    "CCG" => 'P' as u8,
    "CCC" => 'P' as u8,
    "CCT" => 'P' as u8,
    //Glutamine (Q)
    "CAA" => 'Q' as u8,
    "CAG" => 'Q' as u8,
    //Arginine (R)
    "AGA" => 'R' as u8,
    "AGG" => 'R' as u8,
    "CGT" => 'R' as u8,
    "CGC" => 'R' as u8,
    "CGA" => 'R' as u8,
    "CGG" => 'R' as u8,
    //Serine (S)
    "AGT" => 'S' as u8,
    "AGC" => 'S' as u8,
    "TCT" => 'S' as u8,
    "TCC" => 'S' as u8,
    "TCA" => 'S' as u8,
    "TCG" => 'S' as u8,
    //Threonine (T)
    "ACA" => 'T' as u8,
    "ACG" => 'T' as u8,
    "ACC" => 'T' as u8,
    "ACT" => 'T' as u8,
    //Selenocysteine (U)
    "UGA" => 'U' as u8,
    //Valine (V)
    "GTA" => 'V' as u8,
    "GTG" => 'V' as u8,
    "GTC" => 'V' as u8,
    "GTT" => 'V' as u8,
    //Tryptophan (W)
    "TGG" => 'W' as u8,
    //Tyrosine (Y)
    "TAT" => 'Y' as u8,
    "TAC" => 'Y' as u8,
    //Stop Codons
    "TGA" => '*' as u8,
    "TAA" => '*' as u8,
    "TAG" => '*' as u8,
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
    let file_mmap = Mmap::open_path("/home/dhc-user/transcriptome_translator/test/test-nucleo.FASTA", Protection::Read).unwrap();
    let bytes: &[u8] = unsafe {
        file_mmap.as_slice() };
    if let IResult::Done( o, parsed) = FASTA_Read_AA(bytes) {
        //This comes from Nom. FASTA_Read_AA is taking the mmap'd file and parsing it into a
        //Vec<FASTA_AA>
        let complete = nucleo_to_amino(parsed);
        for write in complete {
            println!("{:?}", write.id);
            print!("{:?}", write.sequence);
            //Now we write to disk.
        }
        // Is for debugging
    //    print!(" o is, {:?}", str::from_utf8(o) );
    }
}

pub struct FASTA_Complete<'a> {
    id: &'a str,
    sequence: String,
}
pub fn nucleo_to_amino(read: Vec<FASTA_AA>) -> Vec<FASTA_Complete> {
    let mut completed_fastas = Vec::<(FASTA_Complete)>::new();
      for s in &read {
          let mut seq = s.sequence.clone();
          let mut amino_seq: Vec<&str> = seq.drain(..).collect::<Vec<&str>>();
          let mut amino_seq = amino_seq.join("").into_bytes();
//          println!("\n>{:?} | No Shift", id);
          let result = no_move(amino_seq);
          let complete = FASTA_Complete {
              id: s.id,
              sequence: result
          };
          completed_fastas.push(complete);


      }
      /*
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
      */
      completed_fastas
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

pub fn no_move<'a>(mut amino_clone: Vec<u8>) -> String {
    let mut done = Vec::<u8>::new();
    if amino_clone.len() % 3 == 0 {
    //Is it possible to do this without the loop.
        while amino_clone.is_empty() == false {
            let mapped = amino_clone.drain(..3).collect::<Vec<u8>>();
            let mapped = String::from_utf8(mapped);
            for map in mapped {
                let mapped = CODONS.get(&*map);
                match mapped {
                    Some(ref p) => done.push(**p),
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
                        Some(ref p) => done.push(**p),
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
                            Some(ref p) => done.push(**p),
                            None => println!("Done!"),
                        }
                    }
                }
            }
        }
    }
    let done = String::from_utf8(done).unwrap();
    done
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
