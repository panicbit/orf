#![feature(plugin)]
#![feature(type_macros)]
#![plugin(phf_macros)]
extern crate phf;

#[macro_use]
extern crate nom;
extern crate memmap;
extern crate scoped_threadpool;

use std::str;
use std::io::prelude::*;
use nom::{not_line_ending,line_ending};
use nom::IResult;
use memmap::{Mmap, Protection};
use std::vec::*;
use std::fs::OpenOptions;
use std::sync::mpsc::{Sender, Receiver};
use std::sync::mpsc;
use std::thread;
use std::mem;
use scoped_threadpool::Pool as ThreadPool;
use std::sync::mpsc::channel;
use std::rc::Rc;

#[derive(Debug)]
pub struct FASTA<'a> {
    pub id: &'a str,
    pub sequence: Vec<&'a str>
}

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

pub fn start_parse() {
    let file_mmap = Mmap::open_path("/home/dhc-user/transcriptome_translator/test/test-nucleo.FASTA", Protection::Read).unwrap();
    let bytes: &[u8] = unsafe {
        file_mmap.as_slice() };
//This mmap technique is extremely fast and extremely efficient on large datasets. +1 for memmap
    let mut file = OpenOptions::new().create(true).read(false).write(true).open("./results.txt").unwrap();
    let mut threadpool = ThreadPool::new(4);
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
                        sequence: no_move(amino_seq.clone())
                    };
                    let sl1 = FASTA_Complete {
                        window: "> Shift Left One|",
                        id: fasta.id,
                        sequence: nucleotide_shift_left_one(amino_seq.clone())
                    };
                    let sl2 = FASTA_Complete {
                        window: "> Shift Left Two|",
                        id: fasta.id,
                        sequence: nucleotide_shift_left_two(amino_seq.clone())
                    };
                    let rnm = FASTA_Complete {
                        window: "> Rev. No Move|",
                        id: fasta.id,
                        sequence: rev_no_move(amino_seq.clone())
                    };
                    let rsl1 = FASTA_Complete {
                        window: "> Rev. Shift Left One|",
                        id: fasta.id,
                        sequence: rev_nucleotide_shift_left_one(amino_seq.clone())
                    };
                    let rsl2 = FASTA_Complete {
                        window: "> Rev. Shift Left Two|",
                        id: fasta.id,
                        sequence: rev_nucleotide_shift_left_two(amino_seq.clone())
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
                file.write(results.window.as_bytes());
                file.write(results.id.as_bytes());
                file.write(results.sequence.as_bytes());
            }
        }
        file.sync_all();

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

pub fn rev_nucleotide_shift_left_two(amino_clone: Vec<u8>) -> String {
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
    let mut done = Vec::<u8>::new();
    let mut amino_clone = amino_clone;
    done.push(b'\n');
    amino_clone.reverse();
    amino_clone.remove(0);
    amino_clone.remove(0);
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
    done.push(b'\n');
    let done = String::from_utf8(done).unwrap();
    done
}

pub fn rev_nucleotide_shift_left_one(amino_clone: Vec<u8>) -> String {
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
    let mut done = Vec::<u8>::new();
    let mut amino_clone = amino_clone;
    done.push(b'\n');
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
    done.push(b'\n');
    let done = String::from_utf8(done).unwrap();
    done
}

pub fn rev_no_move(amino_clone: Vec<u8>) -> String {
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
    let mut done = Vec::<u8>::new();
    let mut amino_clone = amino_clone;
    done.push(b'\n');
    amino_clone.reverse();
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
    done.push(b'\n');
    let done = String::from_utf8(done).unwrap();
    done
}

pub fn nucleotide_shift_left_two(amino_clone: Vec<u8>) -> String {
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
    let mut done = Vec::<u8>::new();
    let mut amino_clone = amino_clone;
    amino_clone.remove(0);
    amino_clone.remove(0);
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
    done.push(b'\n');
    let done = String::from_utf8(done).unwrap();
    done
}

pub fn nucleotide_shift_left_one(amino_clone: Vec<u8>) -> String {
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
    let mut done = Vec::<u8>::new();
    let mut amino_clone = amino_clone;
    done.push(b'\n');
    amino_clone.remove(0);
    if amino_clone.len() %3 == 0 {
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
    done.push(b'\n');
    let done = String::from_utf8(done).unwrap();
    done
}

pub fn no_move<'a>(amino_clone: Vec<u8>) -> String {
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
    let mut done = Vec::<u8>::new();
    let mut amino_clone = amino_clone;
    done.push(b'\n');
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
    done.push(b'\n');
    let done = String::from_utf8(done).unwrap();
    done
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
fn main() {
    start_parse();
}
