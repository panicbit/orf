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
           
