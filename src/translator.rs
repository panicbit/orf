use codons::{CODONS, REV_CODONS};
use monster::incubation::{SliceDropFirst, SliceDropLast};

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
