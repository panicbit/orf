use phf::Map;

pub static CODONS: Map<&'static [u8], char> = phf_map! {
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
pub static REV_CODONS: Map<&'static [u8], char> = phf_map! {
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
