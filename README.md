# orf - open reading frame
[![Travis](https://img.shields.io/travis/panicbit/orf.svg?style=flat-square)][travis]
[![Coveralls](https://img.shields.io/coveralls/panicbit/orf.svg?style=flat-square)][coveralls]
[![Crates.io Version](https://img.shields.io/crates/v/orf.svg?style=flat-square)][crates.io]
[![Crates.io Downloads](https://img.shields.io/crates/dv/orf.svg?style=flat-square)][crates.io]
[![Crates.io License](https://img.shields.io/crates/l/orf.svg?style=flat-square)][crates.io]

[crates.io]: https://crates.io/crates/orf
[travis]: https://travis-ci.org/panicbit/orf
[coveralls]: https://coveralls.io/github/panicbit/orf

Translates nucleotide sequences to Amino Acid sequences. (RNA and DNA, shifts reading frames and reverses)

Takes a FASTA file of nucleotide sequences and outputs another FASTA file with Amino Acids.

Deserializing is done with nom.

Threadpools are utilized to allow scaling with hardware, after nom has finished parsing tasks are dispatched as threads become available and then the file is written into a fasta format text file.

WARNING: You must edit main.rs to the /home/user/transcriptome_translator/test/test-nucleo.FASTA for this to work properly, it will panic because of the unwrap for loading the file to memmap.


# Compile Time Instructions:
  Included in the repo is test/test-nucleo.FASTA, these are test FASTA encoded sequence data than can be used for your testing purposes.  To compile your data find where the file is loaded in fn start_parse() and modify the Path to suit your needs, currently only works on linux.
  
  The path must be the full path e.x "/home/user/transcriptome_translator/test/test-nucleo.FASTA"
  
# Caveats(for which there were many):
  FASTA encodings must decode to utf-8 compliant format.
  
##TODO:

* [ ] Use fold_many_0! from nom to allow parsing deserialized structures as they are completed instead of weaiting for a vector of types.
* [ ] Remove warnings by fixing camel_casing and other assorted warnings.
* [ ] Remove the FASTA parser from main.rs and move it to it's own folder/repo.
* [ ] Complete FASTA spec to include support for ';'.
* [ ] Create a parser to read FASTA in reverse to eliminate the need for vec.reverse(), this should be done with another thread.

#Contributors

Last but certainly not least:

[panicbit](https://github.com/panicbit)
