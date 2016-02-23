# transcriptome_translator
Translates nucleotide sequences to Amino Acid sequences. (RNA and DNA, shifts reading frames and reverses)


WARNING: You must edit main.rs to the /home/user/transcriptome_translator/test/test-nucleo.FASTA for this to work properly, it will panic because of the unwrap for loading the file to memmap.

WARNING: Currently no_move is the only function enabled.

Requires multirust and nightly-2016-02-01 (For phf_macros)

Requires Vim.

Currently you must use vim to capture stdout.  

binary > results.txt will not work if the file you expect to be generating is greater than 20MB.

To use Vim to capture stdout: 1) Open Vim 2) enter command '! multirust run nightly-2016-02-01 cargo run'
Vim will appear to do nothing but wait a few minutes and once your binary has finished executing you will see that vim has your output.  The first 5 lines will be cargo output so make sure you trim those if you don't need them.

# Compile Time Instructions:
  Included in the repo is test/test-nucleo.FASTA, these are test FASTA encoded sequence data than can be used for your testing purposes.  To compile your data find where the file is loaded in fn start_parse() and modify the Path to suit your needs, currently only works on linux.
  
  The path must be the full path e.x "/home/user/transcriptome_translator/test/test-nucleo.FASTA"
  
# Caveats(for which there are many):
  FASTA encodings must decode to utf-8 compliant format.
  
  Whole thing is fragile and relies on print!("{}"); macros for vim to capture output.
  This works but it's not nice.  You WILL get a FASTA compliant file this way.
  
  The current strategy is to use vim to capture stdout and then rerun re-use the nom FASTA parser to search for sequences.
  
##TODO:

Remove warnings by fixing camel_casing and other assorted warnings.

Save file the right way.

Remove the FASTA parser from main.rs and move it to it's own folder/repo.

Complete FASTA spec to include support for ';'.

Add support for Amino Acid Sequence -> Nucleotide sequences.

Make serializing parralel if possible.
