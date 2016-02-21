# transcriptome_translator
Translates nucleotide sequences to Amino Acid sequences. (RNA and DNA, shifts reading frames and reverses)


Currently you must use vim to capture stdout.  

binary > file.txt will not work if the file you expect to be generating is greater than 20MB.


To use Vim to capture stdout: 1) Open Vim 2) enter command '! multirust run nightly-2016-02-01 cargo run'
Vim will appear to do nothing but wait a few minutes and once your binary has finished executing you will see that vim has your output.  The first 5 lines will be cargo output so make sure you trim those if you don't need them.

# Compile Time Instructions:
  Included in the repo is test/test-nucleo.FASTA, these are test FASTA encoded sequence data than can be used for your testing purposes.  To compile your data find where the file is loaded in fn start_parse() and modify the Path to suit your needs, currently only works on linux.
  
# Caveats(for which there are many):
  FASTA encodings must decode to utf-8 compliant format.
  


##TODO:
Remove warnings by fixing camel_casing and other assorted warnings.
Save file the right way.
Remove the FASTA parser from main.rs and move it to it's own folder/repo.
Complete FASTA spec to include support for ';'.
