# transcriptome_translator
Translates nucleotide sequences to Amino Acid sequences. (RNA and DNA, shifts reading frames and reverses)


WARNING: You must edit main.rs to the /home/user/transcriptome_translator/test/test-nucleo.FASTA for this to work properly, it will panic because of the unwrap for loading the file to memmap.


# Compile Time Instructions:
  Included in the repo is test/test-nucleo.FASTA, these are test FASTA encoded sequence data than can be used for your testing purposes.  To compile your data find where the file is loaded in fn start_parse() and modify the Path to suit your needs, currently only works on linux.
  
  The path must be the full path e.x "/home/user/transcriptome_translator/test/test-nucleo.FASTA"
  
# Caveats(for which there were many):
  FASTA encodings must decode to utf-8 compliant format.
  
##TODO:

Remove warnings by fixing camel_casing and other assorted warnings.

Save file the right way.

Remove the FASTA parser from main.rs and move it to it's own folder/repo.

Complete FASTA spec to include support for ';'.

Make serializing parralel if possible.
