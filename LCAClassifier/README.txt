** LCAClassifier version 2.0 (March 2014) **

LCAClassifier is a program for taxonomic classification of environmental sequence datasets, 
using the  Lowest Common Ancestor (LCA) algorithm to a custom reference database and taxonomy. 
LCAClassifier is a part of  CREST - Classification Resources for Environmnetal Sequence Tags.

NEW FEATURES IN THIS VERSION (2.0):

- Output is now written to a specific output directory, specified using option -o (by default "CREST_Results") 

- OTU-tables in delimeted text or BIOM format can now be imported together with alignments of 
  representative sequences for each OTU (curated unique sequences or clusters thereof) 
  to a CREST reference database. Composition for individual is determined from abundances in this table.
  Further, a taxonomically annotated OTU table or BIOM file is written to the output directory.
  When an OTU table is used, abundance data need not be included as annotation of sequence names.  
  (If not, it still does however, using the default assumption that each Blast-file in XML format 
  represents one dataset with the same name.) 
  
- Some new output files are written by default to the output directory:
  . Relaitve_Abundaces.tsv: relative abundance data normalised to the total number of assigned reads
  . Richness.tsv: number of unique OTUs for each taxon
  . <Sequence-input-file-name> (if given) +"_Assigned.fasta": sequences with assignments added in FASTA header    

- A minimum relative abundance can be specified (option -m) for inclusion in output Relative_Abundance.tsv
  (taxa below this cutoff in all datasets are excluded).
  
- Database alignments in gzipped (XML) format can now be used as input

- Input from the blastn program part of NCBI's (relatively) new BLAST+ package can be used as an alternative to
  legacy Megablast. Recommeded settings are "-task megablast -outfmt 5" (the later, ie XML format, being mandatory)
      

REQUIREMENTS:

LCAClassifier uses pairwise alignments to a reference database of marker genes (such as SSU rRNA) as input. 
Alignment files must be in XML format and produced by the NCBI blastall suit. For best performance, we 
recommend Megablast, as implemented by blastn in the the NCBI BLAST+ package or Megablast from the legacy
blastall program. Platfrom specific pre-combiled binaries of the NCBI Blastall can be downloaded from 
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/. LCAClassifier also requires  python, which 
is included as part of MacOSX and most Linux distributions. 

INSTALLATION:

    python bootstrap.py
    bin/buildout

If you use the system python, boostrapping may require sudo (sudo python bootstrap.py), unless you have installed the python package setuptools, or are logged in as root. The executable file classifiy should now have been generated in the directory bin and is ready to use! To easily access it from anywhere in your file system, either add this directory to your PATH environment variable or create a symbolic link to it in a directory already in your path, for example $HOME/bin:

   cd ~/bin
   ln -s ~/LCAClassifier/bin/classify .

(replace "~" with the path to your where you installed LCAClassifier in case this was not your home-directory) 

USAGE:

Type classify -h for documentation of the command line interface

A full user guide can be found at https://www.ii.uib.no/trac/LCAClassifier/wiki/Userguide
