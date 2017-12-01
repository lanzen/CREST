**CREST - Classification Resources for Environmnetal Sequence Tags**, is a collection of software and databases for taxonomic classification of environmental marker genes from sequencing-based community profiling studies (also known as "meta-" + "-genomics", " -transcriptomics", "-barcoding"; "taxonomic-" or "phylogenetic-" profiling)). The program LCAClassifier is used used for classification of sequences aligned to the reference databases provided.

If you use CREST or the LCAClassifier in your research, please cite: 

Lanzén A , Jørgensen SL, Huson D, Gorfer M, Grindhaug SH, Jonassen I, Øvreås L, Urich T (2012) CREST - Classification Resources for Environmental Sequence Tags, *PLoS ONE* **7**:e49334

## Classification databases supported

SilvaMod was derived by manual curation of the SILVA nr SSU Ref. It supports SSU sequences from bacteria and archaea (16S) as well as eukaryotes (18S), with a high level of manual curation and defined environmental clades.

Greengenes is an alternative reference database for classification of prokaryotic 16S, curated and maintained by The Greengenes Database Consortium. Release supported: May 2013

Unite seeks to maintain the cleanest possible copy of all public fungal ITS sequences. Third-party sequence annotation is supported, and everyone with the expertise to improve the sequence data and their annotation/metadata are welcome to participate in the annotation effort. Release supported: v7.2, 2017-10-10

To update your database files, new versions can be downloaded separately (wihtout the need to install a new version of CREST), from http://services.cbu.uib.no/supplementary/crest

## Installation

LCAClassifier uses pairwise alignments to a reference database of marker genes (such as SSU rRNA) as input. Alignment files must be in XML format and produced by the NCBI blastall suit. For best performance, we recommend the Megablast algorithm. Megablast is both implemented as a separate program in NCBI's legacy Blastall suite, that can be downloaded from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/; or as part of the _blastn_ program of the newer BLAST+ implementation, available at ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/. The later implementation is slightly faster.

LCAClassifier also requires python (v.2.7 or higher), included as part of recent MacOSX and Linux distributions.
Installation procedure. It also requires Python-Dev.

Download the latest stable distribution and expand it using tar:

`tar -xvzf LCAClassifier3.0.6.tar.gz`

Or download the development version:

`git clone https://github.com/lanzen/CREST.git`

Then go to the LCAClassifier directory and run the script install.sh:

`cd LCAClassifier`
`./install.sh`

The script will take a while (up to 30 minutes) as several python packages as well as the CREST reference databases are downloaded automatially. It will also produce lots of output including warnings that you can safely ignore.

The executable file classifiy should then have been generated in the directory `bin` and is ready to use! To easily access it from anywhere in your file system, either add this directory to your PATH environment variable or create a symbolic link to it in a directory already in your path, for example $HOME/bin:

## How to use the LCAClassifier

### Classification algorithm

Taxonomical classification starts with alignment to a reference sequence database (such as SilvaMod or Greengenes) using the NCBI Blast+ or legacy blastall implementation of Megablast. The LCAClassifier requires that output is saved in XML format. The classification is then carried out based on a subset of the best matching alignments using the Lowest Common Ancestor (LCA) of this subset. Briefly, the subset includes sequences that score within x% of the “bit-score” of the best alignment, providing the best score is above a minimum value. Default values for the minimum bit-score is 155 and for the LCA range (x) 2%. Based on cross-validation testing using the non-redundant SilvaMod database, this results in relatively few false positives for most datasets. However, the LCA range can be increased up to about 10% (using the -r option), to increase accuracy with short reads and for datasets with many novel sequences.

In addition to LCA classification, a minimum similarity filter is used, based on a set of taxon-specific requirements, by default depending on their taxonomic rank. By default, a sequence must be aligned with at least 99% nucleotide similarity to the best reference sequence in order to be classified to the species rank. For the genus, family, order, class and phylum ranks the respective default cut-offs are 97%, 95%, 90%, 85% and 80%. These cutoffs can be changed manually by editing the .map file of the respective reference database, or deactivated using option -f. This filter ensures that classification is made to the taxon of the lowest allowed rank, effectively re-assigning sequences to parent taxa until allowed.

### Preparing sequences

For amplicon sequences, we strongly recommend noise reduction using e.g AmpliconNoise for 454 Pyrosequencing data; vsearch, UPARSE, DADA2, SWARM or similar for Illumina MiSeq and other platforms; as well as chimera removal (using e.g. UCHIME) prior to submission.

###Preparing an OTU table

For amplicon sequencing experiments with many replicates or similar samples (>~10), unique noise-reduced sequences may be further clustered using a similarity threshold (often 97% although larger thresholds are probably preferable) into Operational Taxonomic Units (OTUs), prior to classification. Alternatively, the unique sequence variants may be used as OTUs. The LCAClassifier can then use the distribution data of OTUs across datasets, in order to determine weighted relative abundance for each taxon across datasets, as well are OTU richness and a Chao-estimate of total richness.

An OTU table defines the distribution of an OTU across datasets and can be in delimeted text format or in BIOM-format (the later is default in QIIME and the former in many other applications). Make sure to divide datasets into combinations that "make sense", i.e. that you want to analyse together, before OTU clustering. For example, dividing data from a sequencing run of two separate experiements, or merging together data from two sequencing runs belonging to the same experiment. Also make sure sequence / OTU names in the OTU table correspond to those used in the FASTA file with representative sequences, to be aligned.

###Using individual sequence datasets

In some cases, OTU clustering may not be possible, e.g. when classifying data from shotgun transcriptomic or metagenomic sequencing experiments. In this case, one BLAST output file must be prepared per dataset, instead of submitting an OTU table.


###Aligning sequences using Megablast

First download and install the NCBI blast+ from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/. Then use the blastn program with the Silvamod or Greengenes reference databases that are installed automatically with the LCAClassifier. Other marker genes can also be used, but require that the construction of a custom database (see below).

Silvamod and Greengenes are installed in the directory where LCAClassifier was installed, under parts/flatdb.
Example

To align the sequences of fasta file dataset1.fa to the Silvamod database, providing LCAClassifier was installed in the home directory, use, with blastn (BLAST+ sutie):

`blastn -task megablast -query dataset1.fa -db ~/LCAClassifier/parts/flatdb/silvamod/silvamod.fasta -num_alignments 100 -outfmt 5 -out dataset1_silvamod.xml`

We also recommend that you use -num_threads n to enable multi-threading and speed up the alignment. You can also use gzip to compress your blast output in order to conserve disk space since CREST supports gzipped input files.

...or for the older legacy blast:
`megablast -i dataset1.fa -d ~/LCAClassifier/parts/flatdb/silvamod/silvamod.fasta -b100 -v100 -m7 -o dataset1_silvamod.xml`

Alignments will be written to the file dataset1_silvamod.xml.

### Classification

The executable for the LCAClassifier is called classify and uses a very simple command line interface:

`classify [options] alignment-file(s).xml`

To classify the alignment file dataset1_silvamod.xml using default parameters, simply use:

`classify dataset1_silvamod.xml`

Silvamod is the default reference database. In case another database, e.g. UNITE was used, the -d option must be used to specify the database name:

`classify -d unite dataset1_unite.xml`

Let’s say that dataset1.fa contains representative OTU sequneces. Normally, we would be more interested in the taxonomic composition across each individual sample included in that dataset. The information of how OTUs are distributed will then be contained in an OTU-table, for example called dataset1_OTUs.csv. To incorporate this information, use:

`classify -t dataset1_OTUs.csv dataset1_silvamod.xml`

If the OTU table is instead in BIOM-format, we use the -b option:

`classify -b dataset1_OTUs.biom dataset1_silvamod.xml`

In both cases, a new OTU-table, in the same format and with the same name as given, will be written to the output directory, with taxonomic annotations added (for text-format as an extra column and for BIOM as observation metadata using the key taxonomy. We can also write annotated fasta-file, if the representative sequences are provided (option -i):

`classify -i dataset1.fa -b dataset1_OTUs.biom dataset1_silvamod.xml`

To specify output directory name, use option -o:

`classify -o dataset1_CREST_Out -i dataset1.fa -b dataset1_OTUs.biom dataset1_silvamod.xml`

If several datasets, not clustered into common OTUs, were aligned to the reference database, classification can be carried out for alignment XML file at a time, or simultaneously for several files. The advantage of the later option is that rows will always be inserted in the composition output files for all taxa that are present in at least one of the datasets, even if they are not present in all. For example, to classify dataset1_silvamod.xml.gz and dataset2_silvamod.xml.gz simultaneously, use:

`classify dataset1_silvamod.xml.gz dataset2_silvamod.xml.gz`

Note that gzip-compressed alignment files can be used - practical for avoiding to use a lot of disk space!

### Output

The output of LCAClassifier is written to the directory specified by option -o (by default CREST_Results. It produces at least two output files for each sample dataset classified. If an OTU-table is provided, the output files are named after sample names used in it. If not, they are named using the XML alignment files from which they were derived.
Sample-specific files

For example, for a (compressed) alignment file named Dataset.xml.gz and no specified OTU-table, LCAClassifier will write the following files:

**Relaitve_Abundaces.tsv** provides relative abundance data across datasets normalised to the total number of assigned reads in each dataset. Option -m can be used to filter this output from the least abundant taxa. This is recommended particularly when datasets of different sequencing depths (total number of reads), will be compared.

**Dataset_Composition.tsv** contains the taxonomical composition of the dataset in tab-separated text format, designed for reading in a spreadsheet editor like Open Office (or Excel). Composition is given separately for each rank from domain to species, starting with a meta-level giving the total number of reads, classified and unclassfied reads. For each taxon, the total number of reads, its relative abundance, number of unique sequences and a chao estimate of total unique sequences is given. For each rank, the total number of reads assigned to that rank or better is also given.

**Dataset_Tree.txt** shows the taxonomical composition in a simple space-indented tree (in plaintext format), each taxon annotated with the number of reads.

**Dataset_Assignments.fasta** written when using option -a, this file gives all taxonomically annotated sequences in FASTA format. In the FASTA header, the original sequence name is first given, then the predicted taxonomic assignment from the root of the taxonomic tree to the best rank at which classification was possible (separated by semicolons). If the Minimum similarity filter prevented assignment at a higher level for this sequence, indicating it to represent a novel taxa, the last taxon is prefixed with the word "Unknown".

**Dataset_Assignments.tsv** written by using option -p, this file lists all sequence names and their assignments in a simpe semicolon-separated format. If the Minimum similarity filter prevented assignment at a higher level for this sequence, indicating it to represent a novel taxa, the last taxon is prefixed with the word "Unknown".
Common / dataset-specific output files

**All_Assignments.tsv** provides counts of the number of assignments at each taxonomic rank for all datasets in a tab-separated format. Note that only assignments to the taxon node itself are counted, not to child taxa at lower ranks. For each taxon, the full taxonomic path from root to the taxon itself is also given. This file is more suitable for parsing than Dataset_Composition.txt.

**All_Cumulative.tsv** provides cumulative counts for the number of assignments at each taxonomic rank for all datasets in a tab-separated format. As opposed to All_Assignments.tsv, assignments to child taxa (lower ranks) are also counted. This file is more suitable for parsing than Dataset_Composition.txt.

**Richness.tsv** provides number of unique OTUs for each taxon

**<Sequence-input-file-name> (if given) +"_Assigned.fasta"** sequences with assignments added in FASTA header. Written with option -i

Other options

-h: shows a help message listing all options

-r: LCA bitscore range in percent (given as an integer >0; default = 2)

-s: Minimum bit-score for classification (integer >0; default = 155)

-n: Normalise relative abundances to classified sequences only (ignoring unclassified)

-f: De-activate the minimum similarity filter

-i fasta-file: By default sequences in the taxonomically annotatated FASTA-output (_Assignments.fasta) are taken directly from the alignment XML file. This means that un-aligned sequence stretches or completely un-aligned sequences are ignored. This option instead uses the sequences from the supplied FASTA-file

-q qual-file: This option works like -i, but instead accepts a file with Phred quality scores in fasta.qual format and outputs an annotated .qual file.

-v: Use Verbose mode. Outputs information about interesting assignments to standard out.

-c: Selects which configuration file to use (by default, the file parts/etc/lcaclassifier.conf in the LCAClassifier’s installation directory is used)
Constructing a custom reference database

Starting from a reference alignment of a phylogenetic marker gene annotated with taxonomical information, the script nds2CREST.py can create a custom reference database for use with the LCAClassifier. We recommend using the program ARB to create and annotate such an alignment from a set of sequences.

nds2CREST.py requires at least two different input files:

1. Sequences of the reference marker genes in the alignment, in FASTA format (sequences.fasta in the example below). It is important that the sequences are cropped so that stretches of unaligned sequence is not included. Such sequences may bias the alignments during classification. In ARB, this can be done by using a positional filter when exporting sequences.

2. A tab-separated text file (in ARB language "NDS file") containing the taxonomic annotations of each aligned sequence (taxa.nds in the example below). This file must contain three columns separated by tabs. The first column is the accession number corresponding to the sequence file. The second contains the taxonomic structure with ranks separated by slash, semicolon or underscores (e.g. Domain;Phylum;Class). The third column contains the species or strain name or a unique description of the sequence. This file can be exported from ARB by first setting up NDS display correctly in the menu Tree > NDS (Node Display Setup). For SILVA, check the leaf box for "acc", then for group for "tax_slv" and for leaf for "full name". Then write the file using File > Export > Export fields using NDS. Make sure to tick the checkbox marked "Use TABs for columns"! Example: "A16379 Bacteria/Proteobacteria/Gammaproteobacteria_1/Pasteurellales_Pasteurellaceae/Haemophilus Haemophilus ducreyi"

`nds2CREST.py -o database-name -i sequences.fasta taxa.nds`

A file containing changes to the taxonomic tree for implementation during parsing can also be supplied (option -c) as well as a file containing rank information, in the same format as used by NCBI taxonomy (option -r.). For helo. use option -h.

Database parsing may result in a list of warnings being written to standard.error if for example a sequence included sequences.fasta was not in the taxa.nds file or vice versa, or if duplicate names of taxa appeared in different topological context. In the former case, the sequence will be ignored. In the later, a parent suffix will be added to the sequence, for example "Cryptococcus", which is the name of an insect genus as well as a fungal one, will become "Cryptococcus (Eriococcidae)". The same happens if the same name is used for two different ranks, e.g. "Bacteria;Actinobacteria;Actinobacteria" will become "Bacteria;Actinobacteria (phylum);Actinobacteria (class)".

After successful completion, nds2CREST will create three new files:

- **database-name.tre**: the taxonomy of the reference database in http://en.wikipedia.org/wiki/Newick_format[Newick format], each taxon represented by a unique numeral identifier
- **database-name.map**:  a mapping file, mapping the taxa identifiers to their full names and taxonomic rank
- **database-name.fasta**: the reference database in FASTA-format (simply the same as sequences.fasta but with un-mapped sequences removed).

These files are needed by the LCAClassifier in order to make taxonomic assignments and the .tre and .map files for Silvamod can be found in /LCAClassifier/LCADir/parts/flatdb/silvamod/silvamod.map or silvamod.tre if LCAClassifier was installed in the $HOME directory (). To add your own custom database with the somewhat unimaginative name "database-name" as used in the examples above, create a new directory named "database-name" under the directory flatdb and copy these two files there:

`mkdir ~/LCAClassifier/LCADir/parts/flatdb/database-name`
`cp database-name.* ~/LCAClassifier/LCADir/parts/flatdb/database-name`

You also need to format the FASTA-file for Megablast (or BLAST) using the command formatdb (installed with the blastall suite):

`cd ~/LCAClassifier/LCADir/parts/flatdb/database-name`
`formatdb -i database-name.fasta -pF`

The final change you need to do before using your own reference database is to tell the configuration file of the LCAClassifier where to find it. Edit the file ~/LCAClassifier/parts/etc/lcaclassifier.conf and add a new line. For example:

`database-name = /Home/user/LCAClassifier/parts/flatdb/database-name`

Note that this will be overwritten if you update the LCAClassifier or re-build it. To make the change permanent, also add it in the file ~/LCAClassifier/etc/lcaclassifier.conf.in

As you notice, you can keep your reference database anywhere you want in the file system, not necessarily in the LCAClassifier installation directory, as long as you make the correct changes in lcaclassifier.conf

To use your new reference database:

`megablast -i env.fa -d ~/LCAClassifier/parts/flatdb/database-name/database-name.fasta -b100 -v100 -m7 -o env_custom.xml`
`classifiy -d database-name env_custom.xml`

