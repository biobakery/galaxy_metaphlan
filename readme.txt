============================================
MetaPhlAn: Metagenomic Phylogenetic Analysis
============================================

Nicola Segata: nsegata@hsph.harvard.edu

MetaPhlAn is a computational tool for profiling the composition of microbial communities from metagenomic shotgun sequencing data. MetaPhlAn relies on unique clade-specific marker genes identified from 3,0000 reference genomes, allowing:
- orders of magnitude speedups compared to existing methods;
- unambiguous taxonomic assignments;
- accurate estimation of organismal relative abundance;
- species-level resolution for bacterial and archaeal organisms.

If you use this software, please cite our paper:
"Fast and accurate metagenomic profiling of microbial community composition using unique clade-specific marker genes" 
Nicola Segata, Levi Waldron, Annalisa Ballarini, Vagheesh Narasimhan, Olivier Jousson, Curtis Huttenhower. 
Nature Methods, in press

-------------
PREREQUISITES
-------------
* MetaPhlAn requires python 2.7 or higher with argparse, tempfile and numpy libraries installed (apart for numpy they are usually installed together with the python distribution).

If you provide the output of BLASTN or BowTie2 as input, there are no additional prerequisite.

If you would like to use the integrated BLASTN in MetaPhlAn, you need to have:
* blastn from the NCBI BLAST+ package version 2.2.25 or higher (blastn needs to be in the system path)

If you would like to use the BowTie2 integrated in MetaPhlAn, you need to have:
* BowTie2 version 2.0.0 or higher and perl (bowtie2 needs to be in the system path with execute _and_ read permission)


-----------
BASIC USAGE
-----------

* Profiling a metagenome from raw reads using Blast (requires BLAST 2.2.25+ and the blast marker DB provided with MetaPhlAn):
metaphlan.py metagenome.fasta --blastdb blastdb/mpa

* You can save a lot of computational time if you perform the blasting using BowTie2 instead of Blast
  (requires BowTie2 in the system path with execution and read permissions, Perl installed, and the BowTie2 marker DB provided with MetaPhlAn):
metaphlan.py metagenome.fasta --bowtie2db bowtie2db/mpa

* you can take advantage of multiple CPUs and you can save the blast output for re-running MetaPhlAn extremely quickly:
metaphlan.py metagenome.fasta --blastdb blastdb/mpa --nproc 5 --blastout metagenome.outfmt6.txt
  and the same for BowTie2:
metaphlan.py metagenome.fasta --bowtie2db bowtie2db/mpa --nproc 5 --bowtie2out metagenome.bt2out.txt

* if you already blasted your metagenome against the marker DB (using MetaPhlAn or blastn/BowTie2 alone) you can obtain the results in few seconds:                                                                                                                            
metaphlan.py metagenome.outfmt6.txt                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                               
* the blast/BowTie output can also be provided with a pipe:                                                                                                                                                                                                                    
metaphlan.py < metagenome.outfmt6.txt > profiling_output.txt                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                               
* you can also set advanced options for the BowTie2 step selecting the preset option among  'sensitive','very-sensitive','sensitive-local','very-sensitive-local' (valid for metagenome as input only)                                                                         
metaphlan.py --bt2_ps very-sensitive-local metagenome.fasta                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                               
* for for Blast the main configurable option is the threshold on the evalue (default 1e-5, we strongly recommend not lowering it too much):                                                                                                                                    
metaphlan.py --evalue 1e-7 < metagenome.outfmt6.txt > profiling_output.txt                                                                                                                                                                                                     
                                                                                                                                                                                                                                                                               
* if you suspect that the metagenome contains unknown clades, you may obtain more accurare result with a more sensible blast search lowering the blast word_size (the blasting will be slower):                                                                                
metaphlan.py metagenome.fna --word_size 12 > profiling_output.txt 

-------------------------
FULL COMMAND LINE OPTIONS
-------------------------

usage: metaphlan.py [-h] [-v] [-t ANALYSIS TYPE] [--tax_lev TAXONOMIC_LEVEL]
                    [--blastdb METAPHLAN_BLAST_DB]
                    [--bowtie2db METAPHLAN_BOWTIE2_DB] [--evalue]
                    [--word_size] [--bt2_ps BowTie2 presets] [--tmp_dir]
                    [--min_cu_len]
                    [--input_type {automatic,multifasta,blastout}] [--stat_q]
                    [--blastout FILE_NAME] [--bowtie2out FILE_NAME]
                    [--nproc N]
                    [INPUT_FILE] [OUTPUT_FILE]


DESCRIPTION
 MetaPhlAn version 1.5.1 (13 April 2012): METAgenomic PHyLogenetic ANalysis for taxonomic classification of metagenomic reads.

AUTHORS: Nicola Segata (nsegata@hsph.harvard.edu)

positional arguments:
  INPUT_FILE            the input file can be:
                        * a multi-fasta file containing metagenomic reads
                        OR
                        * a NCBI BLAST output file (-outfmt 6 format) of the metagenome against the MetaPhlAn database. 
                        OR
                        * a BowTie2 output file of the metagenome generated by a previous MetaPhlAn run 
                        The software will recognize the format automatically.
                        If the input file is missing, the script assumes that a BLAST or BowTie2 output is passed to the stdin
  OUTPUT_FILE           the tab-separated output file of the predicted taxon relative abundances 
                        [stdout if not present]

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Prints the current MetaPhlAn version and exit
  -t ANALYSIS TYPE      Type of analysis to perform: 
                         * rel_ab: profiling a metagenomes in terms of relative abundances
                         * reads_map: mapping from reads to clades (only reads hitting a marker)
                         * clade_profiles: normalized marker counts for clades with at least a non-null marker
                        [default 'rel_ab']
  --tax_lev TAXONOMIC_LEVEL
                        The taxonomic level for the relative abundance output:
                        'a' : all taxonomic levels
                        'k' : kingdoms (Bacteria and Archaea) only
                        'p' : phyla only
                        'c' : classes only
                        'o' : orders only
                        'f' : families only
                        'g' : genera only
                        's' : species only
                        [default 'a']
  --blastdb METAPHLAN_BLAST_DB
                        The blast database file of the MetaPhlAn database 
  --bowtie2db METAPHLAN_BOWTIE2_DB
                        The BowTie2 database file of the MetaPhlAn database 
  --evalue              evalue threshold for the blasting
                        [default 1e-6]
  --word_size           word_size value for the blasting
                        [default NCBI BlastN default]
  --bt2_ps BowTie2 presets
                        presets options for BowTie2 (applied only when a multifasta file is provided)
                        The choices enabled in MetaPhlAn are:
                         * sensitive
                         * very-sensitive
                         * sensitive-local
                         * very-sensitive-local
                        [default very-sensitive-local]
  --tmp_dir             the folder used to store temporary files 
                        [default is the OS dependent tmp dir]
  --min_cu_len          minimum total nucleotide length for the markers in a clade for
                        estimating the abundance without considering sub-clade abundances
                        [default 10000]
  --input_type {automatic,multifasta,blastout}
                        set wheter the input is the multifasta file of metagenomic reads or 
                        the blast output (outfmt 6 format) of the reads against the MetaPhlAn db.
                        [default 'automatic', i.e. the script will try to guess the input format]
  --stat_q              Quantile value for the robust average
                        [default 0.1]
  --blastout FILE_NAME  The file for saving the output of the blasting (in outfmt 6 format)
  --bowtie2out FILE_NAME
                        The file for saving the output of BowTie2
  --nproc N             The number of CPUs to use for parallelizing the blasting
                        [default 1, i.e. no parallelism]



---------------
VERSION HISTORY
---------------

======== Version 1.5.2
Bug fix:
- division by zero bug fixed: the bug was occurring when no hits were detected for any of the reads

======== Version 1.5.1
New options:
- option -o added for specifying the output file with a non-positional command line argument

======== Version 1.5
New Features:
- added support for BowTie2 in addition to Blast. Bowtie2 runs at 10,000 reads-per-second (20,000 with 'sensitive' instead of 'very-sensitive-local') with accuracy performances at least comparable to blastn
-- '--bowtie2out' and '--bowtie2db' are the BowTie2 options corresponding to the '--blastout' and '--blastdb' blast options
-- support for the BowTie2 metaphlan output which is is simply a two-column file of reads/marker assignements. If called outside MetaPhlAn, BowTie2 should be called as
   bowtie2 --quiet --sam-no-hd --sam-no-sq --very-sensitive-local --local -x bowtie2db/mpa -f -U metagenome.fna | cut -f 1,3 | grep -P -v "\*$"
-- the BowTie2 MetaPhlAn database in in bowtie2db, is called mpa, and consists of 6 files
-- four preset BowTie2 options are available (--bt2_ps option): --very-sensitive, --very-sensitive, --sensitive-local, --very-sensitive-local
- help message 'metaphlan.py -h' expanded
- '--version' option added

Bug fixes:
- Fixed the error message for missing blastdb file (thanks to Joshua Reyes)

======== Version 1.1
New features and changes in MetaPhlAn version 1.1 with respect to version 1.0:
- no need to provide the MetaPhlAn data file anymore (thus --lib_dir has been removed)
- a robust average is introduced to more accurately estimate taxonomic relative abundances and dramatically lower false positives (--stat_q specifies the quantile for the truncated mean)
- the "word_size" option for blastn has been added in order to handle more accurately metagenomes without closely related reference genomes
- a new analysis type (-t clade_profiles) has been introduced in order to estimate normalized read counts for markers in each clade
- when using MetaPhlAn on blasting outputs, the input file can be provided using a pipe
- a more comprehensive help message has been added (see metaphlan.py -h)
- the code has been re-factored and made more robust for future developments
- a slightly modified blastdb has been developed. As a consequence, the database (both multifasta and blastdb) of version 1.0 is not compatible with version 1.1



