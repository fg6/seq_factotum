# seq_factotum
A factotum Python script to handle genomic sequences: calculate sequences stats (n50,..), select a sequence, get contigs from scaffolds, print fasta from fastq, select sequences whose name does/doesn't match a string, select only sequences longer/shorter than a value, etc..

## Usage: 
    factotum.py -f full_path_to_file [options]
  
    Optional Arguments:
      -h, --help        show this help message and exit
      -q                don't print status messages to stdout

    Required Arguments:
      -f FILENAME       full path to a fasta[.gz]/fastq[.gz] file
  
    Action arguments::

      --list NLIST        List this many seqs (from longest)
      --stats             Print stats: bases, seq_num, longest, mean, n50, n_n50
      --break             Break scaffolds @>3Ns, print ctgs to file and print stats
      --fq2fa             Write fasta from fastq
      --seq SEQNAME       Write in file only this Seq
      --keep KEEPNAME     Write in file only Seqs matching keepname
      --rm RMNAME         Write in file only Seqs not matching rmname
      --min MIN_LENGTH    Write in file only Seqs longer than min_length
      --max MAX_LENGTH    Write in file only Seqs shorter than max_length
      --avoid AVOID.list  Write in file Seqs not matching list in this file
      --nocomments        Write fasta/fastq file with no comments in read name line
 
  
## Requirements
Python 3 
Tested with Python 3.6.1, gcc 4.9.2

## Instructions
Download repository and install utilities/compile tools: 

	$ git clone https://github.com/fg6/seq_factotum.git
	$ cd seq_factotum
	$ ./install.sh
	
To run, add the factotum 'bin' folder location to your PATH, then run:

	$ myfactotum=`pwd`; PATH=$myfactotum/bin/:$PATH   
	$ factotum.py -f file.fasta --stats

	
## External packages
The seq_factotum downloads and installs the gzstream library to handle gzip input files (https://www.cs.unc.edu/Research/compgeom/gzstream/)
