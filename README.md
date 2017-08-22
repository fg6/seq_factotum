# seq_factotum
A Python script to handle genomic sequences: calculate sequences stats, select a sequence, get contigs from scaffolds, print fasta from fastq

## Usage: 
    factotum.py -f full_path_to_file [options]
  
    Optional Arguments:
      -h, --help        show this help message and exit
      -q                don't print status messages to stdout

    Required Arguments:
      -f FILENAME       full path to a fasta/fastq file
  
    Action arguments::

      --list NLIST      List this many seqs (from longest)
      --stats           Print stats: bases, seq_num, longest, mean, n50, n_n50
      --break           Break scaffolds @>3Ns, print ctgs to file and print stats
      --fq2fa           Write fasta from fastq
      --seq SEQNAME     Write in file only this Seq
      --min MIN_LENGTH  Write in file only Seqs longer than min_length
      --max MAX_LENGTH  Write in file only Seqs shorter than max_length
  
## Requirements
Python 3 
Tested with Python 3.6.1, gcc 4.9.2

## Instructions
Download repository and install utilities/compile tools: 

	$ git clone https://github.com/fg6/seq_factotum.git
	$ cd seq_factotum
	$ ./install.sh
	
To run, add the factotum bin location to your PATH, then run:

	$ myfactotum=`pwd`; PATH=$myfactotum/bin/:$PATH   
	$ factotum.py -f file.fasta --stats
