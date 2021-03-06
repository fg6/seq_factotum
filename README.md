# seq_factotum
A factotum Python script to handle genomic sequences:

	1. calculate sequences stats (bases, seq_num, longest, mean, n50,..)
	2. list the longest N sequences and their lengths
	3. select and print out a single sequence
	4. print to file contigs from scaffolds and summary stats
	5. print a fasta from fastq
	6. select sequences whose name does/doesn't match a string
	7. select only sequences longer/shorter than a value
	8. select and printout all sequences not in a file.list
	9. print fasta/q after removing the comments from the sequence name lines


## Usage: 
    factotum.py -f full_path_to_file [options]
  
    Optional Arguments:
      -h, --help        show this help message and exit
      -q                don't print status messages to stdout

    Required Arguments:
      -f FILENAME       full path to a fasta[.gz]/fastq[.gz] file
  
    Action arguments::

      --list NLIST        List information for these many seqs (from longest) plus the shortest sequence
      --stats             Print stats: bases, seq_num, longest, mean, n50, n_n50
      --break             Break scaffolds @>3Ns, print ctgs to file and print stats for contigs
      --fq2fa             Write fasta from fastq
      --seq SEQNAME       Write in file only this Seq
      --keep KEEPNAME     Write in file only Seqs matching keepname
      --rm RMNAME         Write in file only Seqs not matching rmname
      --min MIN_LENGTH    Write in file only Seqs longer than min_length (can be mixed with --max)
      --max MAX_LENGTH    Write in file only Seqs shorter than max_length (can be mixed with --min)
      --avoid AVOID.list  Write in file Seqs not in list of sequences in file AVOID.list
      --nocomments        Write fasta/fastq file with no comments in read name line
 
  
## Requirements
Python 3 

Tested with Python 3.6.1, gcc 4.9.2

## Instructions
Download repository and install utilities/compile tools: 

	$ git clone https://github.com/fg6/seq_factotum.git
	$ cd seq_factotum
	$ ./install.sh
	
To run, add the factotum 'bin' folder location to your PATH before running:

	$ PATH=/full/path/to/seq_factotum/bin/:$PATH   

Example: Check stats for a fasta:

	$ factotum.py -f file.fasta --stats

Example: Print to file only sequences longer than 10Kbp and shorter than 100kbp:

	$ factotum.py -f file.fasta --min 10000 --max 100000


	
## External packages
The seq_factotum downloads and installs the gzstream library to handle gzip input files (https://www.cs.unc.edu/Research/compgeom/gzstream/)
