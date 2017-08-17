#!/usr/bin/env python
from __future__ import print_function  
import argparse
import os.path
import sys
import subprocess

"""
$ launchme.py -f reads.fastq 
"""

inputfile = '' 
ctgname= ''
min_length=0
max_length=0
fq2fa=False
contigs={}
outfile=''
filename=''
out_type=''



def main():

   global inputfile  
   global ctgname
   global min_length
   global max_length
   global fq2fa
   global outfile

   # location of handler:
   myseq_handler = os.getcwd()

   usage="\n %(prog)s -f full_path_to_file [options]" 
   parser = argparse.ArgumentParser(usage=usage)
 

   required = parser.add_argument_group('Required Arguments')
   required.add_argument("-f", dest="filename",
                         help="full path to a fasta/fastq file") #metavar="FILE", 
                         #required=True)
   parser.add_argument("-v", 
                  action="store_true", dest="verbose", default=False,
                  help="don't print status messages to stdout")
 

   writing = parser.add_argument_group('writing arguments')
   writing.add_argument("--fq2fa", dest="fq2fa",  action="store_true", default=False,
                     help="Write fasta from fastq ---NA---")
   writing.add_argument("--write", dest="write",  action="store_true", default=False,
                     help="Write selected contigs to file ---NA---")

   # Write new file input
   # select and print only 1 ctg ### maybe later add list from file? ###
   select  = parser.add_argument_group('select arguments')
   select.add_argument("--ctg", dest="ctgname",
                     help="Print only this contig/read ---NA---")
   select.add_argument("--min", dest="min_length", type=int,
                     help="Print only contigs/reads longer than min_length  ---NA--- ")
   select.add_argument("--max", dest="max_length", type=int,
                     help="Print only contigs/reads shorter than max_length---NA---")

   action  = parser.add_argument_group('calc arguments')
   action.add_argument("--stat", action="store_true", help="Calculate bases,ctg_num,longest,mean,n50,n_n50")

 
   if len(sys.argv)==1:
       parser.print_help()
       sys.exit(1)

   args = parser.parse_args()

   if not os.path.exists(args.filename): 
        print("Sorry, file ", args.filename, "does not exists")
        sys.exit(2)
   inputfile=args.filename

   ## Optional arguments ##
   printout=''
   preout=''  # prefix for file name
   if args.fq2fa or args.write:
       printout='\n Writing to file'
       if args.fq2fa:
           printout+="\n   Change format: fasta from fastq";
           fq2fa=args.fq2fa
           preout+="fq2fa_"
   
   #### Selection arguments ####
   if any(f for f in (args.ctgname, args.min_length, max_length)): 
       printout+='\n Contig Selection:'
   if args.ctgname is not None:
       printout+="\n   Write only contig/read named: " + args.ctgname;
       ctgname=args.ctgname
       preout+=ctgname+"_"
   if args.min_length is not None:
       printout+="\n   Write only contig/reads longer than "+str(args.min_length);
       min_length=args.min_length
       preout+="min"+"_"+min_length+"_"
   if args.max_length is not None:
       printout+="\n   Write only contig/reads shorter than " + str(args.max_length);
       max_length=args.max_length
       preout+="max"+"_"+max_length+"_"


   #print(args)
   if preout=='' and args.write:
       print("\n ******* Error: re-writing the same file? No actions defined! ******* \n")
       parser.print_help()
       sys.exit(2)
       
   if args.stat:
       printout+='\n Calculating:'
       printout+="\n   stats  ";
       
   if printout == '':
       print("\n *** Nothing to be done. ***\n")
       sys.exit(2)

   if args.verbose: print(printout,'\n')


   if args.stat:
      exe=myseq_handler+"/src/n50/n50"
      result = subprocess.check_output([exe, inputfile])
      print(result)

   '''
   from Bio import SeqIO
   from statistics import mean
   read_ok=read();
   if read_ok > 0:
       sys.exit(2)

   if args.stat:
       stats_ok=stats()
       
   # write to file
   if args.write or args.fq2fa:
       outfile=preout+outfile
       write_ok=write()
       if write_ok > 0:
           sys.exit(2)
   '''        




#### My functions ####

def filetype_(filename):
    global outfile
    global out_type


    typefasta=('fasta','fsa','fa')
    typefastq=('fastq','fq')

    if len(filename.split('.gz')) > 1:  # if more than 1 ext (.fa.gz f.i.)
        filext = '.'.join(filename.split('.')[-2:])
    else:
        filext= '.'.join(filename.split('.')[-1:])  # fasta or fastq
    basename=filename.split('.'+filext)[0:]
 
    ## check fastq first, otherwise can confuse "fa"stq for fa in fasta exts
    if any(ext in filext for ext in typefastq):
        file_type='fastq'
        if fq2fa: 
            outfile=basename[0]+'.fasta'
            out_type='fasta'
        else: 
            outfile=basename[0]+'.fastq'
            out_type='fastq'
    elif any(ext in filext for ext in typefasta):
        file_type='fasta'
        outfile=basename[0]+'.fasta'
        out_type='fasta'
    else:
        print("\n  Error: I don't recognize the input file extension,\n   please select a fasta or fastq file only!")
        return 1
    
    return file_type

def read():
    global contigs
    global filename

    filename = os.path.basename(inputfile)
    file_type = filetype_(filename)
    
    # file extension not recognized:
    if file_type == 1:
       return 1
    elif fq2fa and file_type=='fasta':
       print("\n ******* Error: input file is already a fasta!  ******* \n")
       return 2

    contigs = (r for r in SeqIO.parse(inputfile, file_type))
 
    return 0


def stats():
    
    lengths=[]
    for r in contigs: 
        lengths.append(len(r.seq))

    lengths.sort(reverse=True)
    
    ctgs=len(lengths)
    bases=sum(lengths)
    mean_length=int(round(mean(lengths)))
    longest=max(lengths)
    
    ii=0
    done=0
    t50=0
    while done == 0: 
        t50+=lengths[ii]
        if t50 > bases*0.5: done=1
        ii+=1
    n_n50=ii
    n50=lengths[n_n50-1]

    print("Bases=", bases, "Contigs=", ctgs,
          "mean=", mean_length, "longest=", longest, 
          "n50=", n50, "n=", n_n50)
    return 0


def write():
    output=open(outfile,"w")
   
    nn=0 
    for r in contigs:  
        pri=1
        if min_length > 0 and len(r.seq) < min_length: pri=0
        if max_length > 0 and len(r.seq) > max_length: pri=0
        if ctgname is not None and ctgname not in r.id: 
            pri=0
 
        if pri: 
            SeqIO.write(r, output, out_type)
            nn+=1
        
    output.close()

    if nn > 0:
        print(" Done:", nn, "contigs written to", outfile, "\n")
    else:
        print(" Sorry, could not find any contigs with your specifications\n")
        
    return 0


if __name__ == "__main__":
    main()

