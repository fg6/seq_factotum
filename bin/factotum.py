#!/usr/bin/env python
from __future__ import print_function  
import argparse
import os.path
import sys
import subprocess
from shutil import which

inputfile = '' 
seqname= ''
min_length=0
max_length=0
fq2fa=False
contigs={}
outfile=''
filename=''
out_type=''
nlist=0


def main():

   global inputfile  
   global seqname
   global min_length
   global max_length
   global fq2fa
   global outfile
   global write
   global nlist

   write=0

   # location of handler:
   factotum_bin="/".join(which("factotum.py").split('/')[:-1])

   usage="\n %(prog)s -f full_path_to_file [options]" 

   parser = argparse.ArgumentParser(usage=usage)
   parser.add_argument("-q", action="store_true", dest="quite", default=False,
                       help="don't print status messages to stdout")

   required = parser.add_argument_group('Required Arguments')
   required.add_argument("-f", dest="filename",
                         help="full path to a fasta/fastq file") 

   # select and print only 1 ctg ### later add list from file? ###
   # add option write each contg in different file, but question if too many ctgs #
   # jolly: add list chromosomes up to 20, larger and smaller
   action  = parser.add_argument_group('Action arguments:')
   action.add_argument("--list", dest="nlist", type=int, 
                     help="List this many seqs (from longest)")
   action.add_argument("--stats", action="store_true", help="Print stats: bases, seq_num, longest, mean, n50, n_n50")
   action.add_argument("--break", dest="scaffbreak", action="store_true",
                     help="Break scaffolds @>3Ns, print ctgs to file and print stats")
   action.add_argument("--fq2fa", dest="fq2fa",  action="store_true", default=False,
                     help="Write fasta from fastq")
   action.add_argument("--seq", dest="seqname",
                     help="Write in file only this Seq")
   action.add_argument("--min", dest="min_length", type=int,
                     help="Write in file only Seqs longer than min_length")
   action.add_argument("--max", dest="max_length", type=int,
                     help="Write in file only Seqs shorter than max_length")
     
   if len(sys.argv)==1:
       parser.print_help()
       sys.exit(1)

   args = parser.parse_args()

   if not os.path.exists(args.filename): 
        print("Sorry, file ", args.filename, "does not exists")
        sys.exit(2)


   ##################################
   ###########  Arguments ###########
   ##################################
   printout=''
   preout=''  # prefix for output file name
   if args.fq2fa:
      write=1
      printout+="\n Factotum: Changing format: fasta from fastq";
      fq2fa=args.fq2fa
      preout+="fq2fa_"
     
   if args.seqname is not None:
       printout+="\n Factotum: Writing only seq named: " + args.seqname;
       seqname=args.seqname 
       preout+=seqname+"_"
       write=1
   if args.min_length is not None:
       printout+="\n Factotum: Writing only seq longer than "+str(args.min_length);
       min_length=args.min_length
       preout+="min" + "_" + str(min_length) + "_"
       write=1
   if args.max_length is not None:
       printout+="\n Factotum: Writing only seq shorter than " + str(args.max_length);
       max_length=args.max_length
       preout+="max" +"_" + str(max_length) + "_"
       write=1
   if args.scaffbreak:
      printout+="\n Factotum: Breaking scaffolds @ 3Ns";
      preout+="ctgs"+"_"
      write=1
   if args.stats:
      printout+='\n Factotum: Calculating stats'
   if args.nlist:
      if args.nlist == 1: printout += ' Factotum: Listing the longest seqs'
      else:  printout+=' Factotum: Listing the ' + str(args.nlist) + ' longest seqs'
      printout+=' and the shortest one'
      nlist = args.nlist
         
   ##################################
   ######## output file name ########
   ##################################
   inputfile=args.filename
   filename = os.path.basename(inputfile)
   file_type = filetype_(filename)
   outfile=preout+outfile
   if write and filename == outfile:
       print("\n ******* Error: re-writing the same file! Exiting now ******* \n")
       sys.exit(2)
   elif write:
      printout+='\n Output file: ' + outfile
      

   ##################################
   ######### some  printout #########
   ##################################
   if printout == '':
       print("\n *** Nothing to be done. ***\n")
       sys.exit(2)
   elif not args.quite: print(printout) #,'\n')



   ###################################
   ############# Actions #############
   ###################################
   if args.scaffbreak:
      exe=factotum_bin+"/ctgs_from_scaff"
      print("\n",subprocess.check_output([exe, inputfile, outfile]).decode('utf-8').strip())
      
   if args.stats:
      exe=factotum_bin+"/n50"
      print("\n",subprocess.check_output([exe, inputfile]).decode('utf-8').strip())
      

   if args.seqname or args.min_length or args.max_length or args.fq2fa or args.nlist:
      exe=factotum_bin+"/jolly"
      print("\n",subprocess.check_output([exe, inputfile, outfile, out_type, str(min_length), 
                                          str(max_length), seqname, str(nlist)]).decode('utf-8').strip())
   print(' ')



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
        if fq2fa:
           print("Error: file is already a fasta!")
           sys.exit(2)
    else:
        print("\n  Error: I don't recognize the input file extension,\n   please select a fasta or fastq file only!")
        return 1
    

    return file_type







#### My functions: only used with internal calcs: obsolete ####

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

