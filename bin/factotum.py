#!/usr/bin/env python
from __future__ import print_function  
import argparse
import os.path
import sys
import subprocess
from shutil import which

inputfile = '' 
seqname= ''
keepname=''
rmname=''
min_length=0
max_length=0
fq2fa=False
contigs={}
outfile=''
filename=''
out_type=''
avoid=''
remove_comments=''
nlist=0


def main():

   global inputfile  
   global seqname
   global keepname
   global rmname
   global min_length
   global max_length
   global fq2fa
   global outfile
   global write
   global nlist
   global avoid
   global remove_comments
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
   # add reverse complement
   # add gc
   # add length distribution plot
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
   action.add_argument("--keep", dest="keepname",
                     help="Write in file only Seqs matching keepname")
   action.add_argument("--rm", dest="rmname",
                     help="Write in file only Seqs not matching rmname")
   action.add_argument("--min", dest="min_length", type=int,
                     help="Write in file only Seqs longer than min_length")
   action.add_argument("--max", dest="max_length", type=int,
                     help="Write in file only Seqs shorter than max_length")
   action.add_argument("--avoid", dest="avoid",
                    help="Write in file Seqs not matching list in this file")
   action.add_argument("--nocomments", dest="remove_comments", action="store_true",
                    help="Write fasta/fastq file with no comments in read name line")

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
      printout+="\nFactotum: Changing format: fasta from fastq";
      fq2fa=args.fq2fa
      preout+="fq2fa_"
     
   if args.seqname is not None:
       printout+="\nFactotum: Writing only seq named: " + args.seqname;
       seqname=args.seqname 
       preout+=seqname+"_"
       write=1
   if args.keepname is not None:
       printout+="\nFactotum: Writing only seqs matching " + args.keepname;
       keepname=args.keepname
       preout+=keepname+"_"
       write=1
   if args.rmname is not None:
       printout+="\nFactotum: Writing only seqs not matching: " + args.rmname;
       rmname=args.rmname
       preout+="removed_"+rmname+"_"
       write=1
 

   if args.min_length is not None:
       printout+="\nFactotum: Writing only seq longer than "+str(args.min_length);
       min_length=args.min_length
       preout+="min" + "_" + str(min_length) + "_"
       write=1
   if args.max_length is not None:
       printout+="\nFactotum: Writing only seq shorter than " + str(args.max_length);
       max_length=args.max_length
       preout+="max" +"_" + str(max_length) + "_"
       write=1

   if args.avoid is not None:
       printout+="\nFactotum: Writing only seqs not matching list  " + str(args.avoid);
       avoid=args.avoid
       preout+="avoid" +"_" + str(avoid) + "_"
       write=1

   if args.remove_comments:
       printout+="\nFactotum: Writing seqs with no comments" ;
       remove_comments=1
       preout+="nocomment" + "_"
       write=1

   if args.scaffbreak:
      printout+="\nFactotum: Breaking scaffolds @ 3Ns";
      preout+="ctgs"+"_"
      write=1
   if args.stats:
      printout+='\nFactotum: Calculating stats'
   if args.nlist:
      if args.nlist == 1: printout += '\n Factotum: Listing the longest seqs'
      else:  printout+='\nFactotum: Listing the ' + str(args.nlist) + ' longest seqs'
      printout+=' and the shortest one'
      nlist = args.nlist
      preout+="empty"+"_"
         
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
       parser.print_help()
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
      print(subprocess.check_output([exe, inputfile]).decode('utf-8').strip())
      
   if args.seqname or args.min_length or args.max_length or args.fq2fa or args.nlist or args.avoid or args.remove_comments or args.keepname or args.rmname:
      exe=factotum_bin+"/jolly"
      print("\n",subprocess.check_output([exe, inputfile, outfile, out_type, str(min_length), 
                                          str(max_length), str(nlist), seqname, keepname, rmname, avoid, str(remove_comments)]).decode('utf-8').strip())
      print(subprocess.check_output(['rm', '-f', 'empty_*']).decode('utf-8').strip())
   print(' ')



def filetype_(filename):
    global outfile
    global out_type

    typefasta=('fasta','fsa','fa','fna')
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




if __name__ == "__main__":
	main()
