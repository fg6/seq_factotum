#include "readfastaq.h"


int main(int argc, char *argv[])
{ 

  if (argc == 1) {
   fprintf(stderr, "Calculate fasta/fastq stats\n", argv[0]);
   fprintf(stderr, "Usage: %s <file.fq/fa>\n", argv[0]);
   return 1;
  }	
  if((fp = gzopen(argv[1],"r")) == NULL){ 
    printf("ERROR main:: missing input file  !! \n");
    return 1;
  }

  int err=1;

  // File type	
  int isfq=fasttype(argv[1]);

  if(!isfq)
     err=readfasta(argv[1],1); // save info (contig names and length) in vectors
  else
     err=readfastq(argv[1],1);

  if(!err){
    int n, n50;
    long int max, bases, l50;
    float mean;

    std::tie(bases,n,mean, max, l50, n50)= calcstats(rlen);  
    std::cout << std::fixed << std::setprecision(0) <<  "Bases= " << bases << " contigs= "<< n << " mean_length= " 
	<< mean << " longest= " << max << " N50= "<< l50 << " n= " << n50   //counting from 1
	<< std::endl;  
  }
  return 0;
}

