#include "readfastaq.h"

int calc(void);

int main(int argc, char *argv[])
{ 

  if (argc < 3) {
   fprintf(stderr, "Usage: %s <reads.fq/fa> <out.file> [fastq/fasta] [min] [max] [ctg]  \n", 
	   argv[0]); //[ipos] [epos] [len_from_end]\n", argv[0]);
   return 1;
  }	
  if((fp = gzopen(argv[1],"r")) == NULL){ 
    printf("ERROR main:: missing input file  !! \n");
    return 1;
  }

  string seqfile = argv[1];
  string outname = argv[2];
  string ctg="";
  int minl=0;
  int maxl=0;
  int saveinfo=0; // save sequences in vectors
  int readseq=1;
  string otype="same";

  int ipos=0;
  int epos=0;
  int elen=0;  //len_from_end

  if(argc > 3) otype=argv[3];          // write fasta from fastq
  if(argc > 4) minl=atoi(argv[4]); 
  if(argc > 5) maxl=atoi(argv[5]);  // write ctg < maxl
  if(argc > 6) ctg=argv[6];          // write only ctg=ctg

  // these not working yet in new version:
  //if(argc > 4) ipos=atoi(argv[4]);   // write only ctg from ipos base. if ipos==0 write from first base
  //if(argc > 5) epos=atoi(argv[5]);   // write only ctg up to epos base. If epos==0 write up to end of ctg
  //if(argc > 6) elen=atoi(argv[6]);   // write only ctg from ipos base. if ipos==0 write from first base

  if(outname==seqfile){
    cout << " Error! trying to rewrite the same file! " <<endl;
    return 1;
  }

  if(0)cout << argv[0] << " " << seqfile << " " << outname << " " 
	    << otype << " " << minl << " " << maxl << " " << ctg << endl; 


  int isfq=fasttype(argv[1]);
  if(otype=="fastq" && isfq) otype="same";
  
  int err=1;
  outfile.open(outname);
  if(!isfq){
    err=readfasta(argv[1],saveinfo,readseq,minl,maxl,ctg,"same"); 
  }else{
    err=readfastq(argv[1],saveinfo,readseq,minl,maxl,ctg,otype); 
  }
  outfile.close();

  
  if(!err){
    cout << "   Done."<< endl;
  }else{
    cout << "   Error: Something went wrong while reading/writing the file.." << endl;
  }
  
  return 0;
}
