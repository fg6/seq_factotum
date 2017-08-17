
//#include "/nfs/users/nfs_f/fg6/ana/cpp/myinclude/readfaq.h"
#include "/nfs/users/nfs_f/fg6/ana/cpp/myinclude/macro.h"
#include "/nfs/users/nfs_f/fg6/ana/cpp/myinclude/readnwritefaq.h"

int calc(void);

int main(int argc, char *argv[])
{ 

  if (argc == 1) {
   fprintf(stderr, "Usage: %s <reads.fq/fa>  [minL] [ctg] [fastq/fasta]\n", argv[0]);
 //[ipos] [epos] [len_from_end]\n", argv[0]);
   
   return 1;
  }	
  if((fp = gzopen(argv[1],"r")) == NULL){ 
    printf("ERROR main:: missing input file  !! \n");
    return 1;
  }

  string seqfile = argv[1];
  string ctg="";
  int minl=0;
  int saveinfo=0;
  int readseq=0;
  string otype="same";

  int ipos=0;
  int epos=0;
  int elen=0;  //len_from_end

  if(argc > 2) minl=atoi(argv[2]);   // write ctg > minl
  if(argc > 3) ctg=argv[3];          // write only ctg=ctg
  if(argc > 4) otype=argv[4];          // write fasta from fastq

  // these not working yet in new version:
  //if(argc > 4) ipos=atoi(argv[4]);   // write only ctg from ipos base. if ipos==0 write from first base
  //if(argc > 5) epos=atoi(argv[5]);   // write only ctg up to epos base. If epos==0 write up to end of ctg
  //if(argc > 6) elen=atoi(argv[6]);   // write only ctg from ipos base. if ipos==0 write from first base

  
  //string mm="min" + to_string(minl) + "_";
  string mm="max" + to_string(minl) + "_";
  string myname=seqfile;
  string cc=ctg+"_";
  if(minl) myname=myrename(myname,mm);
  if(ctg.size()) myname=myrename(myname,cc);
  if(otype=="fasta")myname=myrename(myname,"",".fasta");
  else if(otype=="fastq")myname=myrename(myname,"",".fastq");

  int err=1;
  if(myname==seqfile){
    cout << " Error! trying to rewrite the same file! " <<endl;
    return 1;
  }


  int isfq=fasttype(argv[1]);
  if(otype=="fastq" && isfq) otype="same";


  myfile.open(myname);
  if(!isfq){
    err=readfasta(argv[1],saveinfo,"same",readseq,minl,ctg); //"same" for writing
  }else{
    err=readfastq(argv[1],saveinfo,otype,readseq,minl,ctg); //"same" for fastq, "fasta" for fasta
  }
  myfile.close();


  if(!err){
    if (minl) cout <<  "  Selected only contigs with length>=" << minl << endl;
    if(ctg.size())
      cout <<  "  Selected only contig named " << ctg << endl;
    //if(ipos || epos)
    // cout << "  Selected only range of bases from position " << ipos
    //	   << " to position " << epos << endl;
    cout << " Output file: " << myname << endl;
  }else{
    cout << " Sorry, something went wrong..." << endl;
  }
  
  return 0;
}
