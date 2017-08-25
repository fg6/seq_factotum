#include "readfastaq.h"
#include <locale>

int calc(void);

int main(int argc, char *argv[])
{ 

  if (argc < 3) {
   fprintf(stderr, "Usage: %s <reads.fq/fa> <out.file> [fastq/fasta] [min] [max] [list] [ctg] \n", 
	   argv[0]); //[ipos] [epos] [len_from_end]\n", argv[0]);
   return 1;
  }	
  if((fp = gzopen(argv[1],"r")) == NULL){ 
    printf("ERROR main:: missing input file  !! \n");
    return 1;
  }

  string seqfile = argv[1];
  string outname = argv[2];
  int writefile = 1;
  string ctg="";
  int minl=0;
  int maxl=0;
  int saveinfo=1; // save sequences in vectors
  int saveseq=0; // save sequences in vectors
  int readseq=1;
  int nlist=0;
  string otype="same";

  int ipos=0;
  int epos=0;
  int elen=0;  //len_from_end

  if(argc > 3) otype=argv[3];          // write fasta from fastq
  if(argc > 4) minl=atoi(argv[4]); 
  if(argc > 5) maxl=atoi(argv[5]);  // write ctg < maxl
  if(argc > 7) ctg=argv[7];          // write only ctg=ctg
  if(argc > 6) nlist=atoi(argv[6]);          


  // these not working yet in new version:
  //if(argc > 4) ipos=atoi(argv[4]);   // write only ctg from ipos base. if ipos==0 write from first base
  //if(argc > 5) epos=atoi(argv[5]);   // write only ctg up to epos base. If epos==0 write up to end of ctg
  //if(argc > 6) elen=atoi(argv[6]);   // write only ctg from ipos base. if ipos==0 write from first base
  std::stringstream ss;
  ss.imbue(std::locale(""));
  
  if(0)cout << argv[0] << " " << seqfile << " " << outname << " " 
	    << otype << " " << minl << " " << maxl << " " << ctg << " " << nlist << endl; 


  if (nlist > 0) writefile = 0;
  if(outname==seqfile){
    cout << " Error! trying to rewrite the same file! " <<endl;
    return 1;
  }


  int isfq=fasttype(argv[1]);
  if(otype=="fastq" && isfq) otype="same";
  int err=1;
  if(writefile)outfile.open(outname);
  if(!isfq){
    err=readfasta(argv[1],saveinfo,readseq,saveseq,minl,maxl,ctg,"same"); 
  }else{
    err=readfastq(argv[1],saveinfo,readseq,saveseq,minl,maxl,ctg,otype); 
  }
  if(writefile)outfile.close();

  
  int n, n50;
  long int max, bases, l50;
  float mean;
 
  if(!err){
    if(rlen.size()>0){  
      if(writefile) ss << "  Selected seq stats: "<< endl;  // for sel case
      std::tie(bases,n,mean, max, l50, n50)= calcstats(rlen);  
      ss << std::fixed << std::setprecision(0) <<  "  Bases= " << bases << " Seqs= "<< n << " Mean_length= " 
	<< mean << " Longest= " << max << " N50= "<< l50 << " N_n50= " << n50   //counting from 1
	<< std::endl;  
     
     if(nlist > 0){
	if (nlist > rlen.size()){
	  cout << " **** Warning: ***** \n There are only rlen.size() Seqs to list " << endl;
	  nlist=rlen.size();
	}

	cout << "\n Global seq stats: "<< endl;  // for list case
	sort(rlen.begin(),  rlen.end(), std::greater<int>());
	if (nlist==1)
	  ss << "\n The Longest Seq is:" << endl;
	else 
	  ss << "\n The " << nlist << " Longest Seqs are:" << endl;

	for (int nn=0; nn < nlist; nn++)
	  ss << " Chr " << rname[nn] << " lenght= " <<  std::fixed << rlen[nn] << " bp" << endl;

	int imin= *min_element(rlen.begin(),rlen.end());
	int ii = std::distance( rlen.begin(), std::find( rlen.begin(), rlen.end(), imin ) );
	ss << "\n The Shortest Seq is " << rname[ii] << " lenght= " << std::fixed << imin << " bp" <<endl;
      }
 
    }else{
      if(writefile) ss << "  Selected seq stats: "<< endl;  // for sel case
      if(writefile) ss << "  No seqs satisfied the requirements: no one selected " << endl;
    }
    if(excluded.size()>0){
      ss << "\n Excluded stats: "<< endl;
      std::tie(bases,n,mean, max, l50, n50)= calcstats(excluded);  
      ss << std::fixed << std::setprecision(0) <<  "  Bases= " << bases << " Seqs= "<< n << " Mean_length= " 
	<< mean << " Longest= " << max << " N50= "<< l50 << " N_n50= " << n50   //counting from 1
	<< std::endl;  
    }else{
      if(writefile) ss << "\n Excluded stats: "<< endl;
      if(writefile) ss << "  No seqs excluded " << endl; // only for selection cases
    }
    
    cout << ss.str() << endl;

  }else{
    cout << "  Error: Something went wrong while reading/writing the file.." << endl;
  }
  
  return 0;
}
