#include "readfastaq.h"

int break_scaffolds(int write);
static  vector<int> ctglen;

int main(int argc, char *argv[])
{ 

  if (argc < 3) {
   fprintf(stderr, "Breaks scaffolds at site with more than 3 Ns\n", argv[0]);
   fprintf(stderr, "Usage: %s <reads.fa>  <command>\n", argv[0]);
   return 1;
  }	
  if((fp = gzopen(argv[1],"r")) == NULL){ 
    printf("ERROR main:: missing input file  !! \n");
    return 1;
  }
  
  string command = argv[2];   
  //if(argc > 3) ctg=argv[3];          
  //if(argc > 4) otype=argv[4];   

  int err=0;
  int saveinfo=1;
  int readseq=1;
  int write=0;
  string otype="same";

  if(command == "break")
    write=1;
    
 
  // Read file //
  int isfq=fasttype(argv[1]);
  if(otype=="fastq" && isfq) otype="same";
  if(!isfq){
    err=readfasta(argv[1],saveinfo,readseq); 
  }else{
    err=readfastq(argv[1],saveinfo,readseq); 
  }
  if(err) {
    cout << " Error while reading file "<< argv[1] << endl;
    return 1;
  }
  

  //if(command == "break"){  // break scaffolds //
  string myname = "ctgs.fasta"; 
  if(write)outfile.open(myname.c_str());
  err=break_scaffolds(write);
  if(write) outfile.close();
  
  if(err) {
    cout << " Error while breaking scaffolds " << endl;
    return 1;
  }else{

  }

  
  //}

  return 0;
}

// ---------------------------------------- //
int break_scaffolds(int write)
// ---------------------------------------- //
{
  char out[5]={">"};
  int scafcounts=0;
  int pri=0;
  vector<int> ctglen;

  for (unsigned i=0; i < rseq.size(); i++) { // loop over scaffolds
  
    string s = rseq[i];
    string name = rname[i];
    int len = rlen[i];
    string delim = "NNN"; // from Zemin suggestion
    
    auto start = 0U;
    auto end = s.find(delim);
    int l=0;
    if (end == std::string::npos){  // if no >= 3 Ns
      if(write) outfile <<  out[0] << name << " ctg_0_starting_at_pos_1" << endl;
      if(write) outfile << s << endl;
      ctglen.push_back(s.size());
    }else{  // if there are Ns:
      scafcounts++;
      while (end != std::string::npos)  // split scaffold at any group of Ns:
	{      
	  string contig= s.substr(start, end - start);
	  int moreNs=0;
	  if(!contig.empty()){
	    l++;
	    string thisname = name + "_ctg_" + to_string(l) + "_starting_at_pos_" + to_string(start+1);
	    if(start!=1){ // not for the first part of the contig
	      for (int mm=0; mm<contig.size(); mm++){
		if(contig[mm]=='N'){   
		  moreNs++;   // can at max be 2, otherwise is already taken out by the delim.length
		}else{
		  contig.erase(0, moreNs);
		  break;
		}
	      }
	    }
	    if(write) outfile <<  out[0] << thisname << endl;
	    if(write) outfile << contig << endl;
	    ctglen.push_back(contig.size());
	  }
	  start = end + delim.length();
	  end = s.find(delim, start);
	}
      if (end == std::string::npos){   // catch the last one, where there are no Ns anymore
	end=s.size();
	l++;
	string contig= s.substr(start, end - start);
	int moreNs=0;
	for (int mm=0; mm<contig.size(); mm++){
	  if(contig[mm]=='N'){
	    moreNs++;
	  }else{
	    contig.erase(0, moreNs);
	    start=start+moreNs;
	    break;
	  }
	}
	string thisname = name + "_ctg_" + to_string(l) + "_starting_at_pos_" + to_string(start+1);
	if(write) outfile <<  out[0] << thisname << endl;
	if(write) outfile << contig << endl;
	ctglen.push_back(contig.size());
      }
    }
  }



  if(scafcounts){
    cout << "\nI've found " << scafcounts << " breaking points of at least 3 Ns";
    if(write) cout << " (contigs written to ctgs.fasta) " << endl;
    else cout << endl;
  }else{
    cout << "\nFound no scaffolds to split!" << endl;
  }

  int n, n50;
  long int max, bases, l50;
  float mean;
  
  cout << "\nScaffold stats:" << endl;
  std::tie(bases,n,mean, max, l50, n50)= calcstats(rlen);  
  cout << std::fixed << std::setprecision(0) <<  " Bases= " << bases << " contigs= "<< n << " mean_length= " 
       << mean << " longest= " << max << " N50= "<< l50 << " n= " << n50   //counting from 1
       << std::endl;  
  
  if(scafcounts){
    cout << "\nContig stats:" << endl;
    std::tie(bases,n,mean, max, l50, n50)= calcstats(ctglen);  
    cout << std::fixed << std::setprecision(0) <<  " Bases= " << bases << " contigs= "<< n << " mean_length= " 
	 << mean << " longest= " << max << " N50= "<< l50 << " n= " << n50   //counting from 1
	 << std::endl;  
  }
  return 0;
}

