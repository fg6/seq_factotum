#include "/nfs/users/nfs_f/fg6/ana/cpp/myinclude/macro.h"
#include "/nfs/users/nfs_f/fg6/ana/cpp/myinclude/readnwritefaq.h"

int break_scaffolds(void);

int main(int argc, char *argv[])
{ 

  if (argc == 1) {
   fprintf(stderr, "Breaks scaffolds at site with more than 3 Ns\n", argv[0]);
   fprintf(stderr, "Usage: %s <reads.fa>\n", argv[0]);
   return 1;
  }	
  if((fp = gzopen(argv[1],"r")) == NULL){ 
    printf("ERROR main:: missing input file  !! \n");
    return 1;
  }
  
  int err=1;
  int saveinfo=1;
  int readseq=1;
  string otype="same";
  int isfq=fasttype(argv[1]);
  if(otype=="fastq" && isfq) otype="same";
  if(!isfq){
    err=readfasta(argv[1],saveinfo,"",readseq); 
  }else{
    err=readfastq(argv[1],saveinfo,"",readseq); 
  }
  if(err) {
    cout << " Error reading file "<< argv[1] << endl;
    return 1;
  }


  string myname = "ctgs.fasta"; 
  int write=1;
  if(write)myfile.open(myname.c_str());
  err=break_scaffolds();
  if(write) myfile.close();

  if(err) {
    cout << " Error while breaking scaffolds " << endl;
    return 1;
  }

  return 0;
}

// ---------------------------------------- //
int break_scaffolds()
// ---------------------------------------- //
{
  char out[5]={">"};
  int scafcounts=0;
  int pri=0;

  for (unsigned i=0; i < rseq.size(); i++) { // loop over scaffold
  
    string s = rseq[i];
    string name = rname[i];
    int len = rlen[i];
    string delim = "NNN"; // from Zemin suggestion
    
    auto start = 0U;
    auto end = s.find(delim);
    int l=0;
    if (end == std::string::npos){  // if no >= 10 Ns
      myfile <<  out[0] << name << " ctg_0_starting_at_pos_1" << endl;
      myfile << s << endl;
    }else{  // if there are Ns:
      scafcounts++;
      while (end != std::string::npos)  // split scaffold at any group of Ns:
	{      
	  string contig= s.substr(start, end - start);
	  int moreNs=0;
	  if(!contig.empty()){
	    l++;
	    string thisname = name + "_ctg_" + to_string(l) + "_starting_at_pos_" + to_string(start+1);
	    //cout << contig << endl;
	    if(start!=1){ // not for the first part of the contig
	      for (int mm=0; mm<contig.size(); mm++){
		if(contig[mm]=='N'){   
		  moreNs++;   // can at max be 2, otherwise is already taken out by the delim.length
		  //cout << contig[mm] << endl;
		}else{
		  contig.erase(0, moreNs);
		  //cout << " found " << moreNs << " more Ns" << endl;
		  break;
		}
	      }
	    }
	    myfile <<  out[0] << thisname << endl;
	    myfile << contig << endl;
	  }
	  start = end + delim.length();
	  end = s.find(delim, start);
	}
      if (end == std::string::npos){   // catch the last one, where there are no Ns anymore
	end=s.size();
	l++;
	string contig= s.substr(start, end - start);
	int moreNs=0;
	//cout << contig << endl;
	for (int mm=0; mm<contig.size(); mm++){
	  if(contig[mm]=='N'){
	    moreNs++;
	    //cout << contig[mm] << endl;
	  }else{
	    contig.erase(0, moreNs);
	    start=start+moreNs;
	    //cout << " found " << moreNs << " more Ns" << endl;
	    break;
	  }
	}
	string thisname = name + "_ctg_" + to_string(l) + "_starting_at_pos_" + to_string(start+1);
	myfile <<  out[0] << thisname << endl;
	myfile << contig << endl;
      }
    }
  }
  

  cout << " I've found " << scafcounts << " scaffold(s) to split" << endl;
  return 0;
}
