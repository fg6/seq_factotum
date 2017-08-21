#include <iomanip>  //setprecision
#include <algorithm>    // sort, reverse
#include <gzstream.h>
#include <vector>  //setprecision
#include <sstream>
#include <numeric> // accumulate
#include <tuple>


using std::cout;
using std::endl;
using std::vector;
using std::string;

static gzFile fp;  // input file
static std::ofstream outfile;  // output_file

static  vector<int> rlen;
static  vector<string> rseq;
static  vector<string> rqual;
static  vector<string> rname;
static  vector<string> rcomment;

static  vector<int> ctgrlen;


// ---------------------------------------- //
int fasttype(char* file)
// ---------------------------------------- //
{ 
  char fq[5]={"@"};
  char fa[5]={">"};
  string ttname;
  int isfq=0;
  igzstream infile(file);

  getline(infile,ttname);
  string ftype=ttname.substr(0,1);
  if(ftype==fa) isfq=0;
  else isfq=1;

  return(isfq);
}


// ---------------------------------------- //
int readfastq(char* file, int saveinfo=0, int readseq=0, int minlen=0, string selctg="")
// ---------------------------------------- //
{ 
  igzstream infile(file);
  char fq[5]={"@"};
  char fa[5]={">"};
  char plus[5]={"+"};
  int nseq=0;

 
  rlen.reserve(100000);
  rname.reserve(100000);
  if(readseq){
    rseq.reserve(100000);
    rqual.reserve(100000);
  }


  string read;
  string lname;
  string lcomment="";   
  string lseq="";
  int seqlen=0;
  int quallen=0;
  string lqual;
  int seqlines=0;


  int stop=1;
  while(stop){
    getline(infile,read);
    
    if(read.substr(0,1)==fq){  // name
      nseq++;

      if(nseq>1){ // previous
	if(seqlen>=minlen){

	  if(saveinfo){
	    rname.push_back(lname);
	    if(lcomment.size())rcomment.push_back(lcomment);
	    rlen.push_back(seqlen);
	    if(readseq)rseq.push_back(lseq);
	    if(readseq)rqual.push_back(lqual);
	  }


	  //cout << fa << lname ;
	  // if(lcomment.size()) cout << lcomment <<endl;
	  // else cout << endl;
	  
	  
	  if(quallen != seqlen)
	    cout << " ERROR! seq length different from quality lenght!! " << endl;
	}
      }
      
      size_t ns=0;
      size_t nt=0;
      ns=read.find(" ");
      nt=read.find("\t");

      if(ns!=std::string::npos) { 
	lname=read.substr(1,ns);
	lcomment=read.substr(ns+1,read.size());
      }else if(nt!=std::string::npos) {
	lname=read.substr(1,nt);
	lcomment=read.substr(nt+1,read.size());
      }else{
	lname=read.erase(0,1);
      }
      if(readseq){
	lseq="";
	lqual="";
      }
      seqlen=0;
      seqlines=0;
      quallen=0;
    }else if(read.substr(0,1)==plus){ // + and qual

      for(int ll=0; ll<seqlines; ll++){
	getline(infile,read);
	if(readseq)lqual.append(read);
	quallen+=read.size();
      }
    }else{ // sequence 
      if(readseq)lseq.append(read);
      seqlines++;
      seqlen+=read.size();
    }
 
    // EOF
    if(infile.eof()){ // previous
      if(seqlen>=minlen){

	if(saveinfo){
	  rname.push_back(lname);
	  if(lcomment.size())rcomment.push_back(lcomment);
	  rlen.push_back(seqlen);
	  if(readseq)rseq.push_back(lseq);
	  if(readseq)rqual.push_back(lqual);
	}
	if(quallen != seqlen)
	  cout << " ERROR! seq length different from quality lenght!! " << endl;
      }
      stop=0;
    }


  }//read loop
 
  return 0;
}


// ---------------------------------------- //
int readfasta(char* file, int saveinfo=0, int readseq=0, int minlen=0, string selctg="")
// ---------------------------------------- //
{ 
  igzstream infile(file);
  char fa[5]={">"};
  int nseq=0;

  rlen.reserve(100000);
  rname.reserve(100000);
  if(readseq)
    rseq.reserve(100000);


  string read;
  string lname;
  string lcomment="";   
  string lseq="";
  int seqlen=0;

  int stop=1;
  while(stop){
    getline(infile,read);
    
    if(read.substr(0,1)==fa){  // name
      nseq++;

      if(nseq>1){ // previous
	if(seqlen>=minlen){

	  if(saveinfo){
	    rname.push_back(lname);
	    rlen.push_back(seqlen);
	    if(readseq)rseq.push_back(lseq);
	    if(lcomment.size())rcomment.push_back(lcomment);
	  }


	}
      }
      
      size_t ns=0;
      size_t nt=0;
      ns=read.find(" ");
      nt=read.find("\t");

      if(ns!=std::string::npos) { 
	lname=read.substr(1,ns);
	lcomment=read.substr(ns+1,read.size());
      }else if(nt!=std::string::npos) {
	lname=read.substr(1,nt);
	lcomment=read.substr(nt+1,read.size());
      }else{
	lname=read.erase(0,1);
      }
      if(readseq)lseq="";
      seqlen=0;

    }else{ // sequence 
      if(readseq)lseq.append(read);
      seqlen+=read.size();
    }
 
    // EOF
    if(infile.eof()){ // previous
      if(seqlen>=minlen){


	if(saveinfo){
	  rname.push_back(lname);
	  rlen.push_back(seqlen);
	  if(readseq)rseq.push_back(lseq);
	  if(lcomment.size())rcomment.push_back(lcomment);
	}      
      }
      stop=0;
    }


  }//read loop
 
  return 0;
}





// numeric to string
template <class T>
inline std::string to_string (const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}
// string of numbers to int 
template <class T1>
inline int to_int (const T1& t)
{
  int ii;
  std::stringstream ss(t);
  ss >> ii;
  return ii;
}
// string of numbers to int 
template <class T2>
inline int to_float (const T2& t)
{
  float ii;
  std::stringstream ss(t);
  ss >> ii;
  float rr=ii*100.;
  return rr;
}


// ---------------------------------------- //
std::tuple<long int, int, float, long int, long int, int>  calcstats(vector<int> lengths)
// ---------------------------------------- //
{
  sort(lengths.begin(),  lengths.end(), std::greater<int>());

  int n=lengths.size();
  long int max=lengths[0];                 	
  long int  bases = accumulate(lengths.begin(), lengths.end(), 0.0);
  float mean = bases*1. / n;

  int n50=0;
  long int l50=0;
  int done=0;
  long int t50=0;
  int ii=0;
  while(done<1){
    t50+=lengths[ii];
    if(t50 > bases*0.5) 
      done=1;
    ii++;
   }

  n50=ii;
  l50=lengths[n50-1];  //counting from 0
  
  if(0)
    std::cout << std::fixed << std::setprecision(0) <<  "Bases= " << bases << " contigs= "<< n << " mean_length= " 
	<< mean << " longest= " << max << " N50= "<< l50 << " n= " << n50   //counting from 1
	<< std::endl;  

  return std::make_tuple(bases,n,mean, max, l50, n50);
}
