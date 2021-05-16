#include "readfastaq.h"
#include <locale>
#include <stdlib.h>     /* srand, rand */

#include <cstdlib>      // std::rand, std::srand
#include <random>

#include <vector>

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


  if (argc < 3) {
    cout << "Missing inputs!"<<endl;
    return 0;
  }

  int err=1;
  // Get readnames
  err=readfastaq(argv[1],1); // save info (contig names and length) in vectors

  // Get subsample of reads
  //srand (time(NULL));


  int initial_num_seq = rname.size();
  int subsample = atoi(argv[3]);

  cout <<endl<< "Select " << subsample << " Random Reads for files " <<argv[1] << " and "<< argv[2] <<endl<<endl;

  std::vector<string> filtered_read_names;
  filtered_read_names = rname;

  std::random_device rd;
  std::mt19937 g(rd());

  std::shuffle ( filtered_read_names.begin(), filtered_read_names.end(), g );
  filtered_read_names.resize(subsample);
 

  // Select reads R1 from readnames and write
  string outname = "subsample_1.fastq";
  outfile.open(outname);
  int sseqs=1;
  sseqs = readfastaq(argv[1], 1, 1, 1,  0,  0,  "",  "", 
                          0,  "",  "", filtered_read_names);
  outfile.close();


  // Select reads R2 from readnames and write
  outname = "subsample_2.fastq";
  outfile.open(outname);
  sseqs = readfastaq(argv[2], 1, 1, 1,  0,  0,  "",  "", 
                          0,  "",  "", filtered_read_names);
  outfile.close();
  

  return 0;
}

