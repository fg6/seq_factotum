#!/bin/bash

myseq_handler=$(dirname $0)
if [[ $myseq_handler == "." ]]; then
    myseq_handler=`pwd`
fi
mysrcs=$myseq_handler/src

# srcs to compile and exes to check:
srcs=( n50 selctgs ctgs_from_scaff)
exes=( mylibs/gzstream/gzstream.o  n50/n50 selctgs/selctgs  ctgs_from_scaff/ctgs_from_scaff)

cd $mysrcs
mkdir -p mylibs

### Intalling gzstream (it needs zlib!)
if [[ ! -d  mylibs/gzstream ]]  || [[ ! -f mylibs/gzstream/gzstream.o ]]; then
    
    rm -rf mylibs
    mkdir mylibs
    cd mylibs
    
    wget https://www.cs.unc.edu/Research/compgeom/gzstream/gzstream.tgz 

    if [[ "$?" != 0 ]]; then
	echo "Error downloading gzstream, try again" 
	rm -rf gzstream gzstream.tgz 
	exit
    else
	tar -xvzf gzstream.tgz &> /dev/null
	if [[ "$?" != 0 ]]; then echo " Error during gzstream un-compressing. Exiting now"; exit; fi
	cd gzstream
	make &> /dev/null
	
	if [[ "$?" != 0 ]]; then echo " Error during gzstream compilation. Exiting now"; exit; fi
	test=`make test | grep "O.K" | wc -l`

	if [[ $test == 1 ]]; then echo " "1. gzstream installed; rm ../gzstream.tgz 
	else  echo  " Gzstream test failed. Exiting now"; exit; fi
    fi
fi

cd $mysrcs

for code in "${srcs[@]}"; do 
    echo $code 
    cd $mysrcs/$code
    make all
done


cd $mysrcs
echo; echo " Checking installations:"

errs=0
for exe in "${exes[@]}"; do
    if [[ ! -f $exe ]]; then 
        echo cannot find $exe: Error! 
        errs=$(($errs+1))
    fi
done
if [  $errs -gt 0 ]; then echo " ****  Errors occurred! **** "; echo; exit; 
else echo " Congrats: installation successful!"; fi




