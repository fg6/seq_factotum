#!/bin/bash

myseq_handler=$(dirname $0)
if [[ $myseq_handler == "." ]]; then
    myseq_handler=`pwd`
fi
mysrcs=$myseq_handler/src
mybin=$myseq_handler/bin

# srcs to compile:
srcs=( n50 ctgs_from_scaff jolly )

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

errs=0
if [[ ! -f mylibs/gzstream/gzstream.o ]]; then
    echo cannot find gzstream: Error! 
    errs=$(($errs+1))
fi

cd $mysrcs
for code in "${srcs[@]}"; do 
    cd $mysrcs/$code
  
    if [[ ! -f $code ]] || [[ $code -ot $code.cpp ]] || [[ $code -ot $mysrcs/myincs/readfastaq.h ]]; then 
	make all 
	rm -f $mybin/$code
	cp $code $mybin/.
    fi
done


cd $mysrcs
echo; echo " All done."; echo " Checking installations:"

for exe in "${srcs[@]}"; do
    if [[ ! -f $mybin/$exe ]]  ||  [[ $mybin/$exe -ot $mysrcs/$exe/$exe ]]; then 

	if  [[ ! -f $mysrcs/$exe/$exe ]]; then 
            echo cannot find $mybin/$exe: Error! 
            errs=$(($errs+1))
	else 
	    cp $mysrcs/$exe/$exe $mybin/.
	fi
    fi
done

PATH=$mybin/:$PATH
if [  $errs -gt 0 ]; then echo " ****  Errors occurred! **** "; echo; exit; 
else echo " Congrats: installations successful!"; fi
