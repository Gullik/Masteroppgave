#!/bin/bash

for d in */; do
    # Will print */ if no directories are available
    echo $d
    cd $d
	for d in */; do
	    # Will print */ if no directories are available
#	    echo $d
	    cd $d
		for FILE in ./*; do
		  echo ${FILE}
		  pdfcrop "${FILE}" "${FILE}" 
		done
	    cd ..
	done
   cd ..
done

