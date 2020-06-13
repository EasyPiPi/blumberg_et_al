#!/bin/bash


##USAGE - sh runRNAanalysis.sh inputdirectory outputdirectory logdirectory
scriptdir=$1
inputdir=$2
out=$3
log=$4

for inputfile in `ls $inputdir`
do
echo "$inputfile start!"
file_name=`echo $inputfile | cut -d '.' -f1`
Rscript $scriptdir/RNA_SEM.R $inputdir/$inputfile $out/${file_name}_out.csv > $log/${file_name}_log.txt
echo "$inputfile done!"
done
