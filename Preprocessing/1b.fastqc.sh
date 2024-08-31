#!/bin/bash

#Joy Sung

echo -e "$(date)"

RAW_SOURCEDIR='../Trimgalore_Results/'
Files=("$RAW_SOURCEDIR"*.fq.gz)
echo "Number samples = ${#Files[*]} trimmed files"

#Make a new directory.
mkdir Trimgalore_Fastqc_Results

#Run Fastqc.
counter=0
for item in ${Files[*]}
do

     echo "Counter variable $counter"

        /vol01/ngs_tools/FastQC/fastqc "$item" --noextract -t 16 -o Trimgalore_Fastqc_Results/ &

        counter=$((counter+1))
	wait
done

echo -e "$(date)"