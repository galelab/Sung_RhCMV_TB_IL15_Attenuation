#!/bin/bash

#Joy Sung

echo -e "$(date)"

bowtie2 --version

BOWTIE2_LIBRARIES="/vol01/genome/rRNA/bowtie2_index_Rhesus_05262021/human_mouse_rhesus_05262021_rRNA_v2"

RAW_SOURCEDIR='./Trimgalore_Results/'
R1_files=("$RAW_SOURCEDIR"*R1*.fq.gz)
R2_files=("$RAW_SOURCEDIR"*R2*.fq.gz)
echo "Number samples = ${#R1_files[*]} R1 files"
echo "Number samples = ${#R2_files[*]} R2 files"

mkdir nohmrRNA_noglobin
mkdir logs
mkdir logs/rRNA_globinfilter

echo -e "$(date)"

#Run bowtie2.
counter=0
for R1file in ${R1_files[*]}
do

     echo "Counter variable $counter"

                sample_name=${R1_files[$counter]}
                sample_name=${sample_name%_R*}
                sample_name=${sample_name#$RAW_SOURCEDIR}

                echo "Sample being processed $sample_name"
                echo "Read 1 file ${R1_files[$counter]}"
                echo "Read 2 file ${R2_files[$counter]}"

                cmd='bowtie2 -p 16 -5 1 -x $BOWTIE2_LIBRARIES --un-conc-gz ./nohmrRNA_noglobin/"$sample_name"_nohmrRNA_noglobin.fastq.gz -1 ${R1_files[$counter]} -2 ${R2_files[$counter]} -S ./nohmrRNA_noglobin/"$sample_name".sam 1>>./logs/rRNA_globinfilter/"$sample_name"_rRNA_globinfilter.log 2>&1 &'

                echo $cmd
                echo $cmd > ./logs/rRNA_globinfilter/"$sample_name"_rRNA_globinfilter.log
                eval $cmd

                counter=$((counter+1))
	wait
done

echo -e "$(date)"