#!/bin/bash

#Joy Sung

echo -e "$(date)"

#Check versions.
/share/lwhitmo/seqsoftware/bin/trim_galore --version
/share/Joy/py3venv_extra/bin/cutadapt --version

#Environment.
source /share/Joy/py3venv_extra/bin/activate

RAW_SOURCEDIR='../final-fastq-set/'
R1_files=("$RAW_SOURCEDIR"*R1*.fastq.gz)
R2_files=("$RAW_SOURCEDIR"*R2*.fastq.gz)
echo "Number samples = ${#R1_files[*]} R1 files"
echo "Number samples = ${#R2_files[*]} R2 files"

mkdir Trimgalore_Results

counter=0
for R1file in ${R1_files[*]}
do

     echo "Counter variable $counter"

     #Trim Galore Parameters#
     #-q: Trims low quality ends from reads (defualt phred score is 20)
     #--phred33: instructs cutadapt to use ASCII+33 quality scores as Phred scores (this is the defualt)
     #--fastqc: run fastqc
     #--output_dir: directory to put trimmed fastq files
     #no adapter sequence or option specified which means it runs autodetect (looks for Illumina universal, Nextera transposase or Illumina small RNA adapter sequences)

        echo "Trim Galore running..."
        echo "${R1_files[$counter]}"
        echo "${R2_files[$counter]}"
        /share/lwhitmo/seqsoftware/bin/trim_galore -q 20 --gzip --cores 8 --path_to_cutadapt /share/Joy/py3venv_extra/bin/cutadapt --paired ${R1_files[$counter]} ${R2_files[$counter]} --output_dir Trimgalore_Results &

        counter=$((counter+1))
	wait
done

echo -e "$(date)"