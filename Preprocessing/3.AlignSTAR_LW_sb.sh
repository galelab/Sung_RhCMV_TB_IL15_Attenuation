#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 39
#SBATCH -t 5800-12
#SBATCH --partition aRNA-Seq
#SBATCH -w roundworm
#AUTHOR: Leanne Whitmore (running for Joy Sung)

##This loads STAR from lwhitmo home directory (Currently it looks like most uptodate star was not loaded)
module load STAR/2.7.10b

####BATCH --partition HoldingPen

STAR --version

mkdir "mapping"
mkdir "mapping/logs"

###Generates variables for paths to raw seq data and reference genome (make sure to have / at end of paths)

RAW_SOURCEDIR='./nohmrRNA_noglobin/'

#Dan Newhouse generate indexes for this genome with option sjdboverhang 99 specified
GENOME_SOURCEDIR='/vol01/genome/Macaca_mulatta/10.109/STAR_2.7.5a_indicies'
GTF_FILE='/vol01/genome/Macaca_mulatta/10.109/Macaca_mulatta.Mmul_10.109.gtf'

##Pulls in all sequencing data from raw source directory (NOTE: *D* ensures no control files are pulled out)
R1_files=("$RAW_SOURCEDIR"*fastq.1.gz)
R2_files=("$RAW_SOURCEDIR"*fastq.2.gz)

#Prints number of samples for reads 1 and reads 2 should be equal
echo "Number samples = ${#R1_files[*]} R1 files"
echo "Number samples = ${#R2_files[*]} R2 files"

echo "STATUS: Aligning TB samples ..."
counter=0

while [ "$counter" -lt  ${#R1_files[*]} ]
do
        echo "Counter variable $counter"
	##--genomeDir - directory to star indexes for reference genome
        ##--clip5pNbases - number of bases to clip off of 5 prime end of reads (both reads 1, and 2): note default is 0
        ##--clip3pNbases - number of bases to clip off of 3 prime end of reads (both reads 1, and 2): note default is 0
        ##--readFilesCommand - Tells STAR read files are compressed zcat is to be specified if files gz and gunzip -c if files are bzip2 files
        ##--readFilesIn - specifies read files (need 2 for paired end)
        ##--outSAMtype - type of alignment file to outline
        ##--outFileNamePrefix - location and name of outputfile
        ##--runThreadN - number of threads/processors for STAR to use in alignment
	##--quantMode=GeneCounts - counts reads per gene and outputs read counts to file ReadsPerGene.out.tab
	##--sjdbGTFfile - specifies path to GTF file
 
        if [ "$counter" -lt ${#R1_files[*]} ]
	then
		##1.Removes read and file information from file name (i.e will remove .fastq.1.gz)
        	sample_name=${R1_files[$counter]}
		sample_name=${sample_name%.fastq.1.gz}
        	##2.Removes /path/to/file/ (these two steps are done so we have an outputfile name for mapping results)
        	sample_name=${sample_name#$RAW_SOURCEDIR}
		echo "Sample being processed $sample_name"
        	echo "Read 1 file ${R1_files[$counter]}"
        	echo "Read 2 file ${R2_files[$counter]}"
        	srun -c 13 STAR  --genomeDir $GENOME_SOURCEDIR --sjdbGTFfile $GTF_FILE --clip5pNbases 1 1 --readFilesCommand zcat --outSAMtype BAM Unsorted --quantMode=GeneCounts --readFilesIn  ${R1_files[$counter]} ${R2_files[$counter]} --outFileNamePrefix ./mapping/"$sample_name" --runThreadN 11 1>./mapping/logs/"$sample_name"_mapping.log 2>&1 &
		counter=$((counter+1))
        fi

        # if [ "$counter" -lt ${#R1_files[*]} ]
        # then
        #         ##1.Removes read and file information from file name (i.e will remove .fastq.1.gz)
        #         sample_name=${R1_files[$counter]}
        #         sample_name=${sample_name%.fastq.1.gz}
        #         ##2.Removes /path/to/file/ (these two steps are done so we have an outputfile name for mapping results)
        #         sample_name=${sample_name#$RAW_SOURCEDIR}
        #         echo "Sample being processed $sample_name"
        #         echo "Read 1 file ${R1_files[$counter]}"
        #         echo "Read 2 file ${R2_files[$counter]}"
        #         srun -c 13 STAR --genomeDir $GENOME_SOURCEDIR --sjdbGTFfile $GTF_FILE --clip5pNbases 1 1  --readFilesCommand zcat --outSAMtype BAM Unsorted --quantMode=GeneCounts --readFilesIn  ${R1_files[$counter]} ${R2_files[$counter]} --outFileNamePrefix ./mapping/"$sample_name" --runThreadN 11 1>./mapping/logs/"$sample_name"_mapping.log 2>&1 &
        #         counter=$((counter+1))
        # fi

        if [ "$counter" -lt ${#R1_files[*]} ]
        then
		##1.Removes read and file information from file name (i.e will remove .fastq.1.gz)
                sample_name=${R1_files[$counter]}
                sample_name=${sample_name%.fastq.1.gz}
        	##2.Removes /path/to/file/ (these two steps are done so we have an outputfile name for mapping results)
        	sample_name=${sample_name#$RAW_SOURCEDIR}
        	echo "Sample being processed $sample_name"
        	echo "Read 1 file ${R1_files[$counter]}"
        	echo "Read 2 file ${R2_files[$counter]}"
        	srun -c 13 STAR --genomeDir $GENOME_SOURCEDIR --sjdbGTFfile $GTF_FILE --clip5pNbases 1 1  --readFilesCommand zcat --outSAMtype BAM Unsorted --quantMode=GeneCounts --readFilesIn  ${R1_files[$counter]} ${R2_files[$counter]} --outFileNamePrefix ./mapping/"$sample_name" --runThreadN 11 1>./mapping/logs/"$sample_name"_mapping.log 2>&1 &
		counter=$((counter+1))
        fi

        if [ "$counter" -lt ${#R1_files[*]} ]
        then

		##1.Removes read and file information from file name (i.e will remove .fastq.1.gz)
                sample_name=${R1_files[$counter]}
                sample_name=${sample_name%.fastq.1.gz}
        	##2.Removes /path/to/file/ (these two steps are done so we have an outputfile name for mapping results)
        	sample_name=${sample_name#$RAW_SOURCEDIR}
        	echo "Sample being processed $sample_name"
		echo "Read 1 file ${R1_files[$counter]}"
        	echo "Read 2 file ${R2_files[$counter]}"
        	srun -c 13 STAR --genomeDir $GENOME_SOURCEDIR --sjdbGTFfile $GTF_FILE --clip5pNbases 1 1 --readFilesCommand zcat --outSAMtype BAM Unsorted --quantMode=GeneCounts --readFilesIn  ${R1_files[$counter]} ${R2_files[$counter]} --outFileNamePrefix ./mapping/"$sample_name" --runThreadN 11 1>./mapping/logs/"$sample_name"_mapping.log 2>&1 &
	fi
	wait
        counter=$((counter+1))

done
# STAR --genomeLoad Remove --genomeDir $GENOME_SOURCEDIR

python STARresultssummary.py --folder mapping/
mutt -s "STAR genome alignment done" lwhitmo@uw.edu -c lwhitmo@uw.edu