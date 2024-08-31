#!/bin/bash

#Joy Sung

echo -e "$(date)"

echo 'starting run'

COUNT_SOURCEDIR='./mapping/'
array=("$COUNT_SOURCEDIR"*ReadsPerGene.out.tab)

echo "Number of count files = ${#array[*]}"
#echo "Number of samples in target file ${#target_array[*]}"

#Read in names from first file (they all should be the same).
cat ${array[0]} | awk '{print $1}' > count_matrix.txt

endoffile='_nohmrRNA_noglobinReadsPerGene.out.tab'
SAMPLEARRAY=()

for item in ${array[@]}
do

    #Read in second column of each file and add it is a new column.
    echo $item
    file_name="${item}"
    sample_name=${file_name#$COUNT_SOURCEDIR}
    sample_name=${sample_name%$endoffile}
    #Fill array with sample names without path.
    SAMPLEARRAY+=($sample_name)
    eval "cat '$file_name' | awk '{print \$4}' | paste count_matrix.txt - > output.txt"

    mv output.txt count_matrix.txt

done

#Delete last 5 lines.
tail -n +5 count_matrix.txt > tmp.txt && mv tmp.txt count_matrix.txt

#Add header.
output=$(printf "\t%s" "${SAMPLEARRAY[@]}")
echo -e "Name${output}" | cat - count_matrix.txt > tmp.txt && mv tmp.txt count_matrix.txt

echo -e "$(date)"