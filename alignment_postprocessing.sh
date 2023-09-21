#!/bin/bash

# usage:
#   alignment_postprocessing.sh -f <sam_filename>

threads=1
while getopts "f:p:" option; do
    case $option in
        f) # sam file
            filename=$OPTARG;;
        p) # number of threads to use
            threads=$OPTARG;;
        \?) # Invalid option
            echo "Error: Invalid option"
            exit;;
    esac
done

prefix=${filename%.*}
if [ $prefix.sam != $filename ]
then
    echo 'not a sam file!'
fi
samtools view -bS $filename | samtools sort -o $prefix.sorted.bam
rm $filename
samtools index $prefix.sorted.bam

samtools view -f66 $prefix.sorted.bam | cut -f9 | awk '{print sqrt($0^2)}' > $prefix.sorted.bam.insert_sizes.txt
samtools view -Mh -L repeats.bed $prefix.sorted.bam | samtools view -f66 | cut -f9 | awk '{print sqrt($0^2)}' > $prefix.sorted.bam.repeat_insert_sizes.txt