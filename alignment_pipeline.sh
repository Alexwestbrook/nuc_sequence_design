#!/bin/bash

#### Download S.Cerevisiae repeats 601 ################################
# prefetch SRR13645559
# fasterq-dump SRR13645559 --outfile MNase_601/SRR13645559_167
# prefetch SRR13645549
# fasterq-dump SRR13645549 --outfile MNase_601/SRR13645549_197
# prefetch SRR13645545
# fasterq-dump SRR13645545 --outfile MNase_601/SRR13645545_237


#### Reference genome building ########################################
# cat genome/167_7_4kbrf.fa genome/YPH499_Stanford_2014_JRIO00000000.fsa > genome/YPH499_167.fa
# cat genome/197_7_4kbrf.fa genome/YPH499_Stanford_2014_JRIO00000000.fsa > genome/YPH499_197.fa
# cat genome/237_7_4kbrf.fa genome/YPH499_Stanford_2014_JRIO00000000.fsa > genome/YPH499_237.fa
# cat genome/167_601_7_4kbrf.fa genome/YPH499_Stanford_2014_JRIO00000000.fsa > genome/YPH499_167_601.fa
# cat genome/197_601_7_4kbrf.fa genome/YPH499_Stanford_2014_JRIO00000000.fsa > genome/YPH499_197_601.fa
# cat genome/237_601_7_4kbrf.fa genome/YPH499_Stanford_2014_JRIO00000000.fsa > genome/YPH499_237_601.fa
# for filename in genome/YPH499*.fa
# do
#     prefix=${filename%.*}
#     bowtie2-build $filename $prefix
# done


#### bowtie2 alignments ##############################################
# threads=16
# data_dir=MNases_06_2020
# ref='genome/YPH499_167'
# for prefix in 'NG-25315_167_1_lib406963_6879_2'
# do
#     cutadapt -j $threads -m 50 -O 1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o $data_dir/$prefix'_1_trimfiltered'.fastq.gz -p $data_dir/$prefix'_2_trimfiltered'.fastq.gz $data_dir/$prefix'_1'.fastq.gz $data_dir/$prefix'_2'.fastq.gz
#     out_prefix=$data_dir/$prefix'_trimfiltered_max250'
#     bowtie2 -p $threads -x $ref -1 $data_dir/$prefix'_1_trimfiltered.fastq.gz' -2 $data_dir/$prefix'_2_trimfiltered.fastq.gz' -X 250 -S $out_prefix.sam
#     bash alignment_postprocessing.sh -p $threads -f $out_prefix.sam
# done
# ref='genome/YPH499_197b'
# for prefix in 'NG-25315_1972_1_lib402286_6844_1'
# do
#     cutadapt -j $threads -m 50 -O 1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o $data_dir/$prefix'_1_trimfiltered'.fastq.gz -p $data_dir/$prefix'_2_trimfiltered'.fastq.gz $data_dir/$prefix'_1'.fastq.gz $data_dir/$prefix'_2'.fastq.gz
#     out_prefix=$data_dir/$prefix'_trimfiltered_max250'
#     bowtie2 -p $threads -x $ref -1 $data_dir/$prefix'_1_trimfiltered.fastq.gz' -2 $data_dir/$prefix'_2_trimfiltered.fastq.gz' -X 250 -S $out_prefix.sam
#     bash alignment_postprocessing.sh -p $threads -f $out_prefix.sam
# done
# ref='genome/YPH499_237'
# for prefix in 'NG-25315_237_1_lib406964_6879_2'
# do
#     cutadapt -j $threads -m 50 -O 1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o $data_dir/$prefix'_1_trimfiltered'.fastq.gz -p $data_dir/$prefix'_2_trimfiltered'.fastq.gz $data_dir/$prefix'_1'.fastq.gz $data_dir/$prefix'_2'.fastq.gz
#     out_prefix=$data_dir/$prefix'_trimfiltered_max250'
#     bowtie2 -p $threads -x $ref -1 $data_dir/$prefix'_1_trimfiltered.fastq.gz' -2 $data_dir/$prefix'_2_trimfiltered.fastq.gz' -X 250 -S $out_prefix.sam
#     bash alignment_postprocessing.sh -p $threads -f $out_prefix.sam
# done

# data_dir=MNase_08_2023
# ref='genome/YPH499_167'
# for prefix in 'NG-34039_167_4_1_lib713576_10294_1' 'NG-34039_167_4_2_lib713577_10294_1'
# do
#     cutadapt -j $threads -m 50 -O 1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o $data_dir/$prefix'_1_trimfiltered'.fastq.gz -p $data_dir/$prefix'_2_trimfiltered'.fastq.gz $data_dir/$prefix'_1'.fastq.gz $data_dir/$prefix'_2'.fastq.gz
#     out_prefix=$data_dir/$prefix'_trimfiltered_max250'
#     bowtie2 -p $threads -x $ref -1 $data_dir/$prefix'_1_trimfiltered.fastq.gz' -2 $data_dir/$prefix'_2_trimfiltered.fastq.gz' -X 250 -S $out_prefix.sam
#     bash alignment_postprocessing.sh -p $threads -f $out_prefix.sam
# done
# ref='genome/YPH499_197b'
# for prefix in 'NG-34039_197_2_27_1_lib713580_10294_1' 'NG-34039_197_2_27_2_lib713581_10294_1'
# do
#     cutadapt -j $threads -m 50 -O 1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o $data_dir/$prefix'_1_trimfiltered'.fastq.gz -p $data_dir/$prefix'_2_trimfiltered'.fastq.gz $data_dir/$prefix'_1'.fastq.gz $data_dir/$prefix'_2'.fastq.gz
#     out_prefix=$data_dir/$prefix'_trimfiltered_max250'
#     bowtie2 -p $threads -x $ref -1 $data_dir/$prefix'_1_trimfiltered.fastq.gz' -2 $data_dir/$prefix'_2_trimfiltered.fastq.gz' -X 250 -S $out_prefix.sam
#     bash alignment_postprocessing.sh -p $threads -f $out_prefix.sam
# done
# ref='genome/YPH499_237'
# for prefix in 'NG-34039_237_44_1_lib713582_10294_1' 'NG-34039_237_44_2_lib713583_10294_1'
# do
#     cutadapt -j $threads -m 50 -O 1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o $data_dir/$prefix'_1_trimfiltered'.fastq.gz -p $data_dir/$prefix'_2_trimfiltered'.fastq.gz $data_dir/$prefix'_1'.fastq.gz $data_dir/$prefix'_2'.fastq.gz
#     out_prefix=$data_dir/$prefix'_trimfiltered_max250'
#     bowtie2 -p $threads -x $ref -1 $data_dir/$prefix'_1_trimfiltered.fastq.gz' -2 $data_dir/$prefix'_2_trimfiltered.fastq.gz' -X 250 -S $out_prefix.sam
#     bash alignment_postprocessing.sh -p $threads -f $out_prefix.sam
# done

# data_dir=MNase_601
# ref='genome/YPH499_167_601'
# for prefix in 'SRR13645559_167'
# do
#     out_prefix=$data_dir/$prefix'_max250'
#     bowtie2 -p $threads -x $ref -1 $data_dir/$prefix'_1.fastq.gz' -2 $data_dir/$prefix'_2.fastq.gz' -X 250 -S $out_prefix.sam
#     bash alignment_postprocessing.sh -p $threads -f $out_prefix.sam
# done
# ref='genome/YPH499_197_601'
# for prefix in 'SRR13645549_197'
# do
#     out_prefix=$data_dir/$prefix'_max250'
#     bowtie2 -p $threads -x $ref -1 $data_dir/$prefix'_1.fastq.gz' -2 $data_dir/$prefix'_2.fastq.gz' -X 250 -S $out_prefix.sam
#     bash alignment_postprocessing.sh -p $threads -f $out_prefix.sam
# done
# ref='genome/YPH499_237_601'
# for prefix in 'SRR13645545_237'
# do
#     out_prefix=$data_dir/$prefix'_max250'
#     bowtie2 -p $threads -x $ref -1 $data_dir/$prefix'_1.fastq.gz' -2 $data_dir/$prefix'_2.fastq.gz' -X 250 -S $out_prefix.sam
#     bash alignment_postprocessing.sh -p $threads -f $out_prefix.sam
# done


#### Getting nucleosome densities with bamCoverage #################################
# threads=16
# for file in MNase_*/*_trimfiltered_max250.sorted.bam
# do
#     bamCoverage -p $threads -b $file -o "${file%%.*}"_140-170_CPM.bw --binSize 1 --minFragmentLength 140 --maxFragmentLength 170 --normalizeUsing CPM
# done

# data_dir=MNase_601
# ref='genome/YPH499_167_601'
# for prefix in 'SRR13645559_167'
# do
#     out_prefix=$data_dir/$prefix'_max250'
#     bamCoverage -p $threads -b $out_prefix.sorted.bam -o $out_prefix'_140-170_CPM'.bw --binSize 1 --minFragmentLength 140 --maxFragmentLength 170 --extendReads 151 --normalizeUsing CPM
# done
# ref='genome/YPH499_197_601'
# for prefix in 'SRR13645549_197'
# do
#     out_prefix=$data_dir/$prefix'_max250'
#     bamCoverage -p $threads -b $out_prefix.sorted.bam -o $out_prefix'_140-170_CPM'.bw --binSize 1 --minFragmentLength 140 --maxFragmentLength 170 --extendReads 161 --normalizeUsing CPM
# done
# ref='genome/YPH499_237_601'
# for prefix in 'SRR13645545_237'
# do
#     out_prefix=$data_dir/$prefix'_max250'
#     bamCoverage -p $threads -b $out_prefix.sorted.bam -o $out_prefix'_140-170_CPM'.bw --binSize 1 --minFragmentLength 140 --maxFragmentLength 170 --extendReads 153 --normalizeUsing CPM
# done