# nuc_sequence_design

## Description:

This is an update of the repository [nuc_sequence_design](https://github.com/etirouthier/nuc_sequence_design). The repository has been cleaned of unnecessary scripts to be more understandable. It is meant as source code for the article *In silico design of DNA sequences for in vivo nucleosome positioning*. You may find notebooks to reproduce the figures, and the code used to generate the sequences.

The repository is essentially aimed at designing a DNA sequence thats attract nucleosome in yeast using the kinetic Monte-Carlo methodology (KMC) ! 

The initial sequence will evolve step by step to a nucleosome positioning sequence (i.e a sequence where a nucleosome is set on the first 147 bp). At every step the algorithm will choose a mutation that tends to favor the positioning power of the sequence. To do so it uses internally a deep learning model able to predict the nucleosome density associated with any sequence in yeast.

------------------


## Getting started:

To design a sequence use the script `sequence_design_all_mut_w_rev_and_weights_w_repeat.py` in `Programme`.

The command are as follows :
  - -l or --length : length of the sequence (default is 167)
  - -r or --repeat : number of repeats to predict on (default is 1)
  - -s or --steps : number of KMC iterations (default is 100)
  - -m or --model : a hdf5 file containing the model to predict with (default is weights_with_rev_compl_rep2.hdf5)
  - -t or --temperature : the "temperature" to use in the model (default is 0.1)
  - -d or --directory : the output directory (default is test)
  - -w or --weights : the weights used in the energy function [a, b, c, d] with a for the CG content energy, b for the positioning energy in the direct strand, c for the energy in the reverse strand and d for the mutation energy (default [1, 1, 1, 1])
  - -g or --gauss : amplitude and background of the gaussian target (default [0.4, 0.2])

The model file must be stored in `Results_nucleosome`, where the model used in this study is provided.

  ------------------


## Results:

Two directory are created with the specified name in `Results_nucleosome`.

The first one contains a config.txt file to store the chosen configuration, a energy.txt file with the energy at every steps of the KMC optimisation and the sequences of each step of the optimization. The second directory stores plots of the predicted nucleosome density at each step.


## Other Resources

In `genome`, you may find fasta files of the designed sequences that where kept for further investigation in the study. They are made of 7 repetitions of the designed sequence, flanked on each side by the 4kb region where they were inserted in the genome of *S. Cerevisiae*.

The bash scripts `alignment_pipeline.sh` and `alignment_postprocessing.sh` contain the command lines used to derive the nucleosome densities from the raw MNase fastq files.

Finally, you may also find the script `predict_on_repeats.py` which extracts the designed sequences from the fasta files to make predictions with the deep learning model.

## Environment
All python scripts have been run with python3.8 with the following packages:
- numpy 1.22.3
- pandas 1.3.5
- tensorflow 2.5.0
- matplotlib 3.5.0
- seaborn 0.11.2
- pysam 0.16.0.1
- pyBigWig 0.3.18

For using tensorflow with GPU support, we used cudatoolkit 11.3.1 for a GPU NVIDIA GeForce RTX 2080 Ti

For the bash scripts, the following packages have been used:
- sratoolkit 3.0.0
- cutadapt 3.4
- bowtie 2.2.5
- samtools 1.10
- bamCoverage 3.5.1