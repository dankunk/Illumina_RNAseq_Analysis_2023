#!/usr/bin/env bash

#SBATCH --job-name=making_a_new_transcriptome
#SBATCH --nodes=1
#SBATCH --ntasks=24 # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --partition=amilan
#SBATCH --qos=normal     # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'
#SBATCH --time=7:30:00   # modify this to reflect how long to let the job go. This indicates 2.5 hours.
#SBATCH --output=2023_GFFread_making_novel_transcriptome-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dkunk1@colostate.edu

## This script will take a merged GTF from stringtie and a genome FASTA file and create a new transcriptome based on the novel GTF file made from stringtie. 
# This transcriptome can then be passed to salmon for quantifying gene expression. 

gffread -w /scratch/alpine/dkunk1@colostate.edu/RNAseq_ana_2023/genome/Rpadi_head_stringtie_merged.fa -g /scratch/alpine/dkunk1@colostate.edu/RNAseq_ana_2023/genome/R_padi_v2.fasta /scratch/alpine/dkunk1@colostate.edu/RNAseq_ana_2023/output_newTuxedo/04_string-tie/Rpadi_head_stringtie_merged.gtf