#!/usr/bin/env bash

#SBATCH --job-name=2023_GFF_compare_merged_vs_reference
#SBATCH --nodes=1
#SBATCH --ntasks=24 # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --partition=amilan
#SBATCH --qos=normal     # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'
#SBATCH --time=7:30:00   # modify this to reflect how long to let the job go. This indicates 2.5 hours.
#SBATCH --output=2023_GFFcompare_merged_vs_ref-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dkunk1@colostate.edu


## A script to compare a reference GTF file to a merged GTF file produced by stringtie. The merged file has isoforms that are novel and not found in the original reference. 

gffcompare -V -R -r /scratch/alpine/dkunk1@colostate.edu/RNAseq_ana_2023/genome/Rpadi_v2.gff -o strtie_merged_comparison /scratch/alpine/dkunk1@colostate.edu/RNAseq_ana_2023/output_newTuxedo/04_string-tie/Rpadi_head_stringtie_merged.gtf