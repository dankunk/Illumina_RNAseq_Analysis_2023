#!/usr/bin/env bash

#SBATCH --job-name=execute_hisat2-build
#SBATCH --nodes=1
#SBATCH --ntasks=15 # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --partition=amilan
#SBATCH --qos=normal     # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'
#SBATCH --time=2:30:00   # modify this to reflect how long to let the job go. This indicates 2.5 hours.
#SBATCH --output=log_hisat2-build_Rpadi-%j.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dkunk1@colostate.edu


# Test that you have hisat2 working
#hisat2

# Extract splice sites from intron containing transcript structures
#hisat2_extract_splice_sites.py Rpadi_v2.gff > Rpadi_v2.gff.ss

# Extract the exon records
#hisat2_extract_exons.py  Rpadi_v2.gff > Rpadi_v2.gff.exons

# Build hisat2 indexes for R padi
#hisat2-build -p ${SLURM_NTASKS} --exon Rpadi_v2.gff.exons --ss Rpadi_v2.gff.ss R_padi_v2.fasta rpv2

#Build hisat2 indexes for R padi (HFM index)
hisat2-build -p ${SLURM_NTASKS} R_padi_v2.fasta rpv2

# Check the build
echo -e "\n\nINDEX-BUILD: inspecting indexes:"
hisat2-inspect -s rpv2 

# Capture version number
echo -e "\n\nINDEX-BUILD: version:"
hisat2-build --version
