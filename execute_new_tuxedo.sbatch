#!/usr/bin/env bash

#SBATCH --job-name=2023_new-Tuxedo_rpadi
#SBATCH --nodes=1
#SBATCH --ntasks=24 # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --partition=amilan
#SBATCH --qos=normal     # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'
#SBATCH --time=7:30:00   # modify this to reflect how long to let the job go. This indicates 2.5 hours.
#SBATCH --output=2023_new-Tuxedo_rpadi-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dkunk1@colostate.edu

######### Instructions ###########

# Modify your SLURM entries above to fit your choices
# Below, modify the SECOND argument to point to YOUR metadata.file
# Below, you don't need to change $SLURM_NTASKS. It will automatically populate whatever you put in --ntasks=# above.
# Execute this script using $ sbatch execute_new_tuxedo.sbatch



## Execute the RNA-seq_pipeline to run the pipeline
bash dk_RNAseq_newtuxedo_pipeline.sh /scratch/alpine/dkunk1@colostate.edu/RNAseq_ana_2023/scripts/rpadi_metadata.txt $SLURM_NTASKS


