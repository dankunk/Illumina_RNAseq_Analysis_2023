#!/usr/bin/env bash

#SBATCH --job-name=uzip_samples_Rpadi
#SBATCH --nodes=1
#SBATCH --ntasks=20 # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --partition=amilan
#SBATCH --qos=normal     # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'
#SBATCH --time=5:30:00   # modify this to reflect how long to let the job go. This indicates 2.5 hours.
#SBATCH --output=uzip_samples_Rpadi-%j.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dkunk1@colostate.edu


bash uzip_script.sh $SLURM_NTASKS