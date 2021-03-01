#!/bin/bash
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=15
#SBATCH --time=24:00:00
#SBATCH --job-name variant_calling

module load NiaEnv/2019b
module load NiaEnv/2019b gnu-parallel
module load samtools/1.9 gatk/4.1.7
module load java
NCORES=4
parallel --joblog slurm-$SLURM_JOBID.log -j $NCORES --wd $PWD < commandlist.txt
