#!/bin/bash
#SBATCH --job-name=nps
#SBATCH --chdir=/scratch/lz62
#SBATCH --output=/scratch/lz62/uplands-ci/logs/FILENAME/nps-%A-%a.log
#SBATCH --time=08:00:00
#SBATCH --partition=all
#SBATCH --cpus-per-task=1
#SBATCH --mem=64000
#SBATCH --array=1-LINE_COUNT

module load geos
module load R/3.5.1-old
mkdir -p /scratch/lz62/uplands-ci/logs/FILENAME

name=$(sed -n "$SLURM_ARRAY_TASK_ID"p /scratch/lz62/uplands-ci/output/hpc/filelist-FILENAME)

srun Rscript /scratch/lz62/uplands-ci/dev/analysis.r $name
