#!/bin/bash -l
#SBATCH --time=1-00:00:00             #Time limit for this job
#SBATCH --ntasks=16                #Number of tasks to fire
#SBATCH --job-name=scanfold3_demo   #Name of this job in work queue
#SBATCH --mem=10g


module load micromamba
micromamba activate /lustre/hdd/LAS/wmoss-lab/programs/envs/ScanFold3

python /lustre/hdd/LAS/wmoss-lab/scripts/ScanFold3/ScanFoldFold.py "HIV.fasta" --tsv "AF033819.3.win_120.stp_1.tsv";

