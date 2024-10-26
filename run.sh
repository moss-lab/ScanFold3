#!/bin/bash -l
#SBATCH --partition=biocrunch         #Partition to submit to
#SBATCH --time=3-00:00:00             #Time limit for this job
#SBATCH --nodes=1                   #Nodes to be used for this job during runtime
#SBATCH --ntasks=25                 #Number of tasks to fire
#SBATCH --job-name=scanfold3_test   #Name of this job in work queue
#SBATCH --mail-user=  #Email to send notifications to
#SBATCH --mail-type=ALL             #Email notification type (BEGIN, END, FAIL, ALL)
#SBATCH --mem=25g

module load boost;
module load gcc;

g++ fold.cpp -o fold;
./fold ./test/coronaframeshift/fs.1.win_120.stp_1.tsv;
wait;

