#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name=sample_response_drug_response_CTRPv2_pilot
#SBATCH --time=0-03:00:00
#SBATCH --partition=exacloud
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=./output_reports/slurm.%N.%j.out
#SBATCH --error=./error_reports/slurm.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=karlberb@ohsu.edu
#SBATCH -a 56-100:1

python ctrpv2_20_percent.py $SLURM_ARRAY_TASK_ID
