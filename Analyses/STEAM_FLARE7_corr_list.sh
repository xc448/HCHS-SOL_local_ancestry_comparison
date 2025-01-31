#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=STEAM_corr_FLARE3  # Job name
#SBATCH --ntasks=1             # Run a single task
#SBATCH --time=96:00:00 
#SBATCH --cpus-per-task=8
#SBATCH --output="/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/GDS_corr_list/FLARE7_corr/%A_%a.log"
#SBATCH --error="/home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/GDS_corr_list/FLARE7_corr/%A_%a.out"
#SBATCH --array=[22]  
  
/usr/bin/Rscript /home/ec2-user/EBS4T/Projects/2023_local_ancestry_comparison/code/STEAM_FLARE7_corr_list.R $SLURM_ARRAY_TASK_ID 
