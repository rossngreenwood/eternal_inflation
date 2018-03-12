#!/bin/bash
#
#SBATCH --partition=128x24
#SBATCH --job-name=eternal_sim_test
#SBATCH --output=res.txt
#SBATCH --ntasks=1
#SBATCH --time=04:30
#SBATCH --mem-per-cpu=400
#
#SBATCH --array=0-7

srun python eternal_sim_jobarray.py --task_id $SLURM_ARRAY_TASK_ID --cores 8 --input_file infile_eternal.txt --output_file outfile.txt

#python eternal_sim_cleanup.py --cores 8 --output_file outfile.txt
