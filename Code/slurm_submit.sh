#!/bin/bash
#
#SBATCH --partition=128x24
#SBATCH --job-name=eternal_sim_test
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=400
#SBATCH --time=00:30
#SBATCH --output=res.txt

srun python eternal_sim_leader.py --cores 4 --input_file infile_eternal.txt --output_file outfile.txt
