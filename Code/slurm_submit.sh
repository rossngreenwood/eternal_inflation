#!/bin/batch
#
#SBATCH -p Instruction   # Partition name
#SBATCH --job_name=test
#SBATCH --output=res.txt
#SBATCH --ntasks=2
#SBATCH --time=00:30
#SBATCH --mem-per-cpu=100

srun run_mat_shell
