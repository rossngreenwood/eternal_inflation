#!/bin/bash
#
#SBATCH --partition=128x24
#SBATCH --job-name=etrnlinf
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=5000
#SBATCH --time=08:00:00
#SBATCH --output=../data/etrnlinf_%j.out
#SBATCH --mail-type=begin
#SBATCH --mail-type=fail
#SBATCH --mail-type=end
#SBATCH --mail-user=rossngreenwood@gmail.com

for test_id in "$@"
do
  infile="../data/input/infile_"
  infile+=$test_id
  infile+=".txt"
  
  outfile="outfile_"
  outfile+=$test_id
  outfile+=".txt"
  
  mkdir -p ../data/out_$test_id
  srun python eternal_sim_leader.py --cores 23 --input_file $infile --output_dir ../data/out_$test_id/ --output_file $outfile
done
