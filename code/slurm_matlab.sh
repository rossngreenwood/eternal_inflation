#!/bin/bash
#
#SBATCH --partition=128x24
#SBATCH --job-name=etrnlinf
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=1000
#SBATCH --time=02:00:00
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
  
  cmd="eternal_sim_leader("
  cmd+="24,'"
  cmd+=$infile
  cmd+="','"
  cmd+=$outfile
  cmd+="','../data/out_"
  cmd+=$test_id
  cmd+="/');exit"
  
  mkdir -p ../data/out_$test_id
  srun sh /hb/software/apps/matlab/bin/matlab -nodisplay -nodesktop -nosplash -r $cmd
done
