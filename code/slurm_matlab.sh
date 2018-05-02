#!/bin/bash
#
#SBATCH --partition=128x24
#SBATCH --job-name=etrnlinf
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=1000
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
  
  outdir="../data/out_"
  outdir+=$testid
  outdir+="/"
  
  cmd="eternal_sim_leader("
  cmd+="24,'"
  cmd+=$infile
  cmd+="','"
  cmd+=$outfile
  cmd+="','"
  cmd+=$outdir
  cmd+="');exit"
  
  mkdir -p ../data/out_$test_id
  srun sh /hb/software/apps/matlab/bin/matlab -r $cmd
  python eternal_sim_cleanup.py --cores 24 --output_dir $outdir --output_file $outfile
done
