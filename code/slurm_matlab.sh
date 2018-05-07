#!/bin/bash
#
#SBATCH --partition=128x24
#SBATCH --job-name=etrnlinf
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=1000
#SBATCH --time=12:00:00
#SBATCH --output=../data/etrnlinf_%j.out
#SBATCH --mail-type=begin
#SBATCH --mail-type=fail
#SBATCH --mail-type=end
#SBATCH --mail-user=rossngreenwood@gmail.com

echo $SLURM_TASK_PID

range_flag=0
while getopts ":r" opt; do
  case $opt in
    r)
      range_flag=1
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
  esac
done

if [[ $range_flag -gt 0 ]]; then

  for ((test_id=$2; test_id<=$3; test_id++))
  do
    printf -v id_str '%04d' $test_id
    
    infile="../data/input/infile_"
    infile+=$id_str
    infile+=".txt"
    
    outfile="outfile_"
    outfile+=$id_str
    outfile+=".txt"
    
    outfile_t="outfile_t_"
    outfile_t+=$id_str
    outfile_t+=".txt"
    
    outdir="../data/out_"
    outdir+=$id_str
    outdir+="/"
    
    cmd="eternal_sim_leader("
    cmd+="24,'"
    cmd+=$infile
    cmd+="','"
    cmd+=$outfile
    cmd+="','"
    cmd+=$outdir
    cmd+="');exit"
    
    mkdir -p ../data/out_$id_str
    srun sh /hb/software/apps/matlab/bin/matlab -nodisplay -nodesktop -r $cmd
    python eternal_sim_cleanup.py --cores 24 --output_dir $outdir --output_file $outfile
    python eternal_sim_truncate.py --cores 24 --output_dir $outdir --output_file $outfile_t
  done

else

  for test_id in "$@"
  do
    infile="../data/input/infile_"
    infile+=$test_id
    infile+=".txt"
    
    outfile="outfile_"
    outfile+=$test_id
    outfile+=".txt"
    
    outfile_t="outfile_t_"
    outfile_t+=$id_str
    outfile_t+=".txt"
    
    outdir="../data/out_"
    outdir+=$test_id
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
    python eternal_sim_truncate.py --cores 24 --output_dir $outdir --output_file $outfile_t
  done

fi
