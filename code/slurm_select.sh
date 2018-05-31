#!/bin/bash

OUTPUT="$(squeue --partition 256x44 -h)"
args="-c "
if [[ ${#OUTPUT} == "0" ]]; then
  args+="44 "
  args+="$@"
  comment+="$@"
  sbatch --partition 256x44 --cpus-per-task 44 /hb/home/rngreenw/eternal_inflation/code/slurm_matlab2.sh $args
else
  args+="24 "
  args+="$@"
  comment+="\""
  comment+="$@"
  comment+="\""
  sbatch --partition 128x24 --cpus-per-task 24 /hb/home/rngreenw/eternal_inflation/code/slurm_matlab2.sh $args
fi
