#!/bin/bash

OUTPUT="$(squeue --partition 256x44 -h)"
args="-c "
if [[ ${#OUTPUT} == "0" ]]; then
  args+="44 "
  args+="$@"
  sbatch --partition 256x44 --cpus-per-task 44 slurm_matlab.sh $args
else
  args+="24 "
  args+="$@"
  sbatch --partition 128x24 --cpus-per-task 24 slurm_matlab.sh $args
fi
echo $args
