#!/bin/bash
#
#SBATCH --job-name=truncgbf
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=100
#SBATCH --time=8:00:00
#SBATCH --output=/hb/home/rngreenw/eternal_inflation/data/truncgbf_%j.out
#SBATCH --mail-type=fail
#SBATCH --mail-type=end
#SBATCH --mail-user=rossngreenwood@gmail.com

run_ids =
weights =

cmd = "import numpy as np; import pandas as pd; from jup_funs import massBin2D; massBin2D(run_ids="$run_ids",weights="$weights",cutparam=['n_s','r'],cutbounds=[[0.95,0.975],[0,0.064]],outfile='fractions_310.csv'); exit()"
srun python3 -c $cmd

wait
