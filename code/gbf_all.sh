#!/bin/bash
#
#SBATCH --job-name=truncgbf
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=100
#SBATCH --time=120:00:00
#SBATCH --output=/hb/home/rngreenw/eternal_inflation/data/gbfall_%j.out
#SBATCH --mail-type=fail
#SBATCH --mail-type=end
#SBATCH --mail-user=rossngreenwood@gmail.com

run_ids = "np.arange(115,1216+1)"
weights = "np.concatenate((np.ones((1,100)))"

ind = 315

cmd = "import numpy as np; from jup_funs import getBinnedFractions2D; fractions=getBinnedFractions2D(run_ids="$run_ids",weights="$weights",binParam=['n_s','r'],bins=(np.linspace(0.95,0.975,30),np.logspace(-3,0,30)),massBounds=((-np.inf,np.inf),(-np.inf,np.inf)),cutParams=['n_s','r'],cutBounds=[[0.95,0.975],[0,0.064]],merit=True); fractions.to_csv('~/eternal_inflation/data/fractions_"$ind".csv'); exit()"
srun python3 -c $cmd
((ind++))

wait
