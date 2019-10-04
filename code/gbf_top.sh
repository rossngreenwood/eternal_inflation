#!/bin/bash
#
#SBATCH --job-name=truncgbf
#SBATCH --ntasks=5
#SBATCH --mem-per-cpu=100
#SBATCH --time=8:00:00
#SBATCH --output=/hb/home/rngreenw/eternal_inflation/data/truncgbf_%j.out
#SBATCH --mail-type=fail
#SBATCH --mail-type=end
#SBATCH --mail-user=rossngreenwood@gmail.com

run_ids = "np.arange(115,1216+1)"
weights = "np.ones((1,1216-115+1)"

sh gen_output_t.sh 115 1216 "[record_flag]" "[[3]]" 

ind = 315

cmd = "import numpy as np; from jup_funs import massBin2D; fractions=massBin2D(run_ids="$run_ids"); fractions.to_csv('~/eternal_inflation/data/fractions_"$ind".csv'); exit()"
srun python3 -c $cmd
((ind++))

cmd = "import numpy as np; from jup_funs import massBin2D; fractions=massBin2D(run_ids="$run_ids",cutParams=['n_s'],cutBounds=[[0.95,0.975]]); fractions.to_csv('~/eternal_inflation/data/fractions_"$ind".csv'); exit()"
srun python3 -c $cmd
((ind++))

cmd = "import numpy as np; from jup_funs import massBin2D; fractions=massBin2D(run_ids="$run_ids",cutParams=['r'],cutBounds=[[0,0.064]]); fractions.to_csv('~/eternal_inflation/data/fractions_"$ind".csv'); exit()"
srun python3 -c $cmd
((ind++))

cmd = "import numpy as np; from jup_funs import massBin2D; fractions=massBin2D(run_ids="$run_ids",cutParams=['n_s','r'],cutBounds=[[0.95,0.975],[0,0.064]]); fractions.to_csv('~/eternal_inflation/data/fractions_"$ind".csv'); exit()"
srun python3 -c $cmd
((ind++))

cmd = "import numpy as np; from jup_funs import massBin2D; fractions=massBin2D(run_ids="$run_ids",cutParams=['Q'],cutBounds=[[np.power(10.,-4.5),np.power(10.,-4.1)]]); fractions.to_csv('~/eternal_inflation/data/fractions_"$ind".csv'); exit()"
srun python3 -c $cmd
((ind++))

wait
