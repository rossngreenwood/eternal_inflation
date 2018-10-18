#!/bin/sh
#SBATCH --partition=128x24
#SBATCH --job-name=truncgbf
#SBATCH --ntasks=1
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1000
#SBATCH --output=/hb/home/rngreenw/eternal_inflation/data/truncgbf_%j.out
#SBATCH --mail-type=END
#SBATCH --mail-user=rossngreenwood@gmail.com
sh gen_output_t.sh $1 $2 "$3" "$4"
printf -v id_str '%03d' $7
gbf_cmd="from jup_funs import getBinnedFractions; import numpy as np; fractions = getBinnedFractions(run_ids=np.arange("$1","$2"+1),binParam='"$5"',bins=""$6""); fractions.to_csv('~/eternal_inflation/data/fractions_"$id_str".csv')"
echo $gbf_cmd
python3 -c "$gbf_cmd"
