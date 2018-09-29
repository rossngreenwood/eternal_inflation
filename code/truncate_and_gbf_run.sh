#!/bin/bash
#
#SBATCH --job-name=etrnlout
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1000
#SBATCH --cpus-per-task=4
#SBATCH --time=3:00:00
#SBATCH --output=/hb/home/rngreenw/eternal_inflation/data/etrnlout_%j.out
#SBATCH --mail-type=fail
#SBATCH --mail-type=end
#SBATCH --mail-user=rossngreenwood@gmail.com

ssh rngreenw@hb.ucsc.edu 'cd /hb/home/rngreenw/eternal_inflation/code/; srun --job-name=etrnlout --ntasks=1 --mem-per-task=4 --time=3:00:00 --output=/hb/home/rngreenw/eternal_inflation/data/etrnlout_%j.out sh truncate_and_gbf.sh $1 $2 "$3" "$4" $5 "$6"
scp rngreenw@hb.ucsc.edu:~/eternal_inflation/code/fractions_tmp.csv ./
