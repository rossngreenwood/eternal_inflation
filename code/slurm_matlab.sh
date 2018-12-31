#!/bin/bash
#
#SBATCH --job-name=etrnlinf
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1000
#SBATCH --time=120:00:00
#SBATCH --output=/hb/home/rngreenw/eternal_inflation/data/etrnlinf_%j.out
#SBATCH --mail-type=fail
#SBATCH --mail-type=end
#SBATCH --mail-user=rossngreenwood@gmail.com

# Defaults
cores="24"
range_flag=0

while getopts ":r::c:" opt; do
  echo "${OPTARG}"
  case "${opt}" in
    r)
      range_flag=1
      ;;
    c)
      cores="${OPTARG}"
      ;;
    *)
      echo "Invalid option: -$OPTARG" # >&2
      ;;
  esac
done

if [[ $range_flag -gt 0 ]]; then
  
  # Make sure output directories exist
  for ((test_id=${@:(-2):1}; test_id<=${@: -1}; test_id++))
  do
    printf -v id_str '%04d' $test_id
    mkdir -p /hb/home/rngreenw/eternal_inflation/data/out_$id_str
  done
  
  infile="/hb/home/rngreenw/eternal_inflation/data/input/infile_%04d.txt"
  outdir="/hb/home/rngreenw/eternal_inflation/data/out_%04d/"
  
  printf -v range_str '[%d,%d]' ${@:(-2):1} ${@: -1}  
  
  cmd="cd('/hb/home/rngreenw/eternal_inflation/code');eternal_sim_leader("
  cmd+=$cores
  cmd+=",true,"
  cmd+=$range_str
  cmd+=",'"
  cmd+=$infile
  cmd+="','"
  cmd+=$outdir
  cmd+="');exit"
  
  # Run simulations for all test_ids
  srun sh /hb/software/apps/matlab/2017b/bin/matlab -nodisplay -nodesktop -r $cmd
  
  # Combine output
  #for ((test_id=${@:(-2):1}; test_id<=${@: -1}; test_id++))
  #do
  #  
  #  printf -v id_str '%04d' $test_id
  #  
  #  outfile="outfile_"
  #  outfile+=$id_str
  #  outfile+=".txt"
  #  
  #  outfile_t="outfile_t_"
  #  outfile_t+=$id_str
  #  outfile_t+=".txt"
  #  
  #  outdir="/hb/home/rngreenw/eternal_inflation/data/out_"
  #  outdir+=$id_str
  #  outdir+="/"
  #  
  #  python /hb/home/rngreenw/eternal_inflation/code/eternal_sim_cleanup.py --cores $cores --output_dir $outdir --output_file $outfile
  #  python /hb/home/rngreenw/eternal_inflation/code/eternal_sim_cleanup.py --cores $cores --output_dir $outdir --output_file $outfile_t --truncate 1
  #
  #done

else

  for test_id in "$@"
  do
    
    printf -v id_str '%04d' $test_id
    
    # Make sure output directory exists
    mkdir -p /hb/home/rngreenw/eternal_inflation/data/out_$id_str
    
    infile="/hb/home/rngreenw/eternal_inflation/data/input/infile_"
    infile+=$id_str
    infile+=".txt"
    
    outfile="outfile_"
    outfile+=$id_str
    outfile+=".txt"
    
    outfile_t="outfile_t_"
    outfile_t+=$id_str
    outfile_t+=".txt"
    
    outdir="/hb/home/rngreenw/eternal_inflation/data/out_"
    outdir+=$id_str
    outdir+="/"
    
    cmd="cd('/hb/home/rngreenw/eternal_inflation/code');eternal_sim_leader("
    cmd+=$cores
    cmd+=",false,"
    cmd+=$id_str
    cmd+=",'"
    cmd+=$infile
    cmd+="','"
    cmd+=$outdir
    cmd+="');exit"
    
    # Run simulation for a single test_id
    srun sh /hb/software/apps/matlab/2017b/bin/matlab -nodisplay -nodesktop -r $cmd
    
    # Combine output
    python /hb/home/rngreenw/eternal_inflation/code/eternal_sim_cleanup.py --cores $cores --output_dir $outdir --output_file $outfile
    python /hb/home/rngreenw/eternal_inflation/code/eternal_sim_cleanup.py --cores $cores --output_dir $outdir --output_file $outfile_t --truncate 1
  
  done

fi
