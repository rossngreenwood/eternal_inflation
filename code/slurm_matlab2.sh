#!/bin/bash
#
#SBATCH --job-name=etrnlinf
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1000
#SBATCH --time=48:00:00
#SBATCH --output=/hb/home/rngreenw/eternal_inflation/data/etrnlinf_%j.out
#SBATCH --mail-type=begin
#SBATCH --mail-type=fail
#SBATCH --mail-type=end
#SBATCH --mail-user=rossngreenwood@gmail.com

echo $SLURM_TASK_PID

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
  
  infile="/hb/home/rngreenw/eternal_inflation/data/input/infile_"
  infile+="%04d"
  infile+=".txt"
  
  outdir="/hb/home/rngreenw/eternal_inflation/data/out_"
  outdir+="%04d"
  outdir+="/"
  
  printf -v range_str '[%d,%d]' ${@:(-2):1} ${@: -1}  
  
  cmd="cd('/hb/home/rngreenw/eternal_inflation/code');eternal_sim_leader("
  cmd+=$cores
  cmd+=",1,"
  cmd+=$range_str
  cmd+=",'"
  cmd+=$infile
  cmd+="','','"
  cmd+=$outdir
  cmd+="');exit"
  
  # Make sure output directories exist
  for ((test_id=${@:(-2):1}; test_id<=${@: -1}; test_id++))
  do
    printf -v id_str '%04d' $test_id
    mkdir -p /hb/home/rngreenw/eternal_inflation/data/out_$id_str
  done
  
  # Run simulations for all test_ids
  srun sh /hb/software/apps/matlab/2017b/bin/matlab -nodisplay -nodesktop -r $cmd
  
  # Combine output
  for ((test_id=${@:(-2):1}; test_id<=${@: -1}; test_id++))
  do
    
    printf -v id_str '%04d' $test_id
    
    outfile="outfile_"
    outfile+=$id_str
    outfile+=".txt"
    
    outfile_t="outfile_t_"
    outfile_t+=$id_str
    outfile_t+=".txt"
    
    outdir="/hb/home/rngreenw/eternal_inflation/data/out_"
    outdir+=$id_str
    outdir+="/"
    
    python /hb/home/rngreenw/eternal_inflation/code/eternal_sim_cleanup.py --cores $cores --output_dir $outdir --output_file $outfile
    python /hb/home/rngreenw/eternal_inflation/code/eternal_sim_cleanup.py --cores $cores --output_dir $outdir --output_file $outfile_t --truncate 1
  done

else

  for test_id in "$@"
  do
    infile="/hb/home/rngreenw/eternal_inflation/data/input/infile_"
    infile+=$test_id
    infile+=".txt"
    
    outfile="outfile_"
    outfile+=$test_id
    outfile+=".txt"
    
    outfile_t="outfile_t_"
    outfile_t+=$id_str
    outfile_t+=".txt"
    
    outdir="/hb/home/rngreenw/eternal_inflation/data/out_"
    outdir+=$test_id
    outdir+="/"
    
    cmd="cd('/hb/home/rngreenw/eternal_inflation/code');eternal_sim_leader("
    cmd+=$cores
    cmd+=",'"
    cmd+=$infile
    cmd+="','"
    cmd+=$outfile
    cmd+="','"
    cmd+=$outdir
    cmd+="');exit"
    
    mkdir -p /hb/home/rngreenw/eternal_inflation/data/out_$test_id
    srun sh /hb/software/apps/matlab/2017b/bin/matlab -r $cmd
    python /hb/home/rngreenw/eternal_inflation/code/eternal_sim_cleanup.py --cores $cores --output_dir $outdir --output_file $outfile
    python /hb/home/rngreenw/eternal_inflation/code/eternal_sim_cleanup.py --cores $cores --output_dir $outdir --output_file $outfile_t --truncate 1
  done

fi
