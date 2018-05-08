#!/bin/bash

for ((test_id = $1; test_id <= $2; test_id++))
do
  
  printf -v id_str '%04d' $test_id
  
  outdir="../data/out_"
  outdir+=$id_str
  outdir+="/"
  
  outfile="outfile_t_"
  outfile+=$id_str
  outfile+=".txt"
  
  echo $id_str
  
  python eternal_sim_truncate.py --cores 24 --output_dir $outdir --output_file $outfile
done
