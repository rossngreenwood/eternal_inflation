#!/bin/bash

for ((test_id = $1; test_id <= $2; test_id++))
do
  
  printf -v id_str '%04d' $test_id
  
  echo $id_str 
  
  outdir="../data/out_$id_str/"
  outfile="outfile_t_$id_str.txt"
  
  python eternal_sim_cleanup.py --output_dir $outdir --output_file $outfile --truncate 1

done
