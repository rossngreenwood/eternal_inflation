#!/bin/bash
# Generate truncated output
sh gen_output_t.sh $1 $2 "$3" "$4"
# Get binned fractions
gbf_cmd="from jup_funs import getBinnedFractions; import numpy as np; fractions = getBinnedFractions(run_ids=np.arange("$1","$2"+1),binParam='"$5"',bins=""$6""); fractions.to_csv('fractions_tmp.csv')"
python3 -c "$gbf_cmd"
