#!/bin/bash

ref_dire=$(realpath $1)
sumstat=$(realpath $2)
gwas=$(realpath $3)
out_prefix=$(realpath $4)
bimfile=$(realpath $5)
chr=$(echo $6)
code=$(realpath $7)



sample_size=$(awk -F "," 'NR>1{print $6}' $gwas|awk '{ if ($1 > max) max = $1 } END { print max }')


N_THREADS=1
export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS


python --version
python $code --ref_dir=${ref_dire} --bim_prefix=${bimfile} --sst_file=${sumstat} --n_gwas=$sample_size --out_dir=${out_prefix} --chrom=${chr}

