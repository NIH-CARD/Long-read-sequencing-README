#!/bin/bash
#  sbatch --mem=246g --cpus-per-task=10 --mail-type=ALL --time=4-0 --gres=lscratch:500  2_run_meryl.sh
module load winnowmap
meryl count k=15 output merylDB /fdb/GENCODE/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
