#!/bin/bash
module load guppy/6.0.1
FAST5_PATH=$1
OUT_PATH=$2

guppy_basecaller -i ${FAST5_PATH} -s ${OUT_PATH} -c dna_r9.4.1_450bps_sup_prom.cfg -x cuda:all:100% -r --read_batch_size 250000 -q 250000

#sbatch --partition=gpu --cpus-per-task=14 --mem=10g --gres=gpu:v100x:2,lscratch:200 --time=4-0 --wrap="bash 1_run_guppy_basecaller.sh /fast5_path /out_path"
