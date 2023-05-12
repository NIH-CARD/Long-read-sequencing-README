#!/bin/bash
#AN EXAMPLES SCRIPT FOR SeqKit
#PLEASE POPULATE SAMPLES
declare -a samples=('NABEC_KEN-1095_FTX_PAM08252_pass_joined.fastq.gz' 'NABEC_KEN-1127_FTX_PAM08415_pass_joined.fastq.gz' 'NABEC_KEN-1131_FTX_PAM30772_pass_joined.fastq.gz' 'NABEC_KEN-1142_FTX_PAM32160_pass_joined.fastq.gz' 'NABEC_KEN-1159_FTX_PAM35004_pass_joined.fastq.gz')

mkdir -p seqkit_all_res
i=1
for fastq in "${samples[@]}"
do
	sample=$(echo "$fastq"|rev|cut -d_ -f3-|rev)
	echo "$i","$sample"
	seqkit stats -j 24 -a ../../res_adwb/nabec/"$fastq" -T > seqkit_all_res/"$sample".stats.pass.tsv
	((i=i+1))
done