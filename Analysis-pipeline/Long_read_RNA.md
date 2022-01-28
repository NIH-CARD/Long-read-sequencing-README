## Starting some notes on rough overview of pipeline of long-read RNA data


## Basecalling and QC 
```
FAST5 Basecalling and filtering
Guppy (v4.4.1) documentation
guppy_basecaller --compress_fastq -i /fast5 --save_path /fastq/ --gpu_runners_per_device 30 \
-c dna_r9.4.1_450bps_hac_prom.cfg -x "cuda:0" --qscore_filtering --min_qscore 7 --trim_strategy none

Quality control 
NanoPlot (v1.32.1)

## Primer and adapter trimming, read reorientation
# Pychopper (v2.5.0) documentation 
# https://github.com/nanoporetech/pychopper
cdna_classifier.py -t 30 -x PCS109 -r report.pdf -u unclassified.fq -w rescued.fq classified.fastq input_q7.fastq

```

## Mapping
```
## Minimap 2
# https://github.com/lh3/minimap2
./minimap2 -ax splice:hq -uf ref.fa query.fa > aln.sam    # Final PacBio Iso-seq or traditional cDNA


```


## Gene expression analysis

```
Gene and transcript quantification, identification and annotation pipelines:

## FLAIR (v1.5)  documentation
# https://github.com/BrooksLabUCSC/flair
python flair.py 12346 -r pychopper_output.fq -g hg38.fa -f gencode_v29.lncipedia_v5_2_hc.annotation.sorted.gt -o flair.output --temp_dir temp_flair

## BAMBU (v0.3.0) documentation
# https://github.com/GoekeLab/bambu
run_bambu.r --tag=bambu --ncore=10 --annotation=gencode_v29.lncipedia_v5_2_hc.annotation.sorted.gtf --fasta=hg38.fa flair_output.bam

```
