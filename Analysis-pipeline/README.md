
Deposit space for general scripts used locally for long-read projects:

Includes:

FAST5_subset.sh => this subsets FAST5 files based on read-names

1_run_guppy_basecallersh => basecalling with guppy using fast5s

2_run_meryl.sh => Identify and count the top 0.02% most frequent 15-mers in the hg38 reference genome prior to mapping with winnowmap.

3_map_sort_index.sh => Map long-reads to reference genome with winnowmap. This also converts the sam to a bam file, sorts it and indexes it. 

split_on_adapter_adjust.swarm => example swarm file for splitting chimeric reads in fastq files.

