## SAMPLOT

```
## https://github.com/ryanlayer/samplot

## setup samplot
module load python
source /data/'blauwendraatc'/conda/etc/profile.d/conda.sh
conda activate base
conda install -c bioconda samplot 

## downloading gene tracks
wget ftp://ftp.ensembl.org/pub/release-105/gff3/homo_sapiens/Homo_sapiens.GRCh38.105.chr.gff3.gz 
bedtools sort -i Homo_sapiens.GRCh38.105.chr.gff3.gz  \
| bgzip -c > Homo_sapiens.GRCh38.105.chr.sort.gff3.gz 
tabix Homo_sapiens.GRCh38.105.chr.sort.gff3.gz

## run samplot
samplot plot \
    -n KOLF_N50_30kb \
    -b split_run_merged_bams_sorted.bam \
    -o chr6_kolf_15496930_15722102.png \
    -c chr6 \
    -s 15496930 \
    -e 15722102 \
    -t DEL \
    -d 100 \
    -T Homo_sapiens.GRCh38.105.chr.sort.gff3.gz

samplot plot \
    -n KOLF_N50_30kb \
    -b split_run_merged_bams_sorted.bam \
    -o chr18_kolf_62115075_62262535.png \
    -c chr18 \
    -s 62115075 \
    -e 62262535 \
    -t DUP \
    -d 100 \
    -T Homo_sapiens.GRCh38.105.chr.sort.gff3.gz


### illumina
cd /data/Neuro_Longread/KOLF_cell_line/KOLF_Illumina_data/cram_NIH/

# first convert cram to bam
samtools view -b  -T /fdb/GENCODE/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa -o 22035_6%231.bam 22035_6%231.cram
samtools index 22035_6%231.bam


samplot plot \
    -n KOLF_NIH_Illumina \
    -b 22035_6%231.bam \
    -o /data/Neuro_Longread/KOLF_PROCESSING/chr6_kolf_15496930_15722102_illumina.png \
    -c chr6 \
    -s 15496930 \
    -e 15722102 \
    -t DUP \
    -d 100 \
    -T /data/Neuro_Longread/KOLF_PROCESSING/Homo_sapiens.GRCh38.105.chr.sort.gff3.gz

```
