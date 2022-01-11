### FAST5 subsetting

# Purpose => subsetting FAST5 files to for example investigate a subset of reads further

# Date => January 11 2022

# Authors => Cornelis

# code from => https://github.com/nanoporetech/ont_fast5_api#fast5_subset 

```
Biowulf

module load ont-fast5-api
# [+] Loading singularity  3.8.5-1  on cn2070 
# [+] Loading ont-fast5-api  3.3.0 

cd /data/Neuro_Longread/

fast5_subset --input /data/multi_reads --save_path /data/subset
    --read_id_list read_id_list.txt --batch_size 100 --recursive

fast5_subset --input /data/Neuro_Longread/Reference_cell_lines/HG_GM24385K_CELL/GM24385K/20211201_1252_1B_PAI79310_8980aa66/fast5/ \
--save_path /data/Neuro_Longread/Reference_cell_lines/HG_GM24385K_CELL_chimeras/ \
--read_id_list /data/Neuro_Longread/Reference_cell_lines/HG_GM24385K_CELL_merged_VS_hifiasm_hg002_pat_v2_1.chimeric_reads.txt \
--batch_size 4000 --recursive -t 10

# format of HG_GM24385K_CELL_merged_VS_hifiasm_hg002_pat_v2_1.chimeric_reads.txt is:
# 00002816-8532-4bab-be65-3924601de45a
# 00003584-2976-489d-9670-45990821f2b8
# etc

# optional: 
tar -zcvf HG_GM24385K_CELL_chimeras.tar.gz HG_GM24385K_CELL_chimeras/


```
