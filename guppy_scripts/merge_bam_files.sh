#!/bin/bash
#EXAMPLE SCRIPT TO MERGE ALL BAM AND FASTQ FILES FOR A SAMPLE BASECALLED ON MULTIPLE MIG Machinces on A100 GPU in MIG MODE
#PLEASE NOTE THAT IN OUR CASE, WE SPLIT SAMPLES IN 8 DIFFERENT FOLDERS WHERE FILES IN EACH FOLDER ARE PROCESSED BY A MIG MACHINE, SO WE NEED TO COMBINE THESER RESULTS.
i=1
declare -a samples=('HBCC_82018_FTX/PAK62963' 'HBCC_82020_FTX/PAK62323' 'HBCC_82020_FTX/PAK70073' 'HBCC_82021_FTX/PAK63241' 'HBCC_82021_FTX/PAK69826' 'HBCC_82022_FTX/PAK62188' 'HBCC_82022_FTX/PAK73635' 'HBCC_81948_FTX/PAK94100' 'HBCC_81986_FTX/PAK95500' 'HBCC_81986_FTX/PAK95790' 'HBCC_81987_FTX/PAK64889' 'HBCC_81987_FTX/PAK95613' 'HBCC_81988_FTX/PAK57535' 'HBCC_81988_FTX/PAK64818' 'HBCC_81998_FTX/PAK96599' 'HBCC_81999_FTX/PAK95499' 'HBCC_81999_FTX/PAK96387' 'HBCC_82005_FTX/PAK70222' 'HBCC_82005_FTX/PAK73463' 'HBCC_82015_FTX/PAK70080' 'HBCC_82015_FTX/PAK70256' 'HBCC_82016_FTX/PAK70245')

cwd=$(echo "$PWD")
echo current working directory is "$cwd"
echo started post-processing "${#samples[@]}" samples
for foldern in "${samples[@]}"
do
    cellid=$(ls -l ../"$foldern"/*.fast5|head -1|rev|cut -d/ -f1|rev|cut -d_ -f1)
    sample=$(echo "$foldern" | cut -d'/' -f2)
    res_folder=../new_bam_results/res_"$foldern"
    echo got sample number "$i", sample "$foldern" and flowcell $cellid
    
    cd ../new_bam_results/res_"$foldern"
    #concat all fastq files
    echo merging fastq files for sample ../new_bam_results/res_"$foldern"
    cat run*/pass/*.gz > "$sample"_"$cellid"_pass_joined.fastq.gz
    cat run*/fail/*.gz > "$sample"_"$cellid"_fail_joined.fastq.gz
    echo merging bam files for sample /niadatastore_output/new_bam_results/res_"$foldern"
    samtools merge -f -t 8 "$sample"_"$cellid"_pass_merged.bam run*/pass/*.bam
    samtools merge -f -t 8 "$sample"_"$cellid"_fail_merged.bam run*/fail/*.bam

    echo now merging summary files for sample "$foldern"
    awk 'NR == 1 { print }; FNR > 1 { print }' run*/sequencing_summary* > "$sample"_"$cellid"_sequencing_summary.txt

    cd "$cwd"    
    ((i=i+1))
done

sudo shutdown
