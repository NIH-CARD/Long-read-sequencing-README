#!/bin/bash

# EXAMPLE SCRIPT FOR MAPPING AN UNMAPPED BAM FILE (using minimap2) PREVIOUSLY GENERATED USING GUPPY v6.2 AND PERFORMING VARIANT CALL (pepper_deepvariant) ON THESE MAPPED BAM FILES. PLEASE ADJUST PARAMETERS AS PER YOUR DESIGN

#PLEASE NOTE THAT THIS IS AN EXAMPLE SCRIPT. PLEASE UPDATE ALL PATHS APPROPRIATELY

#Please list Prefix of all BAM files here. In this examples, our BAM files are prefixed with PAH09187 and PAK26016. 
declare -a bams=('PAH09187' 'PAK26016')

batch="umary_2flowcells" #BATCH NAME
BASE="${HOME}/data"
# Set up input data
INPUT_DIR="${BASE}/input"
REF="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" #Referenc Genome
THREADS="48"

# Set up output directory
OUTPUT_DIR="${BASE}/output"
## Create local directory structure
sudo mkdir -p "${OUTPUT_DIR}"
sudo mkdir -p "${INPUT_DIR}"
mkdir -p ../vcall_res/"${batch}"
mkdir -p ../mapped_bams/"${batch}"
#mapping bam files
echo started processing "${#bams[@]}" samples
i=1
for flowcellid in "${bams[@]}"
do
	echo processing sample "$i"
	echo mapping sample "$i"
	#PLEASE Update Path To Your BAM Files and adjust file names accordingly. ALSO ADJUST BAM TAGS AND COMPUTATIONAL RESOURCES AS PER YOUR DESIGN
	samtools fastq  -T Mm,Ml,MM,ML /home/syed/mount/res_adwb/ftx/*"$flowcellid"*pass*.bam | minimap2 -t "$THREADS" -y -x map-ont -a --eqx -k 17 -K 10g GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz - | samtools view -@ "$THREADS" -bh - | samtools sort -@ "$THREADS" - > ../mapped_bams/"$batch"/"$flowcellid".fastq.cpg.hg38.bam
	samtools index ../mapped_bams/"$batch"/"$flowcellid".fastq.cpg.hg38.bam
	
	echo starting vcalling for sample "$i"
	BAM=../mapped_bams/"$batch"/"$flowcellid".fastq.cpg.hg38.bam
	echo copying "$BAM" locally
	cp "$BAM" "${INPUT_DIR}"
	cp "$BAM".bai "${INPUT_DIR}"
	BAM1=$(echo "$BAM"|rev|cut -d/ -f1|rev)
	sudo docker run \
	-v "${INPUT_DIR}":"${INPUT_DIR}" \
	-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
	kishwars/pepper_deepvariant:r0.8 \
	pepper_variant call_variant \
	-b "${INPUT_DIR}/${BAM1}" \
	-f "${INPUT_DIR}/${REF}" \
	-o "${OUTPUT_DIR}" \
	-m "${INPUT_DIR}"/PEPPER_VARIANT_ONT_R941_GUPPY5_SUP_V8.pkl \
	-t ${THREADS} \
	--no_quantized \
	-s "$flowcellid" \
	--ont_r9_guppy5_sup
	
	echo done vcalling, now copying results

	cp -r "${OUTPUT_DIR}" ../vcall_res/"$batch"/"$folder"_"$flowcellid"

	#((i=i+1))

	echo cleaning up

	sudo rm -rf "${OUTPUT_DIR}"
	sudo rm -f "${INPUT_DIR}/$BAM1"
	sudo rm -f "${INPUT_DIR}/$BAM1".bai
	#fi
	((i=i+1))
done
echo done all samples

