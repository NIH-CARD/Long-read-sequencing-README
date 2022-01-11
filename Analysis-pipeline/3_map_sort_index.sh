#generate swarm:
#for filename in /data/Neuro_Longread/KOLF_PROCESSING/guppy_out/pass/*.fastq; do echo "bash /data/Neuro_Longread/KOLF_PROCESSING/scripts/3_map_sort_index.sh $filename" /data/Neuro_Longread/KOLF_PROCESSING/split_run >> /data/Neuro_Longread/KOLF_PROCESSING//scripts/map_sort_index.swarm; done

#run swarm
#swarm -g 246 -t 10 --time 4-0 -f /data/Neuro_Longread/KOLF_PROCESSING/scripts/map_sort_index.swarm

#!/bin/bash
FULL_SAMPLE_PATH=$1
BASE_NAME="$(basename "$FULL_SAMPLE_PATH")"

if [[ "$BASE_NAME" == *".gz"* ]]
then 
    SAMPLE_NAME="${BASE_NAME//.fastq.gz}"
else
    SAMPLE_NAME="${BASE_NAME//.fastq}"
fi;




OUT_PATH=$2

echo $FULL_SAMPLE_PATH
echo "Processing ${SAMPLE_NAME}"
echo "Output going to ${OUT_PATH}"

echo "map with winnowmap"
module load winnowmap

winnowmap -W /data/Neuro_Longread/KOLF_PROCESSING/repetitive_k15.txt -ax map-ont /fdb/GENCODE/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa ${FULL_SAMPLE_PATH} > ${OUT_PATH}/${SAMPLE_NAME}.sam


echo "sort and convert to bam"
module load samtools

samtools sort -o ${OUT_PATH}/${SAMPLE_NAME}.bam -@ '$SLURM_CPUS_PER_TASK' -O bam ${OUT_PATH}/${SAMPLE_NAME}.sam


echo "index bam"

samtools index -b -@ '$SLURM_CPUS_PER_TASK' ${OUT_PATH}/${SAMPLE_NAME}.bam