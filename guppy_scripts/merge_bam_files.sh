#!/bin/bash
# join fastqs
i=0
cwd=$(echo "$PWD")
echo current working directory is "$cwd"
#THIS SCRIPT ASSUMES THAT ALL SAMPLE RESULT FILES ARE STORED IN ONE LEVEL UPSTREAM OF CURRENT FOLDER
#PLEASE UPDATE PATH PROPERLY
declare -a samples=($(ls -d ../*/|cut -d"/" -f2))
i=0
for foldern in "${samples[@]}"
do
        ((i=i+1))
        echo got sample number "$i" and sample "$foldern"
        cd ../"$foldern"
	
	mkdir -p all_pass_sorted_bam
	mkdir -p all_fail_sorted_bam
        for run in run*;
        do
	    runi=$(echo "$run"|rev|cut -d"_" -f 1 | rev)
            echo copying calling all bam files in pass folder in "$run"
            
	    #First sort all bam files in current folder and copy them in all_pass_sorted_bam folder
	    for j in "$run"/pass/*.bam
	    do
        	fn=$(echo $j|rev|cut -d'/' -f1|rev|cut -d'.' -f1)
        	cp $j all_pass_sorted_bam/"$run"_pass_$fn.bam
	    done
	    echo copying all bam files in fail folder in "$run"
	    for j in "$run"/fail/*.bam
            do
                fn=$(echo $j|rev|cut -d'/' -f1|rev|cut -d'.' -f1)
                cp $j all_fail_sorted_bam/"$run"_fail_$fn.bam
            done
        done
	echo done with sorting bam files from all pass and fails folders in all run folders

        samtools merge -t 32 "$foldern"_pass_merged.bam all_pass_sorted_bam/*.bam
	samtools merge -t 32 "$foldern"_fail_merged.bam all_fail_sorted_bam/*.bam

	cd "$cwd"
done


