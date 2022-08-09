#!/bin/bash
#THE ONLY DIFFERENCE (from basecalling only script) IS THE USE OF --bam_out FLAG. ALSO WE ARE NOT USING ANY REFERENCE GENOME HERE BUT YOU CAN ADD IT IF REQUIRED
#POPULATE SAMPLE FOLDERS LIST - Please make sure that sample directories are properly captured in "samples" array. You might need to change "cut" command

declare -a samples=($(ls -d */|cut -d"/" -f1))
#POPULATE MIG ADDRESSES (Here we have used 4 A100 GPUs in MIG2 mode and hence we have 8 MIG address). PLEASE SEE README.md to get MIG address.
declare -a gpus=('MIG-0a87138a-788c-55b6-9745-11fdd2228873' 'MIG-eab0e6aa-4c44-5f29-bdb3-d1d7441c42ab' 'MIG-e42c3298-6dec-5ce0-8539-70512752db55' 'MIG-7fd12259-8ede-5c46-84fd-448cc5ef8d30' 'MIG-2acca5f2-2aa7-5428-9c68-1436863244a1' 'MIG-5d458a56-24a7-5cdb-8fc1-fc9eb20a1b6c' 'MIG-4cbb3e78-74a0-5258-9b20-04ed306ad50d' 'MIG-308d1d58-74bc-53c4-905f-8b4318925861')
i=0
for foldern in "${samples[@]}"
do
	t1=`date +%s` #keep track of total time used by each sample
	mkdir -p "$foldern"
	((i=i+1))
	
	#Please adjust this path if your samples are located at different location
	fn=$(ls ../mountfolder/long/FAST5/$foldern/*.fast5|wc -l) #total number of fast5 files in current sample
	fn1=$((fn/8+1)) #Here we submit 8 (num of MIGs/GPU * num of GPUS) guppy runs, Please change accordingly
	
	echo now starting sample number $i in $foldern with total fast5 files $fn and splitting into two batches of $fn1 files
	
	mkdir -p "$foldern"
	gsutil -q -m rsync -r gs://foundin_long_short/long/FAST5/"$foldern" "$foldern"
	d=1
        mkdir -p run$d
        cd run$d
	j=0
        for ii in ../"$foldern"/*.fast5
        do
                ln -s "$ii"
		((j=j+1))
		if [ $(expr $j % $fn1) == "0" ] #copy upto 1000 files
                then
			((d=d+1))
                        mkdir -p ../run$d
                        cd ../run$d
                fi

        done
        cd ../

	d=1
	
	#PLEASE NOYE THAT SETTINGS USED FOR GUPPY ARE OPTIMIZED FOR A100 GPU, 40 GB WORKING IN MIG2 MODE. YOU CHANGE THESE IF YOU ARE
USING DIFFERENT GPU/MIG mode	
	mkdir -p res_"$foldern"
	for mig in "${gpus[@]}"
	do
	{
                echo processing sample $foldern for folder run$d
                OUTP=res_"$foldern"/run$d
                CUDA_VISIBLE_DEVICES=$mig ~/guppy_basecaller  --chunks_per_runner 768 --disable_pings --compress_fastq -i run$d  -s ${OUTP} -c dna_r9.4.1_450bps_modbases_5mc_cg_sup_prom.cfg -x auto -r --read_batch_size 250000 -q 25000 --bam_out  > res_"$foldern"/dout$d.txt 
		echo finished run$d
	} &
	
	((d=d+1))
	done


	wait

	echo done all calculations for sample "$foldern", now copying results from res_"$foldern" folder to ../mountfolder/all_bam_res/res_"$foldern" folder
        mkdir -p ../mountfolder/all_bam_res/res_"$foldern"

        gsutil -q -m rsync -r res_"$foldern" gs://foundin_long_short/all_bam_res/res_"$foldern"
	for d in {1..8} #Please change this if you are not using 4 GPUs with 2MIGs/GPU
	do
	mkdir -p ../mountfolder/all_bam_res/res_"$foldern"/run$d
        gsutil -q -m rsync -r res_"$foldern"/run$d gs://foundin_long_short/all_bam_res/res_"$foldern"/run$d

	 mkdir -p ../mountfolder/all_bam_res/res_"$foldern"/run$d/pass
        gsutil -q -m rsync -r res_"$foldern"/run$d/pass gs://foundin_long_short/all_bam_res/res_"$foldern"/run$d/pass


	 mkdir -p ../mountfolder/all_bam_res/res_"$foldern"/run$d/fail
        gsutil -q -m rsync -r res_"$foldern"/run$d/fail gs://foundin_long_short/all_bam_res/res_"$foldern"/run$d/fail

	rm -r run$d
	done


	echo copied data from res_"$foldern" folder for sample "$foldern" to gs://foundin_long_short/all_bam_res/res_"$foldern"
	echo NOW deleting all files and preparing for next sample
	rm -r "$foldern"
	rm -r res_"$foldern"

	t2=`date +%s`

	t=$((t2-t1))

	echo TOTAL TIME FOR SAMPLE "$foldern" is: $(($t/3600)) Hours

done

#THIS makes sure that VM is shutdown when all jobs are finished.

sudo shutdown
