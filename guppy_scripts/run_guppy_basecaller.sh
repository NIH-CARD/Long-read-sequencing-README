#!/bin/bash
#Please populate samples to process
declare -a samples=('HBCC_81948_FTX/PAK94100' 'HBCC_81986_FTX/PAK95500' 'HBCC_81986_FTX/PAK95790' 'HBCC_81987_FTX/PAK64889' 'HBCC_81987_FTX/PAK95613' 'HBCC_81988_FTX/PAK57535' 'HBCC_81988_FTX/PAK64818' 'HBCC_81989_FTX/PAK63420' 'HBCC_81989_FTX/PAK95824' 'HBCC_81997_FTX/PAK69358' 'HBCC_81997_FTX/PAK96150' 'HBCC_81998_FTX/PAK63279')
#Please populate MIG hardware addresses here. Detailed description on how to setup A100 GPUs in MIG and create Compute instances, please see here (https://codeyarns.com/tech/2020-12-15-how-to-use-mig.html)
declare -a gpus=('MIG-36490508-d664-5b42-a1a2-c28b6c98cd9e' 'MIG-a45148db-8e85-573d-ac7c-a4081ab6820f' 'MIG-17ad86db-6a10-5e66-9634-087e7bcd7b04' 'MIG-66716305-c5dd-5239-a148-dedcfe6117db' 'MIG-2e37af9d-a8a4-5d77-b19a-56da02211325' 'MIG-de9e7900-2c5b-5103-a0fa-c71701174df0' 'MIG-2efa367a-27ee-5d5a-88e1-bcc8ae1e3fb0' 'MIG-4a046b60-5b4b-5614-98be-125133bd76a6')


nmigs=8
i=1
ts=`date +%s`
for foldern in "${samples[@]}"
do
    
    echo copying /niadatastore_data/"$foldern" locally

	mkdir -p "$foldern"
	t11=`date +%s` #keep track of total time used by each sample
	#echo Copy Start time is "$t11"

	cp /niadatastore_data/"$foldern"/*.fast5 "$foldern"
	#AS LR-Seq samples are often in TB, we first copy sample on local ssd for optimized throughput
	echo done copying /niadatastore_data/"$foldern" folder locally
	sz=$(du -sh "$foldern")
	echo folder size is: $sz
	t22=`date +%s`
	tt=$((t22-t11))
	echo TOTAL TIME TO COPY SAMPLE "$foldern" is: $tt
	t1=`date +%s` #keep track of total time used by each sample
	fn=$(ls "$foldern"/*.fast5|wc -l) #total number of fast5 files in current sample
	fn1=$((fn/8)) #Here we submit 8 (num of MIGs/GPU * num of GPUS) guppy runs, Please change accordingly

	echo starting sample $i with total fast5 files $fn and splitting into 8 batches of $fn1 files
	d=1
    mkdir -p run$d
	cd run$d
	j=0
	for ii in /datastore-8/nguppy/"$foldern"/*.fast5
	do
		ln -s "$ii"
		((j=j+1))
		if [ $(expr $j % $fn1) == "0" ] #copy upto 1000 files
		then
			((d=d+1))
			if [ $d -le $nmigs ] #last folder would get some remainder firles rem=$((fn%8))
            then
    			mkdir -p ../run$d
    			cd ../run$d
    		fi
		fi
	done
	cd ../
	d=1
	mkdir -p res_"$foldern"
	for mig in "${gpus[@]}"
	do
	{
		echo processing sample $foldern for folder run$d
		OUTP=res_"$foldern"/run$d
		CUDA_VISIBLE_DEVICES=$mig ~/guppy_basecaller  --chunks_per_runner 768 --disable_pings --compress_fastq -i run$d  -s ${OUTP} -c dna_r9.4.1_450bps_modbases_5mc_cg_sup_prom.cfg -x auto -r --read_batch_size 250000 -q 25000  > res_"$foldern"/dout$d.txt 
		echo finished run$d
	} &
	((d=d+1))
	done
	wait
	
	echo done sample "$foldern", copying res_"$foldern" to /niadatastore_output/new_bam_results/res_"$foldern" folder
	t2=`date +%s`
	#echo Guppy End time is "$t2"
	t=$((t2-t1))
	echo TOTAL TIME FOR SAMPLE "$foldern" is: $t
        mkdir -p /niadatastore_output/new_bam_results/res_"$foldern"
    t111=`date +%s`
    #echo Results Copy Start time is "$t111"

    rsync -r res_"$foldern"/ /niadatastore_output/new_bam_results/res_"$foldern"/
    
    t222=`date +%s`
    #echo Results Copy End time is "$t222"
    tt=$((t222-t111))
    echo Total Results Copy time is $tt
    echo NOW deleting all files and preparing for next sample
    rm dout*
    rm -r run*
	rm -r "$foldern"
    rm -r res_"$foldern"
    ((i=i+1))
done
te=`date +%s`

tall=$((te-ts))
echo All samples time is: $tall

echo DONE All samples
sudo shutdown #shutdown after all samples are done


