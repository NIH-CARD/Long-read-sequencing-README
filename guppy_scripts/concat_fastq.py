# join fastqs
import subprocess
import os
import sys
def join_fastq(path, out):
    import gzip
    files = os.path.join(path,'*.fastq.gz')
    cmd = 'cat {f} > {o}'.format(f=files,o=out)
    print (cmd)
    subprocess.check_output(cmd,shell=True)
    cmd = 'gzip %s' %out
    #NO NEED TO COMPRESS AS WE ARE CONCATENATING COMPRESSED FASTQ
    #print ('compressing')
    #subprocess.check_output(cmd,shell=True)
    return

join_fastq(sys.argv[1],sys.argv[2])
#join_fastq('../benchmark_res/fastq_RNAL_PPMI3452_7426_da65_v1/pass/', '../benchmark_res/RNA_sup_joined.fastq.gz')
