#!/bin/sh

#############################################################
##  Run get_coveragy.py gets covereage on a trinity output ##
#############################################################

##This will download and process two SRA files from Nicotiana bonariensis.
sra_ids="SRR6918807 SRR6516146" 
for sra_id in $sra_ids; do
	python /mnt/d/scripts/get_coverage_v5.3.py -s $sra_id \
		-t /mnt/d/trinityrnaseq-Trinity-v2.8.5/Trinity \
		-sn Nicotiana_bonariensis \
		-rf /mnt/d/ref_files/Nt_rbcL.fa,/mnt/d/ref_files/Nt_rbcS1.fa \
		-st rbcL,rbcS \
		-at 350
done
## More sra_ids and species can be added here like the one above to automate the process.
## More instructions are below.
## Make sure Trinity, bbmap.sh, fastq-dump, biophython, seqtk, numpy, pandas, matplotlib are properly installed for Ubuntu.
## By default, the data will be saved in /mnt/d/data/sn/sra_id/
## For each python script call, one scientific name (sn) must be provided.
## For each python script call, multiple reference files (rf) can be provided.
## For each reference_file (rf), a folder will be created to save its associated data.
## For each reference_file (rf), one seq_type (st) must be providied for naming purpose.
## For each seq_type (st), a csv file will be generated in /mnt/d/data/
## If a csv file with the same name already exists, the new data will be appended to the existing csv file.
## When ready, stay in the /home/user/ directory of ubuntu and execute this shell script file like "/directory_to_the_file/run_get_cov_4.sh". This is important because fastq-dump downloads sra file and saves the cache files in /home/user/ncbi/public/sra folder. The script will remove those cache files to prevent the hard drive from running out of space.
 

