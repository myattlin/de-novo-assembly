#!/bin/sh

#############################################################
##  Run get_coveragy_vxxx.py gets covereage on a trinity output ##
#############################################################

#sra_ids="SRR6918807 SRR6516146" 
#for sra_id in $sra_ids; do
#	python /mnt/d/in_run_get_cov/get_coverage_v5.3.py -s $sra_id -t Trinity \
#		-sn Nicotiana_bonariensis \
#		-rf /mnt/d/blast/ref_files/Nt_rbcS_T.fa,/mnt/d/blast/ref_files/Nt_rca1a.fa \
#		-st rbcS_T,rca \
#		-at 350
#done
sra_ids="ERR274390" 
for sra_id in $sra_ids; do
	python /mnt/d/in_run_get_cov/get_coverage_v5.3.py -s $sra_id -t Trinity -c false \
		-sn Nicotiana_sylvestris \
		-rf /mnt/d/blast/ref_files/Nt_rbcS1.fa \
		-st rbcS \
		-at 350
	python /mnt/d/in_run_get_cov/get_coverage_v5.3.py -s $sra_id -t Trinity -c false \
		-sn Nicotiana_sylvestris \
		-rf /mnt/d/blast/ref_files/Nt_rbcL.fa \
		-st rbcL \
		-at 350
	python /mnt/d/in_run_get_cov/get_coverage_v5.3.py -s $sra_id -t Trinity -c false \
		-sn Nicotiana_sylvestris \
		-rf /mnt/d/blast/ref_files/Nt_rca1a.fa \
		-st rca \
		-at 350
	python /mnt/d/in_run_get_cov/get_coverage_v5.3.py -s $sra_id -t Trinity \
		-sn Nicotiana_sylvestris \
		-rf /mnt/d/blast/ref_files/Nt_rbcS_T.fa \
		-st rbcS_T \
		-at 350
done
