# -*- coding: utf-8 -*-
"""
created on Wed Jun 12 10:07 2019

@author: Heidi
code fragments taken from Myat

input sra file --> extracts scquences that resemble the small subunit for rubisco using bbmap -->
trinity de novo assembles the sequences --> use pairwise2 (similar to blast) to compair actual small subunit sequences with the assembled ones to determine which are properly assembled -->
bbmap allignes maps the coverage of the trinity outputs --> out coverage maps for the trinity 
transcripts and a master file containing data
"""

import argparse
import fd_bbmap_trinity_funct as fbt
import blast_trinity_results3 as btr3
from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
import os
#import gzip
#import shutil

if __name__ == "__main__":

	#gets user arguments and sets to defults if not inputed

	parser = argparse.ArgumentParser()

	#adding a mandatory input argument for sra's ID "SRR6364674"
	parser.add_argument('-s', '--sra_id', help='SRA to be processed')

	#adding an imput argument for the organisms scientific name
	parser.add_argument('-sn', '--scientific_name', \
			    help="scientific name of the organism with the specified SAR id")
	
	#adding an argument for the reference file for bbmap
	parser.add_argument('-rf','--reference_file', help='reference file for bbmap')

	#adding an argument for output directory
	parser.add_argument('-od','--outdir', help='output directory for all program outputs')

	#adding an optional argument for the Trinity software location
	parser.add_argument('-t','--trinity_exe', help='location of the Trinity program file')

	#adding an optional argument for the actual sequence files "/mnt/d/test/Nt_rbcS_M.fa"
	parser.add_argument('-af','--act_seq_file', help='actual sequence file')

	#adding an optional argment for match score in blast
	parser.add_argument('-mc','--match_score', help='score given for a match in blast')

	#adding anoptional argument for allignment score threshold must be at least 200
	parser.add_argument('-at','--aln_score_th', help='alignment score threshold must \
			    be at least 200')
	
	#adding argument for what type of sequences are being searched for eg. small subumit
	parser.add_argument('-st','--seq_type', help="which sequencces are being searched for \
			    eg. small subunit (rbcS), large subunit (rbcL)")
	
	#adding an optional argument for trimming the fastq files after bbmap
	parser.add_argument('-tr', '--trim', help='to trim the bbsp_out files')
    
    #adding an optional argument for trimming the fastq files after bbmap
	parser.add_argument('-cln', '--clean', help='to clean up large fastq files')
	
	args = parser.parse_args()

	if args.sra_id:
		sra_id = args.sra_id
		#print("\n")		
	else: raise argparse.ArgumentTypeError('missing sra input')

	if args.scientific_name:
		sci_name = args.scientific_name
		#print("SRA is %s \n" % sra_id)
	else: sci_name = "NA"
	
	print("SRA to be analyzed is " + sra_id + " from " + sci_name+". \n")

	if args.reference_file:
		bbmap_ref_file_list = args.reference_file
		bbmap_ref_file_list=bbmap_ref_file_list.split(',')
	else: raise argparse.ArgumentTypeError('missing bbmap reference file input')
	#else: bbmap_ref_file_list = ['/mnt/d/test/Sl_rbcS.fa']
	
	if args.outdir:
		outdir_original = args.outdir
	else: outdir_original = "/mnt/d/data/"

	if args.trinity_exe:
		tri_exe = args.trinity_exe
	else: tri_exe = "/mnt/d/blast/trinityrnaseq-Trinity-v2.8.5/Trinity"
    
	#this is problematic for Trinity installed within Linux
	#if not os.path.isfile(tri_exe):
	#	raise argparse.ArgumentTypeError('tri_exe file location is incorrect.')

	if args.act_seq_file:
		act_seq_file_list = args.act_seq_file
		act_seq_file_list = act_seq_file_list.split(',')
		act_seq_known = True
	else: act_seq_known = False

	if args.match_score:
		match_score = int(args.match_score)
	else: match_score = 2

	if args.aln_score_th:
		aln_score_th = int(args.aln_score_th)
	else: aln_score_th = 500 #alignment score threshold must be at least 200

	if args.seq_type:
		seq_type = args.seq_type
		seq_type=seq_type.split(',')
	else: seq_type = ['rbcS']

	#makes 2 subdirectories in the directory requested by the user using the scientific 
	#name of the organism and the SRA id
	indir = outdir_original + sci_name + "/" + sra_id + "/"
	
	#print(indir)
	#os.mkdir(indir)
	#print("Directory", indir, "made")	
	try:
		os.mkdir(indir)
		print("Directory", indir, "made \n")
	except:
		pass

	
    #extracts the data from SRA for the specific entered sra id
	if not os.path.isfile(indir+sra_id+"_1.fastq"):
		fbt.fastq_dump(sra_id, indir)

	#removes the sra cache file that was saved in Linux folder
	#this will be problematic if the sra folder does not exist
	cache_files = os.listdir("./ncbi/public/sra/")
	for cache_file in cache_files:
		os.remove("./ncbi/public/sra/"+cache_file)
    #if os.path.isfile("./ncbi/public/sra/" + sra_id + ".sra"):
	#	os.remove("./ncbi/public/sra/" + sra_id + ".sra")
	
	#different trinity options to try
	defult_tri_opts = " --seqType fq --max_memory 2G --CPU 2 --trimmomatic "+\
			    "--full_cleanup --SS_lib_type RF  --normalize_max_read_cov 30"
	tri_opts_2 = defult_tri_opts + " --KMER_SIZE 32"
	tri_opts_3 = defult_tri_opts + " --min_kmer_cov 4 --min_glue 4 --min_iso_ratio 0.2 "+\
					 "--glue_factor 0.2 --jaccard_clip"
	tri_opts_4 = defult_tri_opts + " --KMER_SIZE 32 --min_kmer_cov 4 --min_glue 4 "+\
					 "--min_iso_ratio 0.2 --glue_factor 0.2 --jaccard_clip"
	#print(tri_opts_3)
	
	#the actual sequence(s) entered by the user will be in a list that needs to change
	#every iteration of the for loop so this is the current index in the list
	loop_count = 0
	#loops through the reference files entered by the user to extract data from the fastq
	#dump output
	for bbmap_ref_file in bbmap_ref_file_list:
		#extracts the name of the bbmap referance file to be used in naming other files
		index_of_slash = bbmap_ref_file.rindex('/')+1
		bbmap_ref = bbmap_ref_file[index_of_slash:]
		#makes a subdirrectory for the specific bbmap referance
		outdir = indir + bbmap_ref[:-3] + "/"
		try:
			os.mkdir(outdir)
			print("Directory", outdir, "made")
		except:
			pass

		#runs bbmap to extract the sequences similar to the reference sequence
		#outputs 2 files that contan the extracted sequences
		bbsp_out1, bbsp_out2 = fbt.bbmap(sra_id, bbmap_ref, bbmap_ref_file, outdir, indir)
		
		#if bbmap was already ran, use these below
		bbsp_out1 = outdir + "bbmap_" + bbmap_ref[:-3] + "_1.fq"
		bbsp_out2 = outdir + "bbmap_" + bbmap_ref[:-3] + "_2.fq"
		
		#trimming the fq files and reassigning bbsp_out1 and bbsp_out2
		if args.trim:
			trim_lr = args.trim
			trim_lr = trim_lr.split(',')
			bbsp_out1, bbsp_out2 = fbt.bbtrim(outdir, bbsp_out1, bbsp_out2,
									 int(trim_lr[0]), int(trim_lr[1]))
			print("bbsp fq files trimmed at left by " + trim_lr[0] + " and " +\
				 "right by " + trim_lr[1] + "\n")
		else: trim_lr = [0,0]
		
		#finds the length of the reads from the first bbsplitpairs output file
		readLength = len(next(SeqIO.parse(bbsp_out1,'fastq')))
		print("Read length = " + str(readLength) + "\n")

		#bbmap will be skipped and unzipped gz files used as trinity input files
		#bbsp_out1 = indir + sra_id + "_R1.fq"
		#bbsp_out2 = indir + sra_id + "_R2.fq"
	
		#with gzip.open(indir + sra_id +"_R1.fastq.gz", 'rb') as f_in:
		#	with open(bbsp_out1, 'wb') as f_out:
		#		shutil.copyfileobj(f_in, f_out)
		#with gzip.open(indir + sra_id +"_R2.fastq.gz", 'rb') as f_in:
		#	with open(bbsp_out2, 'wb') as f_out:
		#		shutil.copyfileobj(f_in, f_out)
						
		'''
		#code for when dont want to rerun bbmap
		bbmap_file_prefix = outdir + "bbmap_" + bbmap_ref[:-3]
		bbsp_out1 = bbmap_file_prefix + "_1.fq"
		bbsp_out2 = bbmap_file_prefix + "_2.fq"
		seqtk_out1 = bbsp_out1[:-3] + "_5k" + bbsp_out1[-3:]
		seqtk_out2 = bbsp_out2[:-3] + "_5k" + bbsp_out2[-3:]
		'''
		'''
		#code for when skipping trinity
		tri_fasta_list =[]
		for i in range(6):
			tri_fasta_list.append(outdir+"Trinity_"+str(i+1)+"_"+\
						 bbmap_ref[:-3]+".Trinity.fasta")
		'''	
		
		#makes a list for the fasta outputs of all the trinity runs
		tri_fasta_list = []
		#running the trinity program with different oprions then storing the 
		#output in a list
		#skipping trinity run # 1
		#tri_fasta_temp = fbt.trinity(tri_exe, outdir, bbmap_ref, defult_tri_opts ,
        #                       bbsp_out1, bbsp_out2, "1_")		
		#tri_fasta_list.append(tri_fasta_temp)
		tri_fasta_temp = fbt.trinity(tri_exe, outdir, bbmap_ref, tri_opts_2 , 
						bbsp_out1, bbsp_out2, "2_")		
		tri_fasta_list.append(tri_fasta_temp)
		tri_fasta_temp = fbt.trinity(tri_exe, outdir, bbmap_ref, tri_opts_3 , 
					     bbsp_out1, bbsp_out2, "3_")		
		tri_fasta_list.append(tri_fasta_temp)
		tri_fasta_temp = fbt.trinity(tri_exe, outdir, bbmap_ref, tri_opts_4 , 
					     bbsp_out1, bbsp_out2, "4_")		
		tri_fasta_list.append(tri_fasta_temp)
		
		#finds the number of reads in the bbsplitpairs standard out file
		bbsp_stdout = outdir + "bbsp_" + bbmap_ref[:-3] + "_stdout.txt"
		readFile = open(bbsp_stdout)
		lines = readFile.readlines()
		lines_split = filter(None,[line.split() for line in lines])
		num_reads = int([line[1] for line in lines_split if line[0]=='Result:'][0])
		#num_reads = 10001
		print("num_reads with bbmap = " + str(num_reads) + "\n")
		
		#if there are more than 10000 reads only the first 5000 are extracted then
		#used in trinity
		if num_reads > (2 * 5000):
			seqtk_out1, seqtk_out2 = fbt.bb5000(outdir, bbsp_out1, bbsp_out2)
			tri_fasta_temp = fbt.trinity(tri_exe, outdir, bbmap_ref, tri_opts_2 ,seqtk_out1, seqtk_out2, "5_")
			if os.path.isfile(tri_fasta_temp):
				tri_fasta_list.append(tri_fasta_temp)
			tri_fasta_temp = fbt.trinity(tri_exe, outdir, bbmap_ref, tri_opts_4 ,seqtk_out1, seqtk_out2, "6_")
			if os.path.isfile(tri_fasta_temp):
				tri_fasta_list.append(tri_fasta_temp)
		
		#import the actual or known sequences into a list of Sequence Records
		if act_seq_known:
			act_ref_seq_rec_list = btr3.make_seq_record_list(
					act_seq_file_list[loop_count], 
					"actual")
			print("aligning with act seqs ...")
		else:
			ref_seq_rec_list = btr3.make_seq_record_list(
					bbmap_ref_file, 
					"bbmap ref")
			print("missing actual sequences, aligning with ref seq ...")
		#the number trinity sequence that the loop is currently on
		trinity_count = 1
		#loops through all of the trinity files ranto get their base coverage
		for trinity_fasta in tri_fasta_list:
			#import the trinity_fasta sequences into a list of Sequence Records
			if os.path.isfile(trinity_fasta):
				tri_seq_rec_list = btr3.make_seq_record_list(trinity_fasta, "Trinity")
			else:
				trinity_count += 1
				continue

			#obtain the best aligned real sequence to each of the trinity sequences
			if act_seq_known:				
				best_list = btr3.get_best_align(tri_seq_rec_list, 
									act_ref_seq_rec_list, 
									match_score, 
									aln_score_th)
			else:				
				best_list = btr3.get_best_align(tri_seq_rec_list, 
									ref_seq_rec_list, 
									match_score, 
									aln_score_th)				
			print("%i trinity transcripts selected from trinity run %i" %\
		       (len(best_list), trinity_count+1))#add 1 if run 1 is skpped
			
			#finds the read coverage for each trinity sequence that met the threshold
			#all important aspects of data are recorded in the csv file
			print("generating read coverages ...\n")
			for my_list in best_list:
				my_list[0].id = my_list[0].id + "_" +str(trinity_count+1)#add 1 if run 1 is skipped
				basecov_data = btr3.getBaseCov([bbsp_out1, bbsp_out2], 
								   my_list[0], 
								   outdir+my_list[0].id + "_baseCov.txt", 
								   readLength, 
								   outdir)
				object_my_list = my_list[0].id
				my_list.extend([str(my_list[0].seq),my_list[0].id, num_reads, 
					bbmap_ref, readLength, str(trim_lr[0])+","+str(trim_lr[1]), 
					sra_id, sci_name])
				my_list.reverse()
				my_list.append(outdir+ object_my_list + "_baseCov.txt")
				my_list.pop(-2)
				
			trinity_count += 1

			# makes csv files where aln_tri_list and best_list data is stored	
			btr3.write_csv(
					best_list,
					outdir_original,
					act_seq_known,
					seq_type[loop_count])
		
		#generate baseCov for the actual sequences
		if act_seq_known:
			print("generating read coverages for actual sequences ... \n")
			for act_seq in act_ref_seq_rec_list:
				btr3.getBaseCov([bbsp_out1, bbsp_out2], act_seq, outdir+act_seq.id + "_baseCov.txt", readLength, outdir)
					
		#index of the sequence type in the users inputed sequence list is updated
		loop_count += 1
	
	#original sra files are removed if necessary
	clean_up = True
	if args.clean:
		if args.clean.lower() in ["false", "f", "no", "n"]:
			clean_up = False
	if clean_up:
		os.remove(indir + sra_id + "_1.fastq")
		os.remove(indir + sra_id + "_2.fastq")
