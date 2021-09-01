# -*- coding: utf-8 -*-
"""
created on Wed Jun 12 10:07 2019

@author: Heidi
code fragments taken from Myat
"""

import subprocess
import os


def fastq_dump(ncbi_id, fd_outdir):
	"This function uses fastq dump to download the fastq file for the specified sra_id"
	print ("running fastq-dump")
	fd_cmd = "fastq-dump --defline-seq '@$ac.$si.$sg/$ri' --defline-qual" + \
		 " '+' --split-files --outdir " + fd_outdir + " " + ncbi_id
	runSubprocess(fd_cmd, stdout_file = fd_outdir + "fq_stdout.txt")


def bbmap(sra_id, bbmap_ref, bbmap_ref_file, bbmap_outdir, indir):
	"This function uses bbmap to find the reads in the SRA fastq file that" + \
	"resemble small subunits"
	#fastq files containing the RNA-seq data from the specified sra_id
	fq_files = [indir + sra_id + "_1.fastq", indir + sra_id + "_2.fastq"]
	#modified for RNA-Seq run at Cornell
	#fq_files = [indir + sra_id + "_R1.fastq.gz", indir + sra_id + "_R2.fastq.gz"]
	bbmap_file_prefix = bbmap_outdir + "bbmap_" + bbmap_ref[:-3] 
	bbmap_cmd = 'bbmap.sh' + \
                   ' in=' + fq_files[0] + ' in2=' + fq_files[1] + \
                   ' vslow ref=' + bbmap_ref_file + \
                   ' nodisk threads=2 maxindel=100 strictmaxindel=t local=t samplerate=1 ' + \
                   ' outm=' + bbmap_file_prefix + ".fq ignorebadquality"
	print("running bbmap: "+bbmap_cmd+"\n")	
	runSubprocess(bbmap_cmd, stdout_file = bbmap_file_prefix + "_stdout.txt")

	#split the resulting bbmap file into two files for the foward and reverse directions
	bbsp_out1 = bbmap_file_prefix + "_1.fq"
	bbsp_out2 = bbmap_file_prefix + "_2.fq"
	bbsp_cmd = 'bbsplitpairs.sh' + \
		   ' in=' + bbmap_file_prefix + ".fq" +\
		   ' out=' + bbsp_out1 +\
		   ' out2=' + bbsp_out2
	print("running bbsplitpairs: " + bbsp_cmd+"\n")
	runSubprocess(bbsp_cmd, stdout_file = bbmap_outdir + "bbsp_" + bbmap_ref[:-3] + "_stdout.txt")

	#deleat the original bbmap file because same informatio as in the split files
	os.remove(bbmap_file_prefix + ".fq")

	return(bbsp_out1, bbsp_out2)

def bb5000(outdir, bbsp_out1, bbsp_out2):
	"This function takes the bbsplitpairs outputs and extracts the first 5000 lines which" + \
	"which are stored in anouther file"
	seqtk_out1 = bbsp_out1[:-3] + "_5k" + bbsp_out1[-3:]
	seqtk_out2 = bbsp_out2[:-3] + "_5k" + bbsp_out2[-3:]
	seqtk_cmd1 = "seqtk sample -s 100 " + bbsp_out1 + " 5000 > " + seqtk_out1
	seqtk_cmd2 = "seqtk sample -s 100 " + bbsp_out2 + " 5000 > " + seqtk_out2
	seqtk_stdout1 = outdir + "seqtk_stdout1.txt"
	seqtk_stdout2 = outdir + "seqtk_stdout2.txt"
	runSubprocess(seqtk_cmd1, stdout_file = seqtk_stdout1)
	runSubprocess(seqtk_cmd2, stdout_file = seqtk_stdout2)
	return(seqtk_out1, seqtk_out2)
	
def bbtrim(outdir, bbsp_out1, bbsp_out2, left, right):
	"This function trims the bbsp fq files using seqtk."
	bbtrim_out1 = bbsp_out1[:-3] + "t" + bbsp_out1[-3:]
	bbtrim_out2 = bbsp_out2[:-3] + "t" + bbsp_out2[-3:]
	bbtrim_cmd1 = "seqtk trimfq -b " + str(left) + " -e " + str(right) + " " +\
					bbsp_out1 + " > " + bbtrim_out1
	bbtrim_cmd2 = "seqtk trimfq -b " + str(left) + " -e " + str(right) + " " +\
					bbsp_out2 + " > " + bbtrim_out2	
	runSubprocess(bbtrim_cmd1)
	runSubprocess(bbtrim_cmd2)
	return(bbtrim_out1, bbtrim_out2)

def trinity(tri_exe, outdir, bbmap_ref, tri_opts, left_file, right_file, opts_type):
	"This function uses the trinity program to de novo assemble the small subunits"
	tri_file_prefix = outdir + "Trinity_" + opts_type + bbmap_ref[:-3]
	left_file = " --left " + left_file
	right_file = " --right " + right_file
	tri_outdir = " --output " + tri_file_prefix + "/"
	trinity_cmd = tri_exe + tri_opts + left_file + right_file + tri_outdir
	print("running trinity " + opts_type + ": " + trinity_cmd+"\n")
	runSubprocess(trinity_cmd, stdout_file = tri_file_prefix + "_stdout.txt")
	return(tri_file_prefix + ".Trinity.fasta")


def runSubprocess(cmd_str, stdout_file = "stdout.txt"):
	"This function will call subprocess.run with the cmd_str input and save the stderr " + \
	"in the stdout_file."
	run_output = subprocess.run(cmd_str, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	outfile_handle = open(stdout_file,'wb')
	outfile_handle.write(run_output.stderr)
	outfile_handle.write(run_output.stdout)
	outfile_handle.close()
