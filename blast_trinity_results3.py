# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 21:54:37 2019

@author: Myat
edited by Heidi
"""

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio import pairwise2
#more info at http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc86
#http://biopython.org/DIST/docs/api/Bio.pairwise2-module.html
import os
import matplotlib.pyplot as plt
import pandas as pd
import csv
from Bio.Seq import Seq

#makes a list of secquence records from a fasta file
#prints the number of assemblies in the list
def make_seq_record_list(fasta_seq_file, str4print):
    seq_rec_list = list(SeqIO.parse(fasta_seq_file, "fasta", IUPAC.unambiguous_dna))
    print(r"Successfully imported %i %s sequences from %s." %(len(seq_rec_list), str4print, fasta_seq_file))
    return seq_rec_list

def get_best_align(tri_seq_rec_list, act_seq_rec_list, match_score, aln_score_th):
    best_list = [] #store best alignment for each trinity seq
    for tri_seq in tri_seq_rec_list:
        #resetting act_count and perfect_match
        act_count = 0
        for act_seq in act_seq_rec_list:
            #calculate the score for full perfect alignment
            full_score = len(act_seq)*match_score
            alignment = align_2_seq_rec(act_seq, tri_seq, match_score)
            if act_count == 0:
                #checking alignment with reverse complement
                alignmentR = align_2_seq_rec(act_seq, tri_seq.reverse_complement(), match_score)
                #if the reverse complement has a better alignment, it will be used
                #for subsequent alignments with the remaining actual sequences
                if alignmentR[2] > alignment[2]:
                    cur_tri_id = tri_seq.id + "___RC"
                    tri_seq.seq = tri_seq.reverse_complement().seq
                    tri_seq.id = cur_tri_id
                    alignment = alignmentR
                #print("analzying %s" %tri_seq.id)
                cur_best_aln = alignment
                best_act_seq = act_seq
                #if alignment with the first act seq is not good enough, skip the
                #alignment with remaining act seq. The reason is all rbcS genes are quite
                #similar to each other.
                if alignment[2] < aln_score_th-200:
                    break
                #will go to elif below for alignments with the remaining act_seq
                act_count = 1
            elif cur_best_aln[2] < alignment[2]:
                #update the cur_best_aln with a better alignment
                cur_best_aln = alignment
                best_act_seq = act_seq
            #print("score for %s is %f and current best is %s" %\
            #      (act_seq.id, alignment[2], best_act_seq.id))
            #if full alignment score is achieved, perfect_match will be updated and
            #the remaining actual sequences will be skipped
            if cur_best_aln[2]==full_score:
                #print("full alignment score achieved, will stop checking with the "+\
                #      "remaining actual sequences.")
                break
        align_len = len(cur_best_aln[0][cur_best_aln[3]:cur_best_aln[4]].replace("-",""))        
    	#if the best score is greater than the threshold and >95% of the full act seq is aligned,
        #the portion of tri_seq aligned to act_seq will be extracted and saved in aln_tri_list
        if cur_best_aln[2] > aln_score_th-200 and align_len > 0.90*len(best_act_seq):
            full_score = align_len*match_score
            full_tri_seq = str(tri_seq.seq)
            #extracting the trinity sequence aligned and removing the hyphens
            tri_seq.seq = Seq(cur_best_aln[1][cur_best_aln[3]:cur_best_aln[4]].replace("-",""),
                              IUPAC.unambiguous_dna)
            #save the act seq with best score in the list, also if it is perfectly aligned or not
            best_list.append([tri_seq, cur_best_aln[2]==full_score, cur_best_aln[2], best_act_seq.id, full_tri_seq])
    return best_list
	
def align_2_seq_rec(seq_rec1, seq_rec2,match_score):
    return pairwise2.align.localms(seq_rec1.seq.upper(),
                                   seq_rec2.seq.upper(),
                                   match_score, -3, -5, -2,
                                   penalize_extend_when_opening = True,
                                   one_alignment_only = True)[0]

def write_csv(best_list, outdir, act_seq_known, seq_type):
	#has_header = os.path.isfile(outdir + seq_type + 'best_list.csv')
	if act_seq_known:
		headers = ['scientific name','SRA ID','trim','read length',
			'bbmap reference file','num reads in bbsp stdout',
			'assembled trinity seq name','trinity seq orf',
			'whole trinity seq','best aligned act seq',
			'num aligned nucleotides','perfect alignment',
			'read coverage']
		file_name_ext = '_best_list.csv'
	else:
		headers = ['scientific name','SRA ID','trim','read length',
			'bbmap reference file','num reads in bbsp stdout',
			'assembled trinity seq name','trinity seq orf',
			'whole trinity seq','num aligned nucleotides',
			'read coverage']
		file_name_ext = '_predict_list.csv'
	with open(outdir + seq_type + file_name_ext, mode='a') as best_list_file:
		best_list_writer = csv.writer(best_list_file)
		if os.stat(outdir + seq_type + file_name_ext).st_size == 0:
			best_list_writer.writerow(headers)
		for my_list in best_list:
			if not act_seq_known:
				my_list.pop(-2)
				my_list.pop(-3)
			best_list_writer.writerow(my_list)	

def getBaseCov(fq_files, ref_seq_rec, basecov_out_file, rL, outdir):
    "This function will generate starting base coverage for ref_seq_rec " +\
    "using BBMap perfectmode and the two fq input files. It will save the " +\
    "base coverage info in basecov_out_file as well as the plot in png file."+\
    "Make sure ref_seq_rec is not shorter than rL."
    bbmap_sh = "bbmap.sh"
    #This is where files are stored temporarily for bbmap
    #out_dir = "/mnt/d/blast/blast_assemble_files/"
    SeqIO.write(ref_seq_rec, "temp_ref_bbmap_pm.fa", "fasta")
    if len(fq_files)==2:
        bbmap_cmd = bbmap_sh + \
                    ' in=' + fq_files[0] + ' in2=' + fq_files[1] + \
                    ' ref=' + 'temp_ref_bbmap_pm.fa' + \
                    ' nodisk threads=1 perfectmode -Xmx800m' + \
                    ' outm=' + 'temp_bbmap_pm.fa' + \
                    ' basecov=' + basecov_out_file + ' startcov=t' + \
                    ' overwrite=t'
    else:
        print('Make sure fq_files is a list of two fq files.')
        os.remove("temp_ref_bbmap_pm.fa")
        return None
    runSubprocess(bbmap_cmd, stdout_file = outdir + 'bbmap_pm_stdout.txt')
    os.remove("temp_ref_bbmap_pm.fa")
    os.remove('temp_bbmap_pm.fa')
    basecov_data = pd.read_csv(basecov_out_file, sep='\t')
    basecov_data = basecov_data.iloc[:-1*rL+1,-2:].values
    #writing to the new basecov txt file in csv format
    with open(basecov_out_file, 'w', newline='') as csvfile:
        basecov_writer = csv.writer(csvfile)
        basecov_writer.writerow(["Position","Starting Base Coverage"])
        for basecov in basecov_data:
            basecov_writer.writerow(basecov)
    #plotting the basecov_data
    base = basecov_data[:, 0:1]
    freq = basecov_data[:, 1]
    plt.figure()
    plt.ylabel("Frequency")
    plt.xlabel("Base position")
    plt.scatter(base, freq, color = 'red', marker='.')
    plt.savefig(basecov_out_file[:-4]+'.png', dpi=150, transparent=False)
    plt.pause(0.05)
    plt.close()
    return basecov_data

def runSubprocess(cmd_str, stdout_file = "stdout.txt"):
    "This function will call subprocess.run with the cmd_str input and save the stderr " + \
    "in the stdout_file."
    import subprocess
    run_output = subprocess.run(cmd_str, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    outfile_handle = open(stdout_file,'wb')
    outfile_handle.write(run_output.stderr)
    outfile_handle.write(run_output.stdout)
    outfile_handle.close()
    return None




