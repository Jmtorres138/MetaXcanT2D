#!/usr/bin/python -O 
# Jason Matthew Torres
'''
Extract GTEx eQTL values for MetaXcan model SNPs   
Usage: python JTgetGTEx_eQTL.runjobs.py  
'''
# libraries 
import sys, os, gzip 
import subprocess as sp 

# globals 
tarbell_dir = "/group/im-lab/nas40t2/jason/projects/MetaXcan/results/predictors/gtexV6p/"

key_file = tarbell_dir + "model-key.txt"


# function 

def make_mod_key_dic():
	print "Making GTEx model key.."
	dic = {}
	fin = open(key_file, 'r')
	for line in fin:
		l = line.strip().split()
		mtxn, gtex = l[0], l[1]
		dic[mtxn] = gtex 
	return dic

def main():

	fin = open(tarbell_dir+"genemods.txt",'r')
	fin.readline()
	for line in fin:
		l = line.strip().split()
		chromosome, locus, locus_start, locus_end, genename, model, ensid  = l[0],l[1],l[2],l[3],l[4],l[5],l[6]
		arg_list = ["python", tarbell_dir+"JTlocusGTEx_eQTL.py", model, genename, ensid, locus_start, locus_end, chromosome]
		arg_line = " ".join(arg_list)
		fname = tarbell_dir+model+"_"+genename+".temp.sh"
		fout = open(fname,'w')
		script = '''
			#!/bin/bash                                                                                                                              

			#PBS -N %s                                                                                                                            
			#PBS -S /bin/bash                                                                                                                        
			#PBS -l walltime=24:00:00                                                                                                                
			#PBS -l nodes=1:ppn=10                                                                                                                   
			#PBS -l mem=1gb                                                                                                                         
			#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log                                                                                            
			#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err                                                                                            

			cd $PBS_O_WORKDIR

			module load python/2.7.9

			%s 
			''' % (model+"_"+genename, arg_line)
		fout.write(script+"\n")
		fout.close()
		command = ["qsub", fname]
		if model != "TW_Whole_Blood_DGN": 
			sp.check_call(command) 

if (__name__=="__main__"): main()

