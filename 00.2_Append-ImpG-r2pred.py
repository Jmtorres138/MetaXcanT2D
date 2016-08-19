#!/bin/usr/python -O 
#Jasom Matthew Torres
'''
Input model-snps-gained.txt and append r2pred values for each SNP 
'''

# libraries
import sys, os, gzip  
from time import sleep

mysuffix = sys.argv[1]

root_dir = "/group/im-lab/nas40t2/jason/projects/MetaXcan/meta_files/ImpGformat/"
work_dir = root_dir + "merged_output_files/impG_DIAGRAM/"
##fname = work_dir + "model-snps-gained.txt"
fname = work_dir + "model-snps-gained_gtexV6p.txt"


 #chr22/betas_DIAGRAM/"

def build_ref_dic():
	print "Building reference dictionary..."
	dic = {}
	fin = open(fname,'r')
	for line in fin:
		snp = line.strip()
		dic[snp] = snp
	return dic

def write_append():
	##fout = gzip.open(work_dir+"model-snps-gained"+mysuffix+".txt.gz",'wb')
	fout = gzip.open(work_dir+"model-snps-gained_gtexV6p"+mysuffix+".txt.gz",'wb')
	head_list = ["SNP", "POS", "Ref", "Alt", "Z", "r2pred"]
	header = "\t".join(head_list)
	fout.write(header+"\n")
	ref_dic = build_ref_dic()
	for i in range(1,23):
		print "\nProcessing chromosome: %d" % i
		imp_dir = root_dir+"chr"+str(i)+"/imp"+mysuffix+"/"
		flist = os.listdir(imp_dir)
		flist = [f for f in flist if ".map.impz" in f]
		flist = sorted(flist)
		count = 0
		for f in flist:
			count += 1 
			#print "File: %d" % count 
			sys.stdout.write("\r")
			sys.stdout.write("File: %d" % count)
			sys.stdout.flush()
			sleep(0.25)
			fin = open(imp_dir+f,'r')
			head = fin.readline() 
			for line in fin:
				l = line.strip().split()
				snp, pos, ref, alt, z, r2pred = l[0],l[1],l[2],l[3],l[4],l[5]
				if snp in ref_dic:
					write_list = [snp,pos,ref,alt,z,r2pred]
					fout.write("\t".join(write_list)+"\n")
			fin.close()
	fout.close()
	print("\nProcess complete")


def main():
	write_append()


if (__name__=="__main__"): main()
