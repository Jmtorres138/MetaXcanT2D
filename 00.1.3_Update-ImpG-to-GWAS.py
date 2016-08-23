#!/bin/usr/python -O 
#Jasom Matthew Torres
'''
Takes the gwas input file (for MetaXcan analysis) and appends ImpG-Summary imputed model/preditor SNPs
that meet a user-specified threshold  
'''

# libraries
import sys, os, gzip  
#from time import sleep

mysuffix = sys.argv[1]
threshold = sys.argv[2] 

work_dir = "/group/im-lab/nas40t2/jason/projects/MetaXcan/meta_files/ImpGformat/merged_output_files/impG"+mysuffix+"/"


#SNP	POS	Ref	Alt	Z	r2pred
#rs10399749	55299	C	T	-1.142090	0.140102
#rs12184279	717485	C	A	-0.131335	0.763391

def build_ref_dic():
	print "Building reference dictionary..."
	dic = {}
	fname = work_dir + "model-snps-gained_gtexV6p" + mysuffix + ".txt.gz"
	fin = gzip.open(fname,'rb')
	fin.readline() # header line 
	for line in fin:
		l = line.strip().split()
		snp, pos, ref, alt, z, r2pred = l[0], l[1], l[2], l[3], l[4], l[5] 
		val_list = [snp, pos, ref, alt, z, r2pred]
		if snp not in dic:
			dic[snp] = val_list
		elif snp in dic:
			if float(r2pred) > float(dic[snp][5]):
				dic[snp][4] = z
				dic[snp][5] = str(r2pred)
		else:
			print ("Something is wrong with SNP: %s" % snp)
	fin.close()
	return(dic)

def create_updated_gwasfile():
	print "Updating GWAS file..."
	if mysuffix == "_DIAGRAM":
		f = "/group/im-lab/nas40t2/jason/projects/MetaXcan/meta_files/diagram/hm2ceu/diag3_Z/diag3.z.txt.gz"
		out_name = "diag3.z.impG_"+str(threshold)+".gtexV6p.txt.gz"
	elif mysuffix == "_GERA":
		f = "/group/im-lab/nas40t2/jason/projects/MetaXcan/meta_files/gera/gera_Z/gera.z.txt.gz" # may need to manually create this file (i.e. add Z score column)
		# Z = log(OR) / SE 
		out_name = "gera.z.impG"+str(threshold)+".gtexV6ptxt.gz"
	else:
		print "Suffix must be either _DIAGRAM or _GERA" 
	track_dic = {} 
	fin = gzip.open(f,'rb')
	fin.readline() # header_line
	#SNP	A1	A2	FRQ	INFO	OR	SE	P	Z
	output_dir = work_dir+"gwas_file/"
	try:
		os.mkdir(output_dir)
	except:
		pass
	fout = gzip.open(output_dir+out_name,'wb')
	head_list = ["SNP", "A1", "A2", "Z", "r2pred"]
	fout.write("\t".join(head_list)+"\n")
	for line in fin:
		l = line.strip().split()
		snp, a1, a2, z = l[0], l[1], l[2], l[8]
		write_list = [snp, a1, a2, z, "1"]
		fout.write("\t".join(write_list)+"\n")
		track_dic[snp] = snp 
	fin.close()
	print ("Appended Imputed model SNPs...")
	ref_dic = build_ref_dic()
	for key in ref_dic:
		if key not in track_dic:
			klist = ref_dic[key] 
			snp, pos, ref, alt, z, r2pred = klist[0],klist[1],klist[2],klist[3],klist[4],klist[5]
			if float(r2pred) >= float(threshold):
				write_list = [snp,ref,alt,z,r2pred]
				fout.write("\t".join(write_list)+"\n")
	fout.close()
	print "Process complete"


def main():
	create_updated_gwasfile()


if (__name__=="__main__"): main()
