#!/usr/bin/python -O 
# Jason Matthew Torres
'''
Extract GTEx eQTL values for MetaXcan model SNPs   
Usage: python JTgetGTEx_eQTL.py  model genename ensid locus_start locus_end chromosome 
'''
# libraries 
import sys, os, gzip 
import subprocess as sp 

# globals 
tarbell_dir = "/group/im-lab/nas40t2/jason/projects/MetaXcan/results/predictors/gtexV6p/"
out_dir = tarbell_dir + "gtex_files/"
gtex_eqtl_dir = "/group/im-lab/nas40t2/Data/dbGaP/GTEx/V6p-gtexportal/GTEx_Analysis_v6p_all-associations/"

key_file = tarbell_dir + "model-key.txt"
res_file = tarbell_dir + "mtxnsig-gene-pred.txt"
gtex_snp_file = "/group/im-lab/nas40t2/Data/dbGaP/GTEx/V6-gtexportal/GTEx_OMNI_genot_1KG_imputed_var_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT_release_v6.txt.gz"

model = sys.argv[1]
genename = sys.argv[2] 
ensid = sys.argv[3]
locus_start = sys.argv[4]
locus_end = sys.argv[5]
chromosome = sys.argv[6]

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

def make_snp_key_dic():
	print "Building SNP key..."
	fin = gzip.open(gtex_snp_file,'rb')
	header = fin.readline()
	dic = {}
	for line in fin:
		l = line.strip().split()
		chromo, pos, gtex_id, rsid = l[0], l[1], l[2], l[6]
		if int(chromo) == int(chromosome) and int(pos) >=  int(locus_start) and int(pos) <= int(locus_end):
			dic[gtex_id] = [rsid, pos, chromo]
	return dic

def main():
	print model
	mod_key = make_mod_key_dic()
	tissue = mod_key[model]
	snp_key = make_snp_key_dic()
	fname = gtex_eqtl_dir + tissue + "_Analysis.v6p.all_snpgene_pairs.txt.gz"
	# Write header file 
	headlist = ["gene.id","gtex.id","tss.distance","pval.nominal","slope","slope.se"]
	fout1 = open(out_dir+"header.txt",'w')
	fout1.write("\t".join(headlist)+"\n")
	print fname 
	tempname = out_dir+tissue+"_"+genename+".temp.txt"
	command1 = "LC_ALL=C zfgrep " + ensid + " " +  fname + " > "  + tempname
	sp.check_call(command1, shell=True)
	fin = open(tempname,'r')
	fout = gzip.open(out_dir+model+"_"+genename+".gtex.txt.gz",'wb')
	header = ["rsid","pos","chrom","model","ensid","gtex.id","tss.distance","pval.nominal","slope","slope.se"]
	#fout.write("\t".join(header)+"\n")	
	for line in fin:
		l = line.strip().split()
		eid, snp, tssd, p, beta, se = l[0],l[1],l[2],l[3],l[4],l[5]
		if eid == ensid:
			try:
				rsid = snp_key[snp][0]
				pos = snp_key[snp][1]
				chromo = snp_key[snp][2]
				write_list = [rsid, pos, chromo, tissue, ensid, snp, tssd, p, beta, se]
				fout.write("\t".join(write_list)+"\n")
			except:
				print "SNP %s is not in dictionary" % snp
	fout.close() 
	fin.close()
	os.remove(tempname)
	os.remove(out_dir+"header.txt")     	







if (__name__=="__main__"): main()

