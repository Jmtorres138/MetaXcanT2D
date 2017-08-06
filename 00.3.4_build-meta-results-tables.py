#!/usr/bin/python -O
# Jason Matthew Torres
'''
Script to combine MetaXcan meta-analysis results into summary table
Usage: python JTbuildTables.py meta_name alpha
Example: python 00.3.4_build-meta-results-tables.py DIAGRAM-GERA_ImpG_0.80_gtexV6p 0.5
'''
# libraries
import sys,os,gzip
import subprocess as sp

meta_name = sys.argv[1]
alpha = sys.argv[2]

# globals
root_dir = "/group/im-lab/nas40t2/jason/projects/MetaXcan/results/meta-analyses/DIAGRAM-GERA_ImpG_0.80_gtexV6p/"
out_dir = root_dir+"tables/"
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
res_dir = root_dir+"/output/alpha_"+str(alpha)+"/"
if not os.path.exists(res_dir):
    os.makedirs(res_dir)

file_list = os.listdir(res_dir)
tiss_list = sorted([x.split(".meta.txt")[0] for x in file_list])

def build_gene_dic():
    print("Building gene dictionary...")
    gene_list = []
    gene_dic = {}
    for tiss in tiss_list:
        fin = open(res_dir+tiss+".meta.txt",'r')
        fin.readline()
        for line in fin:
            l = line.strip().split()
            gene,ss_zscore,ss_pval,moresig = l[0],l[1],l[2],l[7]
            gene_list.append(gene)
            try:
                gene_dic[gene].append([tiss,ss_zscore,ss_pval,moresig])
            except:
                gene_dic[gene] = [[tiss,ss_zscore,ss_pval,moresig]]
        fin.close()
    gene_list = sorted(list(set(gene_list)))
    return gene_list, gene_dic

def write_outfile1(gene_list,gene_dic):
    print("Writing Zscore result table...")
    fout = gzip.open(out_dir+meta_name+"."+str(alpha)+".zscore.table.csv.gz",'wb')
    head_list = ["Gene"] + tiss_list
    fout.write(",".join(head_list)+"\n")
    for gene in gene_list:
        write_list = []
        write_list.append(gene)
        for tiss in tiss_list:
            check = False
            for l in gene_dic[gene]:
                if l[0] == tiss:
                    zscore = l[1]
                    check = True
            if check != True:
                zscore = "NA"
            write_list.append(zscore)
        fout.write(",".join(write_list)+"\n")
    fout.close()

def write_outfile2(gene_list,gene_dic):
    print("Writing Pvalue result table...")
    fout = gzip.open(out_dir+meta_name+"."+str(alpha)+".pvalue.table.csv.gz",'wb')
    head_list = ["Gene"] + tiss_list
    fout.write(",".join(head_list)+"\n")
    for gene in gene_list:
        write_list = []
        write_list.append(gene)
        for tiss in tiss_list:
            check = False
            for l in gene_dic[gene]:
                if l[0] == tiss:
                    pvalue = l[2]
                    check = True
            if check != True:
                pvalue = "NA"
            write_list.append(pvalue)
        fout.write(",".join(write_list)+"\n")
    fout.close()

def write_outfile3(gene_list,gene_dic):
    print("Writing 'More Significant' Evaluation result table...")
    fout = gzip.open(out_dir+meta_name+"."+str(alpha)+".moresig.table.csv.gz",'wb')
    head_list = ["Gene"] + tiss_list
    fout.write(",".join(head_list)+"\n")
    for gene in gene_list:
        write_list = []
        write_list.append(gene)
        for tiss in tiss_list:
            check = False
            for l in gene_dic[gene]:
                if l[0] == tiss:
                    moresig = l[3]
                    check = True
            if check != True:
                moresig = "NA"
            write_list.append(moresig)
        fout.write(",".join(write_list)+"\n")
    fout.close()


def main():
    gl, gd = build_gene_dic()
    write_outfile1(gl,gd)
    write_outfile2(gl,gd)
    write_outfile3(gl,gd)

if (__name__=="__main__"): main()
