#!/usr/bin/python -O 
# Jason Matthew Torres
'''
Script to combine MetaXcan results into summary table 
Usage: python JTbuildTables.py meta_name alpha
Example: python JTbuildTables.py DIAGRAM1 0.5
'''
# libraries 
import sys,os,gzip 
import subprocess as sp 

meta_name = sys.argv[1] 
alpha = sys.argv[2] 

# globals 
root_dir = "/group/im-lab/nas40t2/jason/projects/MetaXcan/results/"
out_dir = root_dir + "tables/GTExV6p/"
res_dir = root_dir+meta_name+"/alpha_"+str(alpha)+"/"
file_list = os.listdir(res_dir)
tiss_list = sorted([x.split(".zscore")[0] for x in file_list])

def build_gene_dic():
    print("Building gene dictionary...")
    gene_list = []    
    gene_dic = {}
    for tiss in tiss_list:
        fin = open(res_dir+tiss+".zscores.csv",'r')
        fin.readline()
        for line in fin:
            l = line.strip().split(',')
            ensid,gene,zscore,pval,cvr2,n,ntot = l[0],l[1],l[2],l[4],l[6],l[9],l[11]
            gene_list.append(gene)
            try: 
                gene_dic[gene].append([tiss,zscore,pval,cvr2,n,ntot])
            except:
                gene_dic[gene] = [[tiss,zscore,pval,cvr2,n,ntot]]
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
    print("Writing Predicted R2 result table...")
    fout = gzip.open(out_dir+meta_name+"."+str(alpha)+".r2.table.csv.gz",'wb')
    head_list = ["Gene"] + tiss_list
    fout.write(",".join(head_list)+"\n")
    for gene in gene_list:
        write_list = []
        write_list.append(gene)
        for tiss in tiss_list:
            check = False 
            for l in gene_dic[gene]:
                if l[0] == tiss:
                    r2 = l[3]
                    check = True 
            if check != True:
                r2 = "NA"    
            write_list.append(r2)
        fout.write(",".join(write_list)+"\n")
    fout.close() 

def write_outfile4(gene_list,gene_dic):
    print("Writing number of SNPs used result table...")
    fout = gzip.open(out_dir+meta_name+"."+str(alpha)+".nsnps.table.csv.gz",'wb')
    head_list = ["Gene"] + tiss_list
    fout.write(",".join(head_list)+"\n")
    for gene in gene_list:
        write_list = []
        write_list.append(gene)
        for tiss in tiss_list:
            check = False 
            for l in gene_dic[gene]:
                if l[0] == tiss:
                    nsnps = l[4]
                    check = True 
            if check != True:
                nsnps = "NA"    
            write_list.append(nsnps)
        fout.write(",".join(write_list)+"\n")
    fout.close() 

def write_outfile5(gene_list,gene_dic):
    print("Writing number of SNPs in model result table...")
    fout = gzip.open(out_dir+meta_name+"."+str(alpha)+".nsnpsmod.table.csv.gz",'wb')
    head_list = ["Gene"] + tiss_list
    fout.write(",".join(head_list)+"\n")
    for gene in gene_list:
        write_list = []
        write_list.append(gene)
        for tiss in tiss_list:
            check = False 
            for l in gene_dic[gene]:
                if l[0] == tiss:
                    nsnpsmod = l[5]
                    check = True 
            if check != True:
                nsnpsmod = "NA"    
            write_list.append(nsnpsmod)
        fout.write(",".join(write_list)+"\n")
    fout.close() 

def main():
    gl, gd = build_gene_dic()
    write_outfile1(gl,gd)
    write_outfile2(gl,gd)
    write_outfile3(gl,gd)
    write_outfile4(gl,gd)
    write_outfile5(gl,gd)

if (__name__=="__main__"): main()


