#!/usr/bin/python -O
# Jason Matthew Torres
'''
Script to combine COLOC results from appended Metaxcan results files (see 00.2.2 script) into summary tables
Usage: python 00.2.4_Build-COLOC-Tables.py meta_name alpha
Example: python 00.2.3_BuildResultsTables.py DIAGRAM1 0.5
'''
# libraries
import sys,os,gzip
import subprocess as sp

# Note: appended results tables are located here:
#/group/im-lab/nas40t2/jason/projects/MetaXcan/results/DIAGRAM_ImpG_0.8_gtexV6p/alpha_0.5/combine_coloc
#meta_name = sys.argv[1]
#alpha = sys.argv[2] # 0.5


# globals
root_dir = "/group/im-lab/nas40t2/jason/projects/MetaXcan/"
res_dir = root_dir + "results/DIAGRAM_ImpG_0.8_gtexV6p/alpha_0.5/combine_coloc/"
out_dir = root_dir + "results/tables/GTExV6p/"
file_list = os.listdir(res_dir)
tiss_list = sorted([x.split(".zscore")[0] for x in file_list if ".zscore" in x])
output_prefix = out_dir+"DIAGRAM_ImpG_0.8_gtexV6p.0.5"
#gene,gene_name,zscore,effect_size,pvalue,VAR_g,pred_perf_R2,pred_perf_p,pred_perf_q,n_snps_used,n_snps_in_cov,n_snps_in_model,P.H0,P.H1,P.H2,P.H3
#,P.H4
def build_gene_dic():
    print("Building gene dictionary...")
    gene_list = []
    gene_dic = {}
    for tiss in tiss_list:
        fin = open(res_dir+tiss+".zscores.csv",'r')
        fin.readline()
        for line in fin:
            l = line.strip().split(',')
            gene,ph0,ph1,ph2,ph3,ph4 = l[1],l[12],l[13],l[14],l[15],l[16]
            gene_list.append(gene)
            try:
                gene_dic[gene].append([tiss,ph0,ph1,ph2,ph3,ph4])
            except:
                gene_dic[gene] = [[tiss,ph0,ph1,ph2,ph3,ph4]]
        fin.close()
    gene_list = sorted(list(set(gene_list)))
    return gene_list, gene_dic

def write_outfile1(gene_list,gene_dic):
    print("Writing probabilty H0 table...")
    fout = gzip.open(output_prefix+".ph0.table.csv.gz",'wb')
    head_list = ["Gene"] + tiss_list
    fout.write(",".join(head_list)+"\n")
    for gene in gene_list:
        write_list = []
        write_list.append(gene)
        for tiss in tiss_list:
            check = False
            for l in gene_dic[gene]:
                if l[0] == tiss:
                    ph0 = l[1]
                    check = True
            if check != True:
                ph0 = "NA"
            write_list.append(ph0)
        fout.write(",".join(write_list)+"\n")
    fout.close()

def write_outfile2(gene_list,gene_dic):
    print("Writing probabilty H1 table...")
    fout = gzip.open(output_prefix+".ph1.table.csv.gz",'wb')
    head_list = ["Gene"] + tiss_list
    fout.write(",".join(head_list)+"\n")
    for gene in gene_list:
        write_list = []
        write_list.append(gene)
        for tiss in tiss_list:
            check = False
            for l in gene_dic[gene]:
                if l[0] == tiss:
                    ph1 = l[2]
                    check = True
            if check != True:
                ph1 = "NA"
            write_list.append(ph1)
        fout.write(",".join(write_list)+"\n")
    fout.close()

def write_outfile3(gene_list,gene_dic):
    print("Writing probabilty H2 table...")
    fout = gzip.open(output_prefix+".ph2.table.csv.gz",'wb')
    head_list = ["Gene"] + tiss_list
    fout.write(",".join(head_list)+"\n")
    for gene in gene_list:
        write_list = []
        write_list.append(gene)
        for tiss in tiss_list:
            check = False
            for l in gene_dic[gene]:
                if l[0] == tiss:
                    ph2 = l[3]
                    check = True
            if check != True:
                ph2 = "NA"
            write_list.append(ph2)
        fout.write(",".join(write_list)+"\n")
    fout.close()

def write_outfile4(gene_list,gene_dic):
    print("Writing probabilty H3 table...")
    fout = gzip.open(output_prefix+".ph3.table.csv.gz",'wb')
    head_list = ["Gene"] + tiss_list
    fout.write(",".join(head_list)+"\n")
    for gene in gene_list:
        write_list = []
        write_list.append(gene)
        for tiss in tiss_list:
            check = False
            for l in gene_dic[gene]:
                if l[0] == tiss:
                    ph3 = l[4]
                    check = True
            if check != True:
                ph3 = "NA"
            write_list.append(ph3)
        fout.write(",".join(write_list)+"\n")
    fout.close()

def write_outfile5(gene_list,gene_dic):
    print("Writing probabilty H4 table...")
    fout = gzip.open(output_prefix+".ph4.table.csv.gz",'wb')
    head_list = ["Gene"] + tiss_list
    fout.write(",".join(head_list)+"\n")
    for gene in gene_list:
        write_list = []
        write_list.append(gene)
        for tiss in tiss_list:
            check = False
            for l in gene_dic[gene]:
                if l[0] == tiss:
                    ph4 = l[5]
                    check = True
            if check != True:
                ph4 = "NA"
            write_list.append(ph4)
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
