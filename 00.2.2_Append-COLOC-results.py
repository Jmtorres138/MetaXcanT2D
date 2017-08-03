#!/usr/bin/python -O
# Jason Matthew Torres
'''
Script to combine MetaXcan results with COLOC results
Prerequisite: Run COLOC on the GWAS dataset
Usage: python 00.2.2_Append-COLOC-results.py path_to_metaxcan_results path_to_coloc_results path_to_output_directory
'''
# libraries
import sys,os,gzip
import subprocess as sp
import re

#path_to_metaxcan_results = sys.argv[1]
#path_to_coloc_results = sys.argv[2]
#path_to_output_directory = sys.argv[3]
path_to_metaxcan_results = "/group/im-lab/nas40t2/jason/projects/MetaXcan/results/DIAGRAM_ImpG_0.8_gtexV6p/alpha_0.5/"
path_to_coloc_results = "/group/im-lab/nas40t2/jason/projects/MetaXcan/coloc/coloc/results/coloc_all/"
path_to_output_directory = "/group/im-lab/nas40t2/jason/projects/MetaXcan/results/DIAGRAM_ImpG_0.8_gtexV6p/alpha_0.5/combine_coloc/"

# globals

file_list1 = os.listdir(path_to_metaxcan_results)
tiss_list1 = sorted([x.split(".zscores.csv")[0] for x in file_list1 if ".zscores.csv" in x])
tiss_list1 = [re.sub(r"TW_","",x) for x in tiss_list1]

file_list2 = os.listdir(path_to_coloc_results)
tiss_list2 = sorted([x.split(".txt")[0] for x in file_list2 if ".txt" in x])
tiss_list2 = [re.sub(r"DIAGRAM_T2D_TRANS_ETHNIC_eQTL_","",x) for x in tiss_list2]


def build_coloc_dic(coloc_file):
    dic = {}
    # gene_id	P_H0	P_H1	P_H2	P_H3	P_H4
    fin = open(coloc_file,'r')
    fin.readline()
    for line in fin:
        l = line.strip().split()
        gene,ph0,ph1,ph2,ph3,ph4 = l[0],l[1],l[2],l[3],l[4],l[5]
        gene = gene.split(".")[0]
        dic[gene] = [ph0,ph1,ph2,ph3,ph4]
    fin.close()
    return dic

def append_coloc(mtxn_file,coloc_dic):
    file_name = mtxn_file.split("/")[-1]
    fin = open(mtxn_file,'r')
    fout = open(path_to_output_directory+file_name,'w')
    head_list = fin.readline().strip().split()
    head_list = head_list + ["P.H0","P.H1","P.H2","P.H3","P.H4"]
    fout.write(",".join(head_list)+"\n")
    count=0
    for line in fin:
        count+=1
        sys.stdout.write("\r%d"%count)
        sys.stdout.flush()
        l = line.strip().split(",")
        gene = l[0]
        try:
            write_list = l + coloc_dic[gene]
        except:
            write_list = l + ["NA","NA","NA","NA","NA"]
        fout.write(",".join(write_list)+"\n")
    fin.close()
    fout.close()

def main():
    for tiss in tiss_list1:
        if tiss in tiss_list2:
            sys.stdout.write("\nBuilding dictionary for tissue: %s\n" % tiss)
            coloc_file = path_to_coloc_results+"DIAGRAM_T2D_TRANS_ETHNIC_eQTL_"+tiss+".txt"
            coloc_dic = build_coloc_dic(coloc_file)
            mtxn_file = path_to_metaxcan_results+"TW_"+tiss+".zscores.csv"
            sys.stdout.write("\nWriting appended output for tissue: %s\n" % tiss)
            append_coloc(mtxn_file,coloc_dic)
        else:
            sys.stdout.write("\nNo COLOC results available for tissue: %s\n" % tiss)
    sys.stdout.write("\nProcess Complete\n")
if (__name__=="__main__"): main()
