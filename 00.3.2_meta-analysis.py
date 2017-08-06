#!/python/bin/python -O
# Jason Matthew Torres
'''
Python script to run Gene association with GREX predictors
Usage: This script is called from within 00.3.3_meta-analysis.py scirpt, no need to run it on your own

'''
#libraries
import sys, os
import scipy, numpy
from scipy.stats import norm
import math
from directory_functions import *

#globals
cout = sys.stdout.write
cerr = sys.stderr.write

run_dir = "/group/im-lab/nas40t2/jason/projects/MetaXcan/results/meta-analyses/DIAGRAM-GERA_ImpG_0.80_gtexV6p/"
job_dir = run_dir+"jobs/"
if not os.path.exists(job_dir):
    os.makedirs(job_dir)
out_dir = run_dir+"output/"
if not os.path.exists(out_dir):
    os.makedirs(out_dir )

results_dir = "/group/im-lab/nas40t2/jason/projects/MetaXcan/results/"


cohort_file =  str(sys.argv[1])
tissue = str(sys.argv[2])
alpha = str(sys.argv[3])
size_file = str(sys.argv[4])

def phi(x):
    #'Cumulative distribution function for the standard normal distribution'
    return (1.0 + math.erf(x / math.sqrt(2.0))) / 2.0

def rational_approximation(t):

    # Abramowitz and Stegun formula 26.2.23.
    # The absolute value of the error should be less than 4.5 e-4.
    c = [2.515517, 0.802853, 0.010328]
    d = [1.432788, 0.189269, 0.001308]
    numerator = (c[2]*t + c[1])*t + c[0]
    denominator = ((d[2]*t + d[1])*t + d[0])*t + 1.0
    return t - numerator / denominator

def normal_cdf_inverse(p):

    assert p > 0.0 and p < 1

    # See article above for explanation of this section.
    if p < 0.5:
        # F^-1(p) = - G^-1(p)
        return -rational_approximation( math.sqrt(-2.0*math.log(p)) )
    else:
        # F^-1(p) = G^-1(1-p)
        return rational_approximation( math.sqrt(-2.0*math.log(1.0-p)) )

def make_assoc_dic(assoc_file):
    '''
    Reads a MetaXcan association file and returns a dictionary
    with gene key and value list of beta, se, and pvalue
    '''
    fin = open(assoc_file,'r')
    header = fin.readline()
    dic = {}
    for line in fin:
        l = line.strip().split(",")
        gene,zscore,pval = l[1],l[2],l[4]
        if zscore != "NA" and pval != "NA":
            dic[gene] = [zscore,pval]
        else:
            pass
    fin.close()
    return dic

def meta_sample_size(size_list,pval_list,sign_list):
    '''
    This function allows the user to meta-analyzie a single test (e.g. gene, snp)
    with lists of effective study sample sizes, pvals, and signs across studies
    '''
    if len(size_list) == len(pval_list) == len(sign_list):
        length = len(size_list)
    else:
        raise Exception("The provided lists have differing lengths")
    Z_list = []
    w_list = []
    for i in range(0, length):
        N_i = float(size_list[i])
        P_i = float(pval_list[i])
        delta_i = float(sign_list[i])
        Z_i = normal_cdf_inverse(P_i/float(2)) * delta_i
        w_i = math.sqrt(N_i)
        Z_list.append(Z_i)
        w_list.append(w_i)
    numerator = 0
    for i in range(0,length):
        numerator += float(Z_list[i]) * float(w_list[i])
    denominator = 0
    for i in range(0,length):
        denominator += float(w_list[i]) ** 2
    denominator = math.sqrt(denominator)
    Z = -1 * float(numerator)/float(denominator) # The addition of the -1 is a bit of a hack to get the signs to agree, NEED TO INVESTIGATE THIS FURTHER
    # Believe that there was a typo in the METAL paper (Willer et al.); will go with alternative formula
    P = 2 * (1 - norm.cdf(abs(-Z)))
    return  str(Z), str(P)

def meta_sample_size_version2(size_list,z_list):
    '''
    This function allows the user to meta-analyzie a single test (e.g. gene, snp)
    with lists of effective study sample sizes, pvals, and signs across studies
    '''
    if len(size_list) == len(z_list):
        length = len(size_list)
    else:
        raise Exception("The provided lists have differing lengths")
    Z_list = []
    w_list = []
    for i in range(0, length):
        N_i = float(size_list[i])
        Z_i = z_list[i]
        w_i = math.sqrt(N_i)
        Z_list.append(Z_i)
        w_list.append(w_i)
    numerator = 0
    for i in range(0,length):
        numerator += float(Z_list[i]) * float(w_list[i])
    denominator = 0
    for i in range(0,length):
        denominator += float(w_list[i]) ** 2
    denominator = math.sqrt(denominator)
    Z = float(numerator)/float(denominator)
    P = 2 * (1 - norm.cdf(abs(-Z)))
    return  str(Z), str(P)



def meta_inverse_var(beta_list,se_list):
    '''
    This function allows the user to meta-analyzie a single test (e.g. gene, snp)
    with lists of betas, and standard errors across studies
    NOTE: Will not use this approach for MetaXcan, for simplicity
    will focus on sample size method above
    '''
    if len(beta_list) == len(se_list):
        length = len(beta_list)
    else:
        raise Exception("The provided lists have differing lengths")
    w_list = []
    for i in range(0, length):
        B_i = float(beta_list[i])
        se_i = float(se_list[i])
        w_i = float(1) / (float(se_i) ** 2)
        w_list.append(w_i)
    w_sum = 0
    for i in range(0,length):
        w_sum += float(w_list[i])
    beta_numer = 0
    for i in range(0,length):
        beta_numer += float(beta_list[i]) * float(w_list[i])
    SE = math.sqrt(float(1)/float(w_sum))
    Beta = float(beta_numer) / float(w_sum)
    Z = float(Beta)/float(SE)
    P = 2 * (1 - norm.cdf(abs(-Z)))
    return str(Beta), str(SE), str(Z), str(P)

def make_meta_file(cohort_list, size_list, tissue, alpha):
    global assoc_dir, out_dir
    dict_list = []
    coh_len = len(cohort_list)
    cout("Creating Dictionaries for Tissue: %s\tAlpha: %s\n" % (tissue, str(alpha)))
    for coh in cohort_list:
        #assoc_file_name  =  assoc_dir+coh+"/"+tissue+"_"+str(alpha)+"_assoc.txt"
        assoc_file_name = results_dir + coh+"/alpha_"+str(alpha)+"/"+tissue+".zscores.csv"
        d_coh = make_assoc_dic(assoc_file_name)
        dict_list.append(d_coh)
    gene_list = []
    for dic in dict_list:
        glist = dic.keys()
        gene_list = gene_list + glist
    cout("Determing Gene Set for Tissue: %s\n" % tissue)
    gene_set = []
    for gene in gene_list:
        if gene_list.count(gene) == coh_len:
            gene_set.append(gene)
    gene_set = list(set(gene_set))
    cout("Writing Output file...\n")
    #fout = open(out_dir+tissue+"_"+str(alpha)+"_meta.txt", 'w')
    if not os.path.exists(out_dir+"alpha_"+str(alpha)+"/"):
        os.makedirs(out_dir+"alpha_"+str(alpha)+"/")
    fout = open(out_dir+"alpha_"+str(alpha)+"/"+tissue+".meta.txt",'w')#
    #head_list = ["Gene","ss_Z","ss_P","iv_Beta","iv_SE","iv_Z","iv_P"]
    head_list = ["Gene","ss.zscore","ss.pvalue"]
    for coh in cohort_list:
        zname = coh+".zscore"
        pname = coh+".pvalue"
        append_list = [zname,pname]
        head_list = head_list + append_list
    head_list.append("more.significant")
    fout.write("\t".join(head_list)+"\n")
    for gene in gene_set:
        pval_list = []
        sign_list = []
        z_list = []
        #beta_list = []
        #se_list = []
        for dic in dict_list:
            #b_i, se_i, p_i, sign_i = dic[gene][0], dic[gene][1], dic[gene][2], numpy.sign(float(dic[gene][0]))
            p_i, z_i, sign_i = dic[gene][1], dic[gene][0] ,numpy.sign(float(dic[gene][0]))
            #beta_list.append(b_i)
            #se_list.append(se_i)
            pval_list.append(p_i)
            z_list.append(z_i)
            sign_list.append(sign_i)
        #iv_Beta, iv_SE, iv_Z, iv_P = meta_inverse_var(beta_list,se_list)
        #ss_Z, ss_P = meta_sample_size(size_list,pval_list,sign_list)
        ss_Z, ss_P = meta_sample_size_version2(size_list,z_list)
        #write_list = [gene,ss_Z,ss_P,iv_Beta,iv_SE,iv_Z,iv_P]
        if float(ss_P) < float(pval_list[0]) and float(ss_P) < float(pval_list[1]):
            moresig = "TRUE"
        else:
            moresig = "FALSE"
        write_list =[gene,ss_Z,ss_P]
        for i in range(0,len(z_list)):
            append_list = [z_list[i],pval_list[i]]
            write_list = write_list + append_list
        write_list.append(moresig)
        fout.write("\t".join(write_list)+"\n")
    fout.close()



if __name__ == "__main__":

    cohort_list = list_from_file(cohort_file,1)
    size_list = list_from_file(size_file,1)

    make_meta_file(cohort_list, size_list, tissue, alpha)
