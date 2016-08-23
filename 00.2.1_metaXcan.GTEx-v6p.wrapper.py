#!/usr/bin/python -O 
# Jason Matthew Torres
'''
Wrapper for MetaXcan SNP check and MetaXcan gene associaton test (i.e. scripts M03,M04) 
Usage: python JTmetaXcan_ImpG_gtexV6p.wrapper.py  
'''
# libraries 
import sys,os
import subprocess as sp 

# globals 
met_dir = "/group/im-lab/nas40t2/jason/projects/MetaXcan/"
pydir = met_dir + "MetaXcan-Aug2016/software/"
pyfile1 = pydir + "M03_betas.py"
pyfile2 = pydir + "M04_zscores.py"
temp_dir = met_dir + "temp/"
db_dir = "/group/im-lab/nas40t2/abarbeira/Projects/model_dbs/v6p/0_5/"
cov_dir = "/group/im-lab/nas40t2/abarbeira/Projects/MetaXcanResults/intermediate/cov_tgf_eur/v6p/0_5/"
meta_analysis_dir = met_dir+ "meta_files/ImpGformat/merged_output_files/impG_DIAGRAM/gwas_file_v6p/" 
meta_name =   "DIAGRAM_ImpG_0.8_gtexV6p" # "DIAGRAM_ImpG_0.8_gtexV6p"

# functions 
def run_jobs(tislist,alpha,meta_name,threshold=None):
    for tissue in tislist:
        arg_list1 = ["python",
            pyfile1,
            "--weight_db_path",
            db_dir+str(tissue)+"_"+str(alpha)+".db",
            "--gwas_folder",
            meta_analysis_dir,
            "--output_folder",
            "gtex_metaXcan_betas/"+str(alpha)+"/"+tissue+"/"+meta_name,
            "--beta_zscore_column", "Z",
            "--a1_column", "A1", "--a2_column", "A2",
            "--snp_column", "SNP", 
            "--compressed"]
        arg_line1 = " ".join(arg_list1)
        arg_list2 = ["python",
            pyfile2,
            "--weight_db_path",
            db_dir+str(tissue)+"_"+str(alpha)+".db",
            "--covariance",
            #met_dir+"gtex_covariance/unscaled/"+str(alpha)+"/"+tissue+"/covariance.txt.gz",
            cov_dir+tissue+".txt.gz",
            "--beta_folder",
            "gtex_metaXcan_betas/"+str(alpha)+"/"+tissue+"/"+meta_name,
            "--zscore_scheme", "beta_z_and_ref",
            "--normalization_scheme", "none", 
            "--throw"]
        if threshold:
            out_dir = "results/"+meta_name+"/alpha_"+str(alpha)+"/"+str(threshold)+"/"
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            out_file = out_dir + tissue +".zscores.csv"
            arg_list2 = arg_list2 + ["--output_file",out_file]
        else:
            out_dir = "results/"+meta_name+"/alpha_"+str(alpha)+"/"
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            out_file = out_dir + tissue+".zscores.csv"
            arg_list2 = arg_list2 + ["--output_file",out_file]
        arg_line2 = " ".join(arg_list2)
        fname = temp_dir+tissue+".temp.sh"
        print "Writing script file: %s" % fname
        fout = open(fname,'w')
        script = '''
                 #!/bin/bash                                                                                                                              

                 #PBS -N %s                                                                                                                            
                 #PBS -S /bin/bash                                                                                                                        
                 #PBS -l walltime=24:00:00                                                                                                                
                 #PBS -l nodes=1:ppn=10                                                                                                                   
                 #PBS -l mem=2gb                                                                                                                         
                 #PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log                                                                                            
                 #PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err                                                                                            

                 cd $PBS_O_WORKDIR

                 module load python/2.7.9

                 %s
                 %s 
                 ''' % ("MtsXn"+"."+meta_name+"."+tissue+"."+str(alpha), arg_line1, arg_line2)
        fout.write(script+"\n")
        fout.close() 
        command = ["qsub",fname] 
        sp.check_call(command)

def main():
    dbfiles = os.listdir(db_dir)
    dbfiles_sub = [x for x in dbfiles if x.endswith("_0.5.db")]
    tislist = [x.split("_0.5.db")[0] for x in dbfiles_sub] 
    #tislist = [x for x in tislist if "TS_" not in x]
    #tislist = [x for x in tislist if "Organ_" not in x]
    tislist = ["TW_Muscle_Skeletal"]#["TW_WholeBloodDGN"]#
    alpha_list = ["0.5"]# "1"
    for alpha in alpha_list:
        run_jobs(tislist,alpha,meta_name,threshold=None)

if (__name__=="__main__"): main()
