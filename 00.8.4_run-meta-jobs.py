#!/python/bin/python -O
# Jason Matthew Torres
'''
Python script to meta-analyses of MetaXcan results from DIAGRAM and GERA 
'''
#libraries 
import sys, os
import subprocess as sp
from JTfunctions import * 
#globals 
cout = sys.stdout.write
cerr = sys.stderr.write

meta_dir = "/group/im-lab/nas40t2/jason/projects/MetaXcan/results/meta-analyses/DIAGRAM-GERA_ImpG_0.80_gtexV6p/"
job_dir = meta_dir+"jobs/"


cohort_file =  str(sys.argv[1]) # file with unique name of cohort per line ie T2D NBS 58C
tiss_file = str(sys.argv[2]) # file with name of tissue per line
alpha_file = str(sys.argv[3]) # file with alpha per line ie 0.5 1 
size_file = str(sys.argv[4])

def run_job_script(tissue,alpha):
    '''
    Generate shell script to submit to cluster
    '''
    global cohort_file, size_file, job_dir 
    tissue, alpha = str(tissue), str(alpha) 
    cout("Writing job script for Tissue: %s\tAlpha: %s\n" % (tissue,alpha))
    fout = open("temp_job.sh", 'w')
    call = "python  %sJTmeta.py %s %s %s %s" % (meta_dir, meta_dir+cohort_file, tissue, alpha, meta_dir+size_file) 
            
    script = '''
#!/bin/bash
#PBS -N meta_%s_%s 
#PBS -S /bin/bash
#PBS -l mem=10gb
#PBS -o %sjob_%s_%s_meta.out
#PBS -e %sjob_%s_%s_meta.err

module load python/2.7.6

%s   
    ''' % (tissue,alpha,job_dir,tissue,alpha,job_dir,tissue,alpha,call)
    
    fout.write(script)
    fout.close()
    cout("Running job....\n")
    job = ["sh","temp_job.sh"] #job = call.strip().split()
    #job = ["qsub","temp_job.sh"] #job = call.strip().split()
    sp.check_call(job)

    #os.remove("temp_job.sh")
        
if (__name__ == "__main__"):
    
    tiss_list = list_from_file(tiss_file,1) 
    alpha_list = list_from_file(alpha_file,1) 

    for tiss in tiss_list:
        for alp in alpha_list:
            run_job_script(tiss,alp)

                
    
    
