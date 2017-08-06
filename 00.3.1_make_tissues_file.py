#!/python/bin/python -O
# Jason Matthew Torres
'''
Python script to create tissue file for running meta-analysis scripts
Usage: python 00.3.1_run-meta-jobs.py 

'''


import os

dir = "/group/im-lab/nas40t2/jason/projects/MetaXcan/results/DIAGRAM_ImpG_0.8_gtexV6p/alpha_0.5/"
cwd = "/group/im-lab/nas40t2/jason/projects/MetaXcan/results/meta-analyses/DIAGRAM-GERA_ImpG_0.80_gtexV6p/"
outfile = cwd+"tissues.txt"

l = os.listdir(dir)
l = [x.split(".zscores.csv")[0] for x in l if ".zscores.csv" in x]
l.sort()

fout = open(outfile,'w')
for tiss in l:
    print tiss
    fout.write(tiss+"\n")
fout.close()
