# MetaXcanT2D

* Summary statistics from two meta-analyses of GWAS on T2D (DIAGRAM trans-ethnic and GERA-T2D)
* Z-scores for GTEx predictor model SNPs (elastic net; 45 tissue models) imputed to 1000 Genomes phase1 reference panel
* Imputation was performed with ImpG-Summary Software (Version 1.0 July 2013) 
* MetaXcan _zscore_ scheme was used to generate gene-level disease association statistics   

Note: for 04.4_LocusPlotsWindows.Rmd script, you must manually determine the loci corresponding to each group and ajust the locus six array command calls accordingly

Command line aguments used for running 00.2.3 scripts: 
python MetaXcanT2D/00.2.3_BuildResultsTables.gtexV6p.py DIAGRAM_ImpG_0.8_gtexV6p 0.5
python MetaXcanT2D/00.2.3_BuildResultsTables.gtexV6p.py GERA_ImpG_0.8_gtexV6p 0.5 

python 00.2.4_Build-COLOC-Tables.gtexV6p.py


python 00.3.4_build-meta-results-tables.py DIAGRAM-GERA_ImpG_0.80_gtexV6p 0.5
