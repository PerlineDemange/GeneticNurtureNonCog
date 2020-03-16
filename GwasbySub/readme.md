Script to run the GWAS-by-subtraction

1. *run_gwas_by_subtraction.sh*: Run the preparatory steps 1-3 for Genomic SEM: munge the summary statistics, calculate the multivariate LDSC and prepare the summary statistics for GWAS. Run the Gwas-by-subtraction without individual SNP effects. 
2. *CogNonCog2.R* and *run_gbs.sh*: Run the Gwas-by-subtraction for every SNPs.  
3. *gwasbysub_processing.txt*: clean summary statistics, identify independent significant SNPs, calculate effective sample size, mean chi2, SNP h2 and LD score intercept. 
