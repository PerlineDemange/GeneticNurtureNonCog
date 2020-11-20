# SNPs LD in 1 MB radius.

python ./Program/LDpred.py --coord=COG_CEU.coord --ld_radius=270 --local_ld_file_prefix=COG_CEU_LD --N=223819 --PS=1.00,0.75,0.50,0.30,0.20,0.10,0.05,0.01 --out=COG_CEU_WTS > COG_CEU_LDPred.log
python ./Program/LDpred.py --coord=NONCOG_CEU.coord --ld_radius=270 --local_ld_file_prefix=NON_COG_LD --N=458211 --PS=1.00,0.75,0.50,0.30,0.20,0.10,0.05,0.01 --out=NONCOG_CEU_WTS > NONCOG_CEU_LDPred.log
