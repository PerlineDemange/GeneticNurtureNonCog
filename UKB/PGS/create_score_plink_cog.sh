#!/bin/bash
#SBATCH -t 03:30:00
#SBATCH -N 1


cd $TMPDIR

cp  /project/ukbaumc/UKBGWAS/genotypes/HM3/ukb_imp.HM3.EUR.v3.* .
cp  /home/pdemange/GeneticNurtureNonCog/LDpred_weights/Cog/* .

ValGf="$TMPDIR/ukb_imp.HM3.EUR.v3" 	## Path to the genotype file in plink binary format to be used in score construction - leave out the plink file extensions (*.bed, *.bim, *.fam).

echo "Creating working directories."
if [ -a ./scores ]
	then
		echo "'scores' directory already exists. Using the existing directory."
	else
		mkdir $TMPDIR/scores
fi


echo "Creating scores."
module load pre2019
module load plink2
plink --bfile $ValGf --score $TMPDIR/COG_CEU_WTS_LDpred-inf.txt 3 4 7 --out $TMPDIR/scores/COG_LDpred-inf_scores 


cp -r $TMPDIR/scores /home/pdemange/UKB/PGS/Cog

