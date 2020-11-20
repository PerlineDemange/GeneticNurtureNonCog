# NOTES.
# - Requires build 37 and + stand data for now in all cases, if half of the data are left after matching please check strand of sumstats.
# - For LDpred needed chr pos ref alt reffrq info rs pval effalt=beta N
# - Minimal needed from summary stats : chr,pos,ref(nea),alt(ea),beta,se,p,N
# - RS_IDs taken from genotype map files as markernames to make sure thay are consistent.
# - No X chromosome (should have "X" as a chromsome ID otherwise).

# Needed order: "CHR,POS,REF,ALT,FRQREF,BETA-ALT,SE,P,N"
# Header sumstats: SNP	CHR	BP	A1	A2	MAF	Beta	SE	Z	P.

INSUMSTATS=./Sumstats/Cog_GWAS_110320.txt
OUTLDPRED=./COG_11032020_ldpred_in.dat

cat $INSUMSTATS | sed '1,1d' | awk '{print $2" "$3" "$4" "$5" 0.50 "$7" "$8" "$10" 223819"}' > _Sum1.dat

# Sanity check p-values ao.
# Note the sumstats give maf, not effect allele frequency. It seems that all SNPs are connected to REF=A1 and ALT=A2.
# thus MAF is not related the A1 alelle. Ignoring frequencies and assuming palindromic are correct in Meta.
# Fixing allele frequencies to 0.50 in above statement. For now considering A2 to be effect allele.

awk '{if ($8>=0 && $8<=1) print $0}' _Sum1.dat | tr ' ' ',' > _Table1.dat

sqlite3 /dev/shm/_temp.db < ProcessSumstats.sql
rm /dev/shm/_temp.db

# Joined table
#RSID	JOINID	CHR	BP	REF	ALT	CHR	BP	REF*/ALT	ALT/REF*	AFREF*	BETA	SE	P	N	
#1	2	3	4	5	6	7	8	9		10		11	12	13	14	15	

echo "CHR POS SNP_ID REF ALT REF_FRQ PVAL BETA SE N" > $OUTLDPRED

sed 's/,/ /g' /dev/shm/_OK.dat | awk '{print "chr"$3" "$4" "$1" "$5" "$6" "$11" "$14" "$12" "$13" "$15}' >> $OUTLDPRED
sed 's/,/ /g' /dev/shm/_INV.dat | awk '{print "chr"$3" "$4" "$1" "$5" "$6" "(1-$11)" "$14" "(-1*$12)" "$13" "$15}' >> $OUTLDPRED

rm /dev/shm/_*
rm _*

# Coordinating.

SUMST=$OUTLDPRED
OUT=COG_CEU.coord
GT=./LDPred/EUROPEAN_1KG_P3V5A_B37_FixJJ_C1-22_LDPred
echo $SUMST 
cut -d' ' -f3 $SUMST | sed '1,1d' > _Needed_snps.dat
plink2 --bfile $GT --extract _Needed_snps.dat --maf 0.01 --make-bed --out _GTdata
ldpred coord --gf=_GTdata --ssf-format=LDPRED --max-freq-discrep=0.50 --ssf=$SUMST --out=$OUT > COG_coord.log
rm _*

