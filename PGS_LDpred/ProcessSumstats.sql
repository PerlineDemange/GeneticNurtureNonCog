.mode csv
.separator ","
create table tab1(RSID TEXT, JOINID TEXT, CHR TEXT, POS INTEGER, REF TEXT, ALT TEXT);
create table tab2(CHR TEXT, POS INTEGER, REF TEXT, ALT TEXT, REFFREQ REAL, BETA REAL, SE REAL, P REAL, N INTEGER);
.import ./LDPred/1KG_LDPred_Markermap.csv tab1
.import ./_Table1.dat tab2
.output /dev/shm/_OK.dat
select * from tab1 inner join tab2 on tab1.CHR=tab2.CHR and tab1.POS=tab2.POS and tab1.REF=tab2.REF and tab1.ALT=tab2.ALT;
.output /dev/shm/_INV.dat
select * from tab1 inner join tab2 on tab1.CHR=tab2.CHR and tab1.POS=tab2.POS and tab1.REF=tab2.ALT and tab1.ALT=tab2.REF;
.output stdout
.exit

