#!/bin/bash
#SBATCH -t 35:00:00
#SBATCH -N 1


cd $TMPDIR
#Copying data to scratch for a single node job. Make reading data faster. 
cp $HOME/UKB/GWASbySub/* . # copy content of dir in the current folder (. is telling current folder)

#load necessary modules 
module load pre2019
module load R/3.4.3 

#Execute program located in current folder # no save to not save the envrionment and slave to run as few background processes as possible
{
  Rscript --no-save --slave CogNonCog2.R 1 100000 cognitive1m noncog1m 1000
}&{
  Rscript --no-save --slave CogNonCog2.R 100001 200000 cognitive2m noncog2m 1000
}&{
  Rscript --no-save --slave CogNonCog2.R 200001 300000 cognitive3m noncog3m 1000
}&{
  Rscript --no-save --slave CogNonCog2.R 300001 400000 cognitive4m noncog4m 1000
}&{
  Rscript --no-save --slave CogNonCog2.R 400001 500000 cognitive5m noncog5m 1000
}&{
  Rscript --no-save --slave CogNonCog2.R 500001 600000 cognitive6m noncog6m 1000
}&{
  Rscript --no-save --slave CogNonCog2.R 600001 700000 cognitive7m noncog7m 1000
}&{
  Rscript --no-save --slave CogNonCog2.R 700001 800000 cognitive8m noncog8m 1000
}&{
  Rscript --no-save --slave CogNonCog2.R 800001 900000 cognitive9m noncog9m 1000
}&{
  Rscript --no-save --slave CogNonCog2.R 900001 1000000 cognitive10m noncog10m 1000
}&{
  Rscript --no-save --slave CogNonCog2.R 1000001 1071804 cognitive107m noncog107m 1000

}
wait

#Copy output folder back from scratch
cp *.Rda $HOME/UKB/GWASbySub/output/