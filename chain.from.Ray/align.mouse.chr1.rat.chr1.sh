#!/bin/bash
#BSUB -L /bin/bash
#BSUB -R "rusage[mem=2000]"
#BSUB -R "rusage[tmp=20000]"
#BSUB -M 100000000
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -J mouse.chr1.rat.chr1

CHR1=/home/rmarin/archive/Genomics/Mouse/Ensembl/Softmasked/mouse.chr1.fa
CHR2=/home/rmarin/archive/Genomics/Rat/Ensembl/Softmasked/rat.chr1.fa
ARCHIVE=/home/rmarin/archive/Genomics/Mouse/Multiz/Rat/Chains
ARCHIVE_SCRIPTS=/home/rmarin/archive/Genomics/Mouse/Multiz/Rat/Jobs
MATRIX=/home/rmarin/archive/Genomics/Mouse/Multiz/Rat/Jobs/substitution.matrix
WORKING_DIR=/scratch/local/weekly/rmarin/Genomics/Mouse/Multiz/Rat/Chains
LOG=align.mouse.chr1.rat.chr1.log

# create directory structure
if [ ! -d $WORKING_DIR ]
then
	mkdir -p $WORKING_DIR
fi

if [ ! -d $ARCHIVE ]
then
	ssh $archiveIP 'mkdir -p $ARCHIVE' 
fi

cd $WORKING_DIR
echo "---------------------| `date` | ---------------------" > $LOG

if [ ! -s $ARCHIVE_SCRIPTS/$LOG ]
then
	echo "Running lastz and chaining..." >> $LOG
	$ARCHIVE_SCRIPTS/send.lastz.sh $CHR1 $CHR2 $ARCHIVE no $MATRIX >> $LOG 2>&1 
	echo "---------------------| `date` | ---------------------" >> $LOG
	rsync -av $LOG $cigserver:$ARCHIVE_SCRIPTS/
	rm $LOG &
else
	echo "$ARCHIVE_SCRIPTS/$LOG already exists!" 
fi

date
