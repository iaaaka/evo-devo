#!/bin/bash
#BSUB -L /bin/bash
#BSUB -R "rusage[mem=10000]"
#BSUB -R "rusage[tmp=20000]"
#BSUB -M 100000000
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -J mouse.rat.over.chain


SP1=Mouse
SP2=Rat
sp1=`echo $SP1 | awk '{print tolower($0)}'`
sp2=`echo $SP2 | awk '{print tolower($0)}'`
ARCHIVE_FASTA1=/home/rmarin/archive/Genomics/$SP1/Ensembl/Softmasked
ARCHIVE_FASTA2=/home/rmarin/archive/Genomics/$SP2/Ensembl/Softmasked
ARCHIVE=/home/rmarin/archive/Genomics/$SP1/Multiz/$SP2/Chains
ARCHIVE_SCRIPTS=/home/rmarin/archive/Genomics/$SP1/Multiz/$SP2/Jobs
WORKING_DIR=/scratch/local/weekly/rmarin/Genomics/$SP1/Multiz/$SP2/Chains
BASE=$sp1"To"$sp2
LOG=chain.liftover.$BASE.log


# create directory structure
if [ ! -d $WORKING_DIR ]
then
	mkdir -p $WORKING_DIR
fi

cd $WORKING_DIR
echo "---------------------| `date` | ---------------------" > $LOG

if [ ! -s $ARCHIVE_SCRIPTS/$LOG ]
then
	echo "chainMergeSort..." >> $LOG
	chainMergeSort $ARCHIVE/$sp1.chr*$sp2.chr*.chain > $BASE.chain.raw 2>> $LOG 
	echo "---------------------| `date` | ---------------------" >> $LOG
	echo "chainNet..." >> $LOG
	chainNet $BASE.chain.raw $ARCHIVE_FASTA1/$sp1.sizes $ARCHIVE_FASTA2/$sp2.sizes $BASE.net /dev/null >> $LOG 2>&1
	echo "---------------------| `date` | ---------------------" >> $LOG
	echo "netChainSubset..." >> $LOG
	netChainSubset $BASE.net $BASE.chain.raw $BASE.over.chain >> $LOG 2>&1
	echo "---------------------| `date` | ---------------------" >> $LOG
	rsync -av $BASE.chain.raw $cigserver:$ARCHIVE/
	rsync -av $BASE.over.chain $cigserver:$ARCHIVE/
	rsync -av $BASE.net $cigserver:$ARCHIVE/
	rsync -av $LOG $cigserver:$ARCHIVE_SCRIPTS/
	rm -rf $BASE* $LOG &
else
	echo "$ARCHIVE_SCRIPTS/$LOG already exists!" 
fi

