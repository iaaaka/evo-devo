#!/bin/bash

FASTA1=$1
FASTA2=$2
ARCHIVE=$3
MULTI=$4
MATRIX=$5
BASE1=`basename $FASTA1 | awk '{print substr($1,1,length($1)-3)}'`
BASE2=`basename $FASTA2 | awk '{print substr($1,1,length($1)-3)}'`
BASE=$BASE1.$BASE2

# Parameters for a rapid test
PARAM_LASTZ="--notransition --step=20 --nogapped --progress --verbosity=10 --runtime"

# Parameters specific to align rat and mouse
# http://hgdownload.cse.ucsc.edu/goldenPath/mm10/vsRn5/README.txt
PARAM_LASTZ="K=3000 L=3000 H=2000 Y=5000 E=55 T=2 O=600 --progress --verbosity=10 --runtime --format=axt Q=$MATRIX"
PARAM_CHAIN="-minScore=5000 -linearGap=medium"

echo "---------------------| `date` | ---------------------" 
echo
echo "------------------           Lastz          -------------------"
if [ $MULTI == yes ]
then
	echo "lastz $FASTA1[multiple] $FASTA2 $PARAM_LASTZ > $BASE.axt "
	lastz $FASTA1[multiple] $FASTA2 $PARAM_LASTZ > $BASE.axt 
else
	echo "lastz $FASTA1 $FASTA2 $PARAM_LASTZ > $BASE.axt "
	lastz $FASTA1 $FASTA2 $PARAM_LASTZ > $BASE.axt 
fi

echo "---------------------| `date` | ---------------------" 
echo 
echo "------------------          axtChain         -------------------"
echo "axtChain -verbose=2 -faQ -faT $BASE.axt $FASTA1 $FASTA2 $BASE.chain $PARAM_CHAIN"
axtChain -verbose=2 -faQ -faT $BASE.axt $FASTA1 $FASTA2 $BASE.chain $PARAM_CHAIN

# copy to archive
echo "---------------------| `date` | ---------------------" 
echo "copying OUTPUT to $ARCHIVE"
rsync -av $BASE.axt $archiveIP:$ARCHIVE/
rsync -av $BASE.chain $archiveIP:$ARCHIVE/

# removing all files
echo "---------------------| `date` | ---------------------" 
echo "Removing all files in local HD"
rm -rf $BASE.* &

