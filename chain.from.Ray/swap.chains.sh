#!/bin/bash


sp1=$1
sp2=$2
ARCHIVE=$3

for CHAIN in $ARCHIVE/$sp1.chr*$sp2.chr*.chain
do
	CHR_sp1=`basename $CHAIN | awk '{split($1,v,"."); print v[2]}'`
	CHR_sp2=`basename $CHAIN | awk '{split($1,v,"."); print v[4]}'`
	chainSwap $CHAIN $sp2.$CHR_sp2.$sp1.$CHR_sp1.swap.chain
	echo `basename $CHAIN`" done"
done
