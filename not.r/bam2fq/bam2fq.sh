#!/bin/bash

#usage bam2fq in tmp out

samtools sort -n -@ 3 -O bam -T $2.parts -m 3G -o $2 $1 2> /dev/null
java -Xmx24g -jar /msa2040/raid5/Pavel.Mazin/projects/evo.devo/code/not.r/bam2fq/bam2fq.jar $2 | gzip -1 > $3.gz
rm $2
