#!/bin/bash
#######################################################################


reads=/home/Pavel.Mazin/link/projects/evo.devo/processed/fq.from.bams/$species
cd /home/Pavel.Mazin/link/projects/evo.devo/processed/mapping/hisat2.f/$species
index=/home/Pavel.Mazin/link/projects/evo.devo/processed/index/$indexn

ids=( $(ls -1 $reads/*fq.gz | cut -d '/' -f10 | sed s/.fq.gz//) )

i=${ids[$PBS_ARRAYID]}

hisat2 \
	--no-softclip \
	--max-intronlen 1000000 \
	--rna-strandness R \
	--novel-splicesite-outfile $i.splicesites \
	-k 20 \
	--no-unal \
	-q \
	--threads 16 \
	--time \
	-x $index \
	-U $reads/$i.fq.gz \
	2> $i.log \
	| samtools view -Sb - > $i.bam
#	samtools sort -@ 10 -m 1G $i.bam $i
#	samtools index $i.bam
