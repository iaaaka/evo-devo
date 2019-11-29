#!/bin/bash
#$ -l m_mem_total=22G
#$ -l slt=16
#######################################################################
export PATH=$SGE_O_HOME/bin/bin:/biosoftware/conda/bin:$PATH

cd $SGE_O_HOME/projects/evo.devo/processed/mapping/chimp
index=$SGE_O_HOME/projects/evo.devo/processed/index/Pan_troglodytes.Pan_tro_3.0.dna.toplevel


hisat2 \
	--no-softclip \
	--max-intronlen 1000000 \
	-k 20 \
	--no-unal \
	-q \
	--threads 16 \
	--time \
	-x $index \
	--sra-acc SRR306825 \
	2> SRR306825.log \
	| samtools view -Sb - > SRR306825.bam
	samtools sort -@ 10 -m 1G SRR306825.bam SRR306825
	samtools index SRR306825.bam
