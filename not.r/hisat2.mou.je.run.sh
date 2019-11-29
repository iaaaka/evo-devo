#!/bin/bash
#$ -l m_mem_total=22G
#$ -l slt=16
#######################################################################
export PATH=$SGE_O_HOME/bin/bin:/biosoftware/conda/bin:$PATH

cd $SGE_O_HOME/projects/evo.devo/processed/mapping/japan.embrio
index=$SGE_O_HOME/projects/evo.devo/processed/index/Mus_musculus.GRCm38.dna.primary_assembly.cleanNames
reads=$SGE_O_HOME/projects/evo.devo/raw/fq.japan.embrio
ids=( $(ls -1 $reads/Mm*fastq.gz | cut -d '/' -f9 | sed s/.fastq.gz//) )

i=${ids[$SGE_TASK_ID-1]}

hisat2 \
	--no-softclip \
	--max-intronlen 1000000 \
	--novel-splicesite-infile $SGE_O_HOME/projects/evo.devo/processed/mapping/junctions/mouse.spMerged.filtered.splicesites \
	--novel-splicesite-outfile $i.splicesites \
	-k 20 \
	--no-unal \
	-q \
	--threads 16 \
	--time \
	-x $index \
	-U $reads/$i.fastq.gz \
	2> $i.log \
	| samtools view -Sb - > $i.bam
	samtools sort -@ 10 -m 1G $i.bam $i
	samtools index $i.bam
