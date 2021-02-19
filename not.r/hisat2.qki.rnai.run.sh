#!/bin/bash
#PBS -N map_qki_rnai
#PBS -l walltime=3:00:00
#PBS -l nodes=1:ppn=5
#PBS -l mem=15gb
#######################################################################
export PATH=/gss/mazin/bin/bin:$PATH

cd /gss/mazin/projects/evo.devo/processed/mapping/QKI.RNAi
index=/gss/mazin/projects/evo.devo/processed/index/Homo_sapiens.GRCh37.73.dna.primary_assembly.cleanNames

ids=(SRR3469452 SRR3469453 SRR3469516 SRR3469517 SRR3469440 SRR3469441 SRR3469464 SRR3469465)

i=${ids[$PBS_ARRAYID-1]}

hisat2 \
	--no-softclip \
	--max-intronlen 1000000 \
	--novel-splicesite-infile $SGE_O_HOME/projects/evo.devo/processed/mapping/junctions/human.spMerged.filtered.splicesites \
	-k 20 \
	--no-unal \
	-q \
	--threads 5 \
	--time \
	-x $index \
	--sra-acc $i \
	2> bams/$i.log \
	| samtools view -Sb - > bams/$i.bam
	samtools sort -@ 5 -m 1G bams/$i.bam bams/$i
	samtools index bams/$i.bam

java -Xmx10g -jar /gss/mazin/bin/sajr.jar count_reads sajr.config \
	-batch_in=bams/$i.bam \
	-ann_in=../../annotation/hqmrboc.subsample/merged/human.sajr \
	-batch_out=sajr/$i \
	-use_mult=false >> sajr.log

