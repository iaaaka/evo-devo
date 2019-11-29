#!/bin/bash
#$ -l mem_free=15G 
#$ -cwd
#######################################################################
export PATH=/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games:/biosoftware/bin:/biosoftware/conda/bin:/biosoftware/vcftools/bin:/biosoftware/bamtools:/biosoftware/STAR/:/biosoftware/samtools/bin:$PATH

fa=(Gallus_gallus.Galgal4.dna.toplevel.cleanNames.fa Homo_sapiens.GRCh37.73.dna.primary_assembly.cleanNames.fa Macaca_mulatta.MMUL_1.dna.toplevel.cleanNames.fa Monodelphis_domestica.BROADO5.dna.toplevel.cleanNames.fa Mus_musculus.GRCm38.dna.primary_assembly.cleanNames.fa Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.cleanNames.fa Rattus_norvegicus.Rnor_5.0.dna.toplevel.cleanNames.fa)
sp=(chicken human macaque opossum mouse rabbit rat)
f=${fa[$SGE_TASK_ID-1]}
s=${sp[$SGE_TASK_ID-1]}

cd $SGE_O_HOME/data/processed/evo.devo/extract.seq/new.ad.200

cat ${s}200.gtf | python $SGE_O_HOME/projects/evo.devo/code/not.r/extractSeqFromFasta.py $SGE_O_HOME/projects/evo.devo/processed/index/$f > ${s}200.fa
