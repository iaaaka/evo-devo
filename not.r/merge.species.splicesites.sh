#!/bin/bash
#$ -l mem_free=10G
#$ -cwd
#$ -soft 
#######################################################################

ins=(mouse chicken rabbit rat human macaque opossum)
#ins=(mou chi rab rat hum mac opo)
fas=(Mus_musculus.GRCm38.dna.primary_assembly.cleanNames.fa Gallus_gallus.Galgal4.dna.toplevel.cleanNames.fa Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.cleanNames.fa Rattus_norvegicus.Rnor_5.0.dna.toplevel.cleanNames.fa Homo_sapiens.GRCh37.73.dna.primary_assembly.cleanNames.fa Macaca_mulatta.MMUL_1.dna.toplevel.cleanNames.fa Monodelphis_domestica.BROADO5.dna.toplevel.cleanNames.fa)
s=${ins[$SGE_TASK_ID-1]}
f=${fas[$SGE_TASK_ID-1]}
cd $HOME/projects/evo.devo/processed/mapping/junctions
#cat ../hisat2.f/$s/*.splicesites | /biosoftware/conda/bin/python $HOME/projects/evo.devo/code/not.r/extract.seq.of.spl.sites.py $HOME/projects/evo.devo/processed/index/$f $s | sort > $s.merged.gff
# run make.junction.liftover.sh
/biosoftware/conda/bin/python $HOME/projects/evo.devo/code/not.r/merge.lo.py $SGE_TASK_ID | gzip > $s.all.sp.merged.gz
