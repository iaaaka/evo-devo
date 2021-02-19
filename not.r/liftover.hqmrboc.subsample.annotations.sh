#!/bin/bash
#$ -l slt=2
#$ -l mem_free=6G
#$ -cwd

#######################################################################

export PATH=/biosoftware/conda/bin:$PATH

#cd $SGE_O_HOME/projects/evo.devo/processed/annotation/hqmrboc.subsample/merged/liftover/first
cd $SGE_O_HOME/projects/evo.devo/processed/annotation/hqmrboc.subsample/merged/liftover/second

python $SGE_O_HOME/projects/evo.devo/code/not.r/liftover.hqmrboc.subsample.annotations.py $SGE_TASK_ID
