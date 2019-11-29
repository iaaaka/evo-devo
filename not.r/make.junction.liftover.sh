#!/bin/bash
#$ -l mem_free=10G
#$ -cwd
#$ -soft 
#######################################################################

cd $HOME/projects/evo.devo/processed/mapping/junctions
/biosoftware/conda/bin/python  $HOME/projects/evo.devo/code/not.r/liftover.junctions.py $SGE_TASK_ID


