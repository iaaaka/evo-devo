#!/bin/bash
#$ -l mem_free=15G 

export PATH=/biosoftware/conda/bin:$PATH
cd $SGE_O_HOME/projects/evo.devo/processed/cov.bias
cmnd=`head -n $SGE_TASK_ID mouse.in | tail -n 1`
eval "java -jar $SGE_O_HOME/java_src/BamUtil/cov.bias.jar $cmnd 2>> $SGE_TASK_ID.log >> $SGE_TASK_ID.out"
