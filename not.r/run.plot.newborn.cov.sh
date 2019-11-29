#!/bin/bash
#$ -l mem_free=16G # not sure it is correct way to request RAM
#$ -l slt=16 #number of cores

R < $SGE_O_HOME/projects/evo.devo/code/newborn.exons.ploc.cov.R --no-save
