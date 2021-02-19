#!/bin/bash
#PBS -N hqmrboc.subsampling.cluster.R
#PBS -l walltime=30:00:00
#PBS -l nodes=1:ppn=16
#PBS -l mem=26gb
#######################################################################

cd /gss/mazin/projects/evo.devo
module load R/4.0.3
 
R < code/hqmrboc.subsampling.cluster.R --no-save