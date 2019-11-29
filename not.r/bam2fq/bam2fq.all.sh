#!/bin/bash
#######################################################################
export PATH=/msa2040/raid5/Pavel.Mazin/bin/bin:$PATH

# mouse
#qsub -l nodes=1:ppn=4,mem=24gb -l walltime=48:00:00 -l feature=hpmsa -t 0-476 bam2fq.all.sh
#cd /msa2040/raid5/Pavel.Mazin/projects/evo.devo/processed/fq.from.bams/mouse
#sams=( $(ls -1 /msa2040/raid5/Pavel.Mazin/projects/evo.devo/raw/bams.from.lausanne/*Mouse*bam | cut -d '/' -f9 | sed s/.sorted.bam//) )
#sam=${sams[$PBS_ARRAYID]}
#/msa2040/raid5/Pavel.Mazin/projects/evo.devo/code/not.r/bam2fq/bam2fq.sh /msa2040/raid5/Pavel.Mazin/projects/evo.devo/raw/bams.from.lausanne/$sam.sorted.bam tmp/$sam.bam $sam.fq

# chicken
#qsub -l nodes=1:ppn=4,mem=24gb -l walltime=48:00:00 -l feature=hpmsa -t 0-258 bam2fq.all.sh
#cd /msa2040/raid5/Pavel.Mazin/projects/evo.devo/processed/fq.from.bams/chicken
#sams=( $(ls -1 /msa2040/raid5/Pavel.Mazin/projects/evo.devo/raw/bams.from.lausanne/*Chicken*bam | cut -d '/' -f9 | sed s/.sorted.bam//) )
#sam=${sams[$PBS_ARRAYID]}
#/msa2040/raid5/Pavel.Mazin/projects/evo.devo/code/not.r/bam2fq/bam2fq.sh /msa2040/raid5/Pavel.Mazin/projects/evo.devo/raw/bams.from.lausanne/$sam.sorted.bam tmp/$sam.bam $sam.fq

# Rabbit
#qsub -l nodes=1:ppn=4,mem=24gb -l walltime=48:00:00 -l feature=hpmsa -t 0-436 bam2fq.all.sh 
#cd /msa2040/raid5/Pavel.Mazin/projects/evo.devo/processed/fq.from.bams/rabbit
#sams=( $(ls -1 /msa2040/raid5/Pavel.Mazin/projects/evo.devo/raw/bams.from.lausanne/*Rabbit*bam | cut -d '/' -f9 | sed s/.sorted.bam//) )
#sam=${sams[$PBS_ARRAYID]}
#/msa2040/raid5/Pavel.Mazin/projects/evo.devo/code/not.r/bam2fq/bam2fq.sh /msa2040/raid5/Pavel.Mazin/projects/evo.devo/raw/bams.from.lausanne/$sam.sorted.bam tmp/$sam.bam $sam.fq

# Rat
#qsub -l nodes=1:ppn=4,mem=24gb -l walltime=48:00:00 -l feature=hpmsa -t 0-455 bam2fq.all.sh 
#cd /msa2040/raid5/Pavel.Mazin/projects/evo.devo/processed/fq.from.bams/rat
#sams=( $(ls -1 /msa2040/raid5/Pavel.Mazin/projects/evo.devo/raw/bams.from.lausanne/*Rat*bam | cut -d '/' -f9 | sed s/.sorted.bam//) )
#sam=${sams[$PBS_ARRAYID]}
#/msa2040/raid5/Pavel.Mazin/projects/evo.devo/code/not.r/bam2fq/bam2fq.sh /msa2040/raid5/Pavel.Mazin/projects/evo.devo/raw/bams.from.lausanne/$sam.sorted.bam tmp/$sam.bam $sam.fq

#Human
#qsub -l nodes=1:ppn=4,mem=24gb -l walltime=48:00:00 -l feature=hpmsa -t 0-363 bam2fq.all.sh
#cd /msa2040/raid5/Pavel.Mazin/projects/evo.devo/processed/fq.from.bams/human
#sams=( $(ls -1 /msa2040/raid5/Pavel.Mazin/projects/evo.devo/raw/bams.from.lausanne/*Human*bam | cut -d '/' -f9 | sed s/.sorted.bam//) )
#sam=${sams[$PBS_ARRAYID]}
#/msa2040/raid5/Pavel.Mazin/projects/evo.devo/code/not.r/bam2fq/bam2fq.sh /msa2040/raid5/Pavel.Mazin/projects/evo.devo/raw/bams.from.lausanne/$sam.sorted.bam tmp/$sam.bam $sam.fq

#Macaque
#qsub -l nodes=1:ppn=4,mem=24gb -l walltime=48:00:00 -l feature=hpmsa -t 0-184 bam2fq.all.sh
#cd /msa2040/raid5/Pavel.Mazin/projects/evo.devo/processed/fq.from.bams/macaque
#sams=( $(ls -1 /msa2040/raid5/Pavel.Mazin/projects/evo.devo/raw/bams.from.lausanne/*Macaque*bam | cut -d '/' -f9 | sed s/.sorted.bam//) )
#sam=${sams[$PBS_ARRAYID]}
#/msa2040/raid5/Pavel.Mazin/projects/evo.devo/code/not.r/bam2fq/bam2fq.sh /msa2040/raid5/Pavel.Mazin/projects/evo.devo/raw/bams.from.lausanne/$sam.sorted.bam tmp/$sam.bam $sam.fq

#Opossum
#qsub -l nodes=1:ppn=4,mem=24gb -l walltime=48:00:00 -l feature=hpmsa -t 0-274 bam2fq.all.sh
cd /msa2040/raid5/Pavel.Mazin/projects/evo.devo/processed/fq.from.bams/opossum
sams=( $(ls -1 /msa2040/raid5/Pavel.Mazin/projects/evo.devo/raw/bams.from.lausanne/*Opossum*bam | cut -d '/' -f9 | sed s/.sorted.bam//) )
sam=${sams[$PBS_ARRAYID]}
/msa2040/raid5/Pavel.Mazin/projects/evo.devo/code/not.r/bam2fq/bam2fq.sh /msa2040/raid5/Pavel.Mazin/projects/evo.devo/raw/bams.from.lausanne/$sam.sorted.bam tmp/$sam.bam $sam.fq
