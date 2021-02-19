#!/bin/bash	

CPU=1
SP1=Mouse
SP2=Rat
sp1=`echo $SP1 | awk '{print tolower($0)}'`
sp2=`echo $SP2 | awk '{print tolower($0)}'`
ARCHIVE_FASTA1=/home/rmarin/archive/Genomics/$SP1/Ensembl/Softmasked
ARCHIVE_FASTA2=/home/rmarin/archive/Genomics/$SP2/Ensembl/Softmasked
ARCHIVE=/home/rmarin/archive/Genomics/$SP1/Multiz/$SP2/Chains
ARCHIVE_SCRIPTS=/home/rmarin/archive/Genomics/$SP1/Multiz/$SP2/Jobs
MATRIX=$ARCHIVE_SCRIPTS/substitution.matrix
WORKING_DIR=/scratch/local/weekly/rmarin/Genomics/$SP1/Multiz/$SP2/Chains

for CHR1 in $ARCHIVE_FASTA1/$sp1*chr*.fa
do
	MULTI=`basename $CHR1 | awk '{if($1~/chrUn/) {print "yes"}else{print "no"}}'`
	for CHR2 in $ARCHIVE_FASTA2/$sp2*chr*.fa  
		do
		BASE1=`basename $CHR1 | awk '{print substr($1,1,length($1)-3)}'`
		BASE2=`basename $CHR2 | awk '{print substr($1,1,length($1)-3)}'`
		BASE=$BASE1.$BASE2
		BSUB=align.$BASE.sh

		echo "#!/bin/bash" > $BSUB
		echo "#BSUB -L /bin/bash" >> $BSUB
		echo -e "#BSUB -R \"rusage[mem=2000]\"" >> $BSUB
		echo -e "#BSUB -R \"rusage[tmp=20000]\"" >> $BSUB
		echo "#BSUB -M 100000000" >> $BSUB
		echo "#BSUB -n $CPU" >> $BSUB
		echo -e "#BSUB -R \"span[ptile=$CPU]\"" >> $BSUB
		echo "#BSUB -J $BASE" >> $BSUB
	
		echo >> $BSUB
	
		echo "CHR1=$CHR1" >> $BSUB
		echo "CHR2=$CHR2" >> $BSUB
		echo "ARCHIVE=$ARCHIVE" >> $BSUB
		echo "ARCHIVE_SCRIPTS=$ARCHIVE_SCRIPTS" >> $BSUB
		echo "MATRIX=$MATRIX" >> $BSUB
		echo "WORKING_DIR=$WORKING_DIR" >> $BSUB
		echo -e "LOG=align.$BASE.log" >> $BSUB
		echo >> $BSUB
	
		echo "# create directory structure" >> $BSUB
		echo -e "if [ ! -d \$WORKING_DIR ]" >> $BSUB
		echo "then" >> $BSUB
		echo -e "\tmkdir -p \$WORKING_DIR" >> $BSUB
		echo "fi" >> $BSUB
		echo >> $BSUB

		echo -e "if [ ! -d \$ARCHIVE ]" >> $BSUB
		echo "then" >> $BSUB
		echo -e "\tssh \$archiveIP 'mkdir -p \$ARCHIVE' " >> $BSUB
		echo "fi" >> $BSUB
		echo >> $BSUB

		echo -e "cd \$WORKING_DIR" >> $BSUB
		echo -e "echo \"---------------------| \`date\` | ---------------------\" > \$LOG" >> $BSUB
		echo >> $BSUB

		echo -e "if [ ! -s \$ARCHIVE_SCRIPTS/\$LOG ]" >> $BSUB
		echo "then" >> $BSUB
		echo -e "\techo \"Running lastz and chaining...\" >> \$LOG" >> $BSUB
		echo -e "\t\$ARCHIVE_SCRIPTS/send.lastz.sh \$CHR1 \$CHR2 \$ARCHIVE $MULTI \$MATRIX >> \$LOG 2>&1 " >> $BSUB
		
		echo -e "\techo \"---------------------| \`date\` | ---------------------\" >> \$LOG" >> $BSUB

	        echo -e "\trsync -av \$LOG \$cigserver:\$ARCHIVE_SCRIPTS/" >> $BSUB
        	echo -e "\trm \$LOG &" >> $BSUB
	        echo "else" >> $BSUB
	        echo -e "\techo \"\$ARCHIVE_SCRIPTS/\$LOG already exists!\" " >> $BSUB
        	echo "fi" >> $BSUB
	        echo >> $BSUB
		echo "date" >> $BSUB
	done
done
