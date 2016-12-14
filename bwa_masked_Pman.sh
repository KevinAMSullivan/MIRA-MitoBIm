#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N BWA_Whole
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q R2D2
#$ -pe fill 19 
#$ -P communitycluster

#This file will process raw Illumina data using Trimmomatic.  This will be followed by mapping to a reference genome to create a new genome assembly.

BASEDIR=/lustre/scratch/kevsulli/Mitobim/Whole_Run
mkdir $BASEDIR/Output
WORKDIR=$BASEDIR/Output

cd $WORKDIR

THREADS=19  #Line 9 sets this up to run 20 processors.  If you want fewer, make sure to change that line as well as this one.

#for RAW_READ_FILE in $BASEDIR/data_links/*.fastq.gz;
#do gunzip -c $BASEDIR/data_links/$RAW_READ_FILE >$RAW_READ_FILE.fastq &;	\
#done 
#wait


RAW_READS_HOME=/lustre/scratch/kevsulli/Mitobim/Output/data_links   #the location of your raw data
#PROCESSED_READS_HOME=$BASEDIR/data_processed	#the location of your processed data
QUALITY_INFO=$RAW_READS_HOME/QC_files	#Where the quality stats will be saved
SUPPORT_FILES=$BASEDIR/support_files	#where the support files like the adapter sequences will be located.

######
#set up alias' for major programs
######
BWA_HOME=/lustre/work/apps/bwa-0.7.12
SAMTOOLS_HOME=/lustre/work/apps/samtools-1.2
SAMTOOLS1_8_HOME=/lustre/work/apps/samtools-0.1.18
PICARD_HOME=/lustre/work/apps/picard-tools-1.91
BCFTOOLS_HOME=/lustre/work/apps/samtools-0.1.18/bcftools
RAY_SOFTWARE=/lustre/work/daray/software
TRIM_HOME=/lustre/work/apps/Trimmomatic-0.27
FASTX_HOME=/lustre/work/apps/fastx_toolkit-0.0.14/bin
VCFTOOLS_HOME=/lustre/work/daray/software/vcftools_0.1.12b/bin
BEDTOOLS_HOME=/lustre/work/apps/bedtools-2.17.0/bin


######
#Set up insert size.  This will be specific to the insert size for your particular taxon's library.
######
#insSize=1000  

######
#make sure your genome file has no blank lines  - ALREADY DONE, NOT BEING USED HERE
######
#sed '/^$/d' $REF_HOME/$refgenome >tempGenome
#cp tempGenome $REF_HOME/$refgenome".clean"

#echo "spaces" |  mailx -s "spaces" kev.am.sullivan@gmail.com


	

#######
###!!!!!!! There are comments after the "\" below.  This won't work with them present.  Make sure to get rid of anything after \ on all lines.
#######

for RAW_READ_FILE in $RAW_READS_HOME/104488_S2_L001_R1_001.fastq
do
	ABBREV=$(basename $RAW_READ_FILE _L001_R1_001.fastq) #This will be the name you use to process your files.  You will need to change the RP1, RP2, RU1, and RU2 slots just below here accordingly.

REFGENOME=$ABBREV-aMon-mt-genome_out_AllStrains.unpadded.fasta  	# your reference genome for the assembly
REF_HOME=$BASEDIR/Genomes	#the location of your reference genome

################################################################################
# Map reads to genome with BWA
#~~~~~~~~~~~

#[1a] Use bwa to index the genome  
$BWA_HOME/bwa index \
		-a is \
		$REF_HOME/$REFGENOME


		
#===================
# [1b] Map the reads to the genome
	$BWA_HOME/bwa mem 			\
	-M	\
		$REF_HOME/$REFGENOME			\
       		$RAW_READS_HOME/$ABBREV"_L001_R1_001.fastq"	\
       		$RAW_READS_HOME/$ABBREV"_L001_R2_001.fastq"	\
	-t $THREADS	\
	> $ABBREV"_aln-pe.sam" 	
		
						
echo $ABBREV"_R1&2_map" |  mailx -s $ABBREV"_R1&2_map" kev.am.sullivan@gmail.com	


	
#===================
# [3] use sampe and SAMtools to create bam files of the mapped reads
#create sam file from paired mapped reads
	$BWA_HOME/bwa sampe 				\
	-f $ABBREV"_SAMPE.sam" 	\
		$REF_HOME/$REFGENOME 				\
		$ABBREV"_aln-pe.sam" 		\
		$ABBREV"_L001_R1_001.fastq"			\
		$ABBREV"_L001_R2_001.fastq"			
		


#convert paired sam file to bam
	$SAMTOOLS_HOME/samtools view 			\
		-Sb 					\
		-o $ABBREV"_aln-pe.bam" 	\
		$ABBREV"_aln-pe.sam" 		



#Not sure what this does - sort bam file?
	$SAMTOOLS_HOME/samtools view 			\
		-F 4 					\
		-q 20 					\
		-b						\
		-o $ABBREV"_R3.bam" 	\
		$ABBREV"_aln-pe.bam"  

	
	$SAMTOOLS_HOME/samtools sort 			\
		$ABBREV"_R3.bam"	\
		-@ $THREADS 		\
		$ABBREV"_R3_sorted"

#### for samtools v1.2		
	$SAMTOOLS_HOME/samtools sort 			\
		-o $ABBREV"_R3_sorted.bam"	\
		-O bam	\
		-T $ABBREV"_sorted" \
		-@ $THREADS 					\
		$ABBREV"_aln-pe.bam" 		

	
		
#===================
# [4] remove sequencing duplicates from the sorted bam file w/ PICARD	
	java 						\
        -Xmx24g 				\
		-Djava.io.tmpdir=tmp 			\
		-jar $PICARD_HOME/MarkDuplicates.jar 	\
        	I=$ABBREV"_R3_sorted.bam" 	\
       		O=$ABBREV"_R3_noDup.bam"	\
        	M=$ABBREV"_R3_dupMetric.out" 	\
        	REMOVE_DUPLICATES=true 			\
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 	\
		VALIDATION_STRINGENCY=SILENT 		\
		ASSUME_SORTED=TRUE 			\
		TMP_DIR=tmp

#==================
#[5] create pileup from noDup.bam

$SAMTOOLS1_8_HOME/samtools mpileup \
	-C50 \
	-f $REF_HOME/$REFGENOME \
	$ABBREV"_R3_noDup.bam" \
	>$ABBREV"_mPileUp_0_1_18.vcf"
	


	
#=======================	
#[6] generate fasta consensus from pileup

perl $RAY_SOFTWARE/pileup2fasta_v1-4.pl \
	-i $ABBREV"_mPileUp_0_1_18.vcf" \
	-o $ABBREV"__masked.fa"	\
	-g $ABBREV"__masked.gff"	\
	-b 8	\
	-s		\
	-V		

$SAMTOOLS_HOME/samtools index \
	-b \
	$ABBREV"_R3_noDup.bam" 	

$BEDTOOLS_HOME/bedtools genomecov \
	-ibam $ABBREV"_R3_noDup.bam" \
	-g $ABBREV"__masked.fa" \
	| grep genome >$ABBREV"_genome_cov.txt" 


#echo $ABBREV_fasta_finished" |  mailx -s $ABBREV"_fasta_finished" kev.am.sullivan@gmail.com


	
#echo $ABBREV_assembly_finished" |  mailx -s $ABBREV"_assembly_finished" kev.am.sullivan@gmail.com

done
sleep 5

