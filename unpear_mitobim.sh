#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N unpear-mitobim
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q Yoda
#$ -pe fill 1
#$ -P communitycluster

#This file will process raw Illumina data using Trimmomatic.  This will be followed by mapping to a reference mitochondrial genome using MIRA4 to create a new mitochondrial genome assembly.

cd /lustre/work/kevsulli/work/Sequences/MitoBim
mkdir UnPear 

BASEDIR=/lustre/work/kevsulli/work/Sequences/MitoBim/UnPear
WORKDIR=$BASEDIR/output

mkdir $WORKDIR
cd $WORKDIR

REFGENOME=Amont_Mito_genome.fa  	# your reference genome for the assembly
REF_HOME=/lustre/work/kevsulli/work/Sequences/MitoBim	#the location of your reference genome
REF=$(basename $REFGENOME _full.fa)

RAW_READS_HOME=/lustre/work/kevsulli/work/Sequences/MitoBim/Trim1   #the location of your raw data
mkdir $BASEDIR/data_raw
UNZIPPED_RAW_HOME=$BASEDIR/data_raw
mkdir $BASEDIR/support_files
SUPPORT_FILES=$BASEDIR/support_files	#where the support files like the adapter sequences will be located.

######
#set up alias' for major programs
######
BWA_HOME=/lustre/work/apps/bwa-0.6.2
SAMTOOLS_HOME=/lustre/work/apps/samtools-1.2
SAMTOOLS1_8_HOME=/lustre/work/apps/samtools-0.1.18
PICARD_HOME=/lustre/work/apps/picard-tools-1.91
BCFTOOLS_HOME=/lustre/work/apps/samtools-0.1.18/bcftools
RAY_SOFTWARE=/lustre/work/daray/software
TRIM_HOME=/lustre/work/apps/Trimmomatic-0.27
FASTX_HOME=/lustre/work/apps/fastx_toolkit-0.0.14/bin
VCFTOOLS_HOME=/lustre/work/daray/software/vcftools_0.1.12b/bin
MIRA_HOME=/lustre/work/apps/mira


for RAW_READ_FILE in $RAW_READS_HOME/*_2_filterPaired.fq
do
	SAMPLE_ID=$(basename $RAW_READ_FILE _2_filterPaired.fq)
#Unzip the raw reads into the processed_reads directory
#	gunzip -c $RAW_READS_HOME/$SAMPLE_ID"_L001_R1_001.fastq.gz" >$UNZIPPED_RAW_HOME/$SAMPLE_ID"_R1.fastq"
#	gunzip -c $RAW_READS_HOME/$SAMPLE_ID"_L001_R2_001.fastq.gz" >$UNZIPPED_RAW_HOME/$SAMPLE_ID"_R2.fastq"

mkdir $WORKDIR/$SAMPLE_ID
cd $WORKDIR/$SAMPLE_ID


#======================
#MIRA4 assembly 
#Create manifest.config for MIRA
echo -e "\n#manifest file for basic mapping assembly with illumina data using MIRA 4\n\nproject = initial-mapping-of-"$SAMPLE_ID"-to-UnPear-mt\n\njob=genome,mapping,accurate\n\nparameters = -NW:mrnl=0 -AS:nop=1 SOLEXA_SETTINGS -CO:msr=no\n\nreadgroup\nis_reference\ndata = $REF_HOME/$REFGENOME\nstrain = UnPear-mt-genome\n\nreadgroup = reads\ndefault_qual = 30\ndata = "$RAW_READS_HOME/$SAMPLE_ID"_1_filterPaired.fq" $RAW_READS_HOME/$SAMPLE_ID"_2_filterPaired.fq" $RAW_READS_HOME/$SAMPLE_ID"_1_filterUnPaired.fq" $RAW_READS_HOME/$SAMPLE_ID"_2_filterUnPaired.fq\ntechnology = solexa\nstrain = "$SAMPLE_ID"\n" > $SUPPORT_FILES/$SAMPLE_ID"_manifest.conf"

#Run MIRA
$MIRA_HOME/bin/mira $SUPPORT_FILES/$SAMPLE_ID"_manifest.conf"

#======================
#MITObim assembly
#Bait and iteratively map to the reference genome using MITObim
perl $RAY_SOFTWARE/MITObim_1.8.pl \
	-start 1 \
	-end 10 \
	-sample $SAMPLE_ID \
	-ref UnPear-mt-genome \
	-readpool $RAW_READS_HOME/$SAMPLE_ID"_1_filterPaired.fq" $RAW_READS_HOME/$SAMPLE_ID"_2_filterPaired.fq" $RAW_READS_HOME/$SAMPLE_ID"_1_filterUnPaired.fq" $RAW_READS_HOME/$SAMPLE_ID"_2_filterUnPaired.fq \
	-maf $WORKDIR/$SAMPLE_ID/"initial-mapping-of-"$SAMPLE_ID"-to-UnPear-mt_assembly"/"initial-mapping-of-"$SAMPLE_ID"-to-UnPear-mt_d_results"/"initial-mapping-of-"$SAMPLE_ID"-to-UnPear-mt_out.maf" \
	&> $SAMPLE_ID".log"

cd ..
done

#-bash-4.1$ /PATH/TO/MITObim.pl -start 1 -end 10 -sample testpool -ref Salpinus_mt_genome -readpool reads.fastq -maf initial-mapping-testpool-to-Salpinus-mt_assembly/initial-mapping-testpool-to-Salpinus-mt_d_results/initial-mapping-testpool-to-Salpinus-mt_out.maf &> log







