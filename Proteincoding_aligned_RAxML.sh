#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N PmexiRAxML
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q Chewie
#$ -pe sm 11
#$ -P communitycluster

/lustre/work/apps/standard-RAxML/raxmlHPC \
    -f a     -x 12368     -p 243531     -# 1000     -s ProteinCoding_aligned.phy     -m GTRGAMMA     -n ProteinCoding_aligned.out -o A_monten >run.log
