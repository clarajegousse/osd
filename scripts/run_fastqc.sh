# ----- SET DIRECTORIES and FILENAMES -----

# raw reads
RAW_READS_DIR=/project/microme/clara/osd/data/miseq/ERR771024/

# working filenames
f_reads=ERR771024_1
r_reads=ERR771024_2

# trimmomatic and FastQC
TRIM_DIR=/project/microme/clara/tools/Trimmomatic-0.38
FASTQC_DIR=/usr/local/bioinfo/FastQC_v0.11.7

# results directory
TRIM_RES_DIR=/project/microme/clara/osd/res/trimmomatic/miseq/ERR771024
FASTQC_RES_DIR=/project/microme/clara/osd/res/fastqc/miseq/ERR771024

mkdir -p $TRIM_RES_DIR
mkdir -p $FASTQC_RES_DIR

# ----- RE-RUN FastQC -----

gunzip $RAW_READS_DIR/$f_reads.fastq.gz $RAW_READS_DIR/$r_reads.fastq.gz

# FastQC does not accept compressed files!
# Only with fastq files not fastq.gz
$FASTQC_DIR/fastqc -o $FASTQC_RED_DIR -f fastq $RAW_READS_DIR/$f_reads.fastq.gz

$FASTQC_DIR/fastqc -o $FASTQC_RES_DIR -f fastq $RAW_READS_DIR/$r_reads.fastq.gz

# download the outputfile from FastQC
# scp -r clara@130.208.252.141:/project/microme/clara/mime/res/fastqc/hiseq/run4789/stdshotgun/* /Users/Clara/Projects/mime/res/fastqc/hiseq/run4789/stdshotgun
