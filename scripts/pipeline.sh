#!/usr/bin/env bash
#$ -S /bin/bash
#$ -M clara@matis.is
#$ -V
#$ -m bae
#$ -cwd

# ----- PROJECT -----

HOME_DIR=/project/microme/clara

project=osd
sequencer=miseq
run=ERR867776
libprep=""

PROJECT_DIR=$HOME_DIR/"$project"
cd $PROJECT_DIR

# ----- RAW DATA -----

# raw reads
RAW_READS_DIR=$PROJECT_DIR/data/"$sequencer"/"$run"

# ----- SOFTWARES AND TOOLS -----

# trimmomatic
TRIM_DIR=/project/microme/clara/tools/Trimmomatic-0.38
BWA_DIR=/project/microme/clara/tools/bwa-0.7.17

# kraken2 directory
KRAKEN2_DIR=/project/microme/clara/softwares/kraken2
KRAKEN2_DB_DIR=/project/microme/clara/databases/kraken2
# pick: silva  standard  custom	greengenes rdp
KRAKEN2_DB=custom

# krona directory
KRONA_DIR=/usr/local/bioinfo/KronaTools-2.6

# quast
QUAST_DIR=/project/microme/clara/tools/quast   

# ----- OUTPUT DIRECTORIES -----

# results directory
TRIM_RES_DIR="$PROJECT_DIR"/res/trimmomatic/"$sequencer"/"$run"
HUMAN_DIR=/project/microme/clara/databases/human/"$sequencer"/"$run"
MAPHUMAN="$PROJECT_DIR"/res/maphuman/"$sequencer"/"$run"
NOTHUMAN=$PROJECT_DIR/res/nothuman/"$sequencer"/"$run"
KRAKEN_RES_DIR=$PROJECT_DIR/res/kraken2/"$sequencer"/"$run"
KRONA_RES_DIR=$PROJECT_DIR/res/krona/"$sequencer"/"$run"
FLASH_RES_DIR=$PROJECT_DIR/res/flashmerged/"$sequencer"/"$run"
MEGAHIT_RES_DIR=$PROJECT_DIR/res/megahit/"$sequencer"/"$run"
PRODIGAL_RES_DIR=$PROJECT_DIR/res/prodigal/"$sequencer"/"$run"
QUAST_RES_DIR=$PROJECT_DIR/res/quast/"$sequencer"/"$run"

mkdir -p $TRIM_RES_DIR
mkdir -p $HUMAN_DIR
mkdir -p $MAPHUMAN
mkdir -p $NOTHUMAN
mkdir -p $KRAKEN_RES_DIR
mkdir -p $KRONA_RES_DIR
mkdir -p $FLASH_RES_DIR
mkdir -p $MEGAHIT_RES_DIR
mkdir -p $PRODIGAL_RES_DIR
mkdir -p $QUAST_RES_DIR

for file in $RAW_READS_DIR/*R1*fastq.gz
do

	#file=$(ls $RAW_READS_DIR/*R1*fastq.gz | head -2 | tail -1 )
	echo $file
	base=${file##*/}
	echo $base
	forward=${base%.*.*}
	echo $forward
	reverse=$(echo $forward | sed 's/R1/R2/g')
	echo $reverse
	smp=$(echo $forward | sed 's/_R1//g' | cut -f 5 -d ".")
	echo $smp

	# ----- RUN TRIMMOMATIC -----

	# java -jar $TRIM_DIR/trimmomatic-0.38.jar \
	# PE $RAW_READS_DIR/"$forward".fastq.gz $RAW_READS_DIR/"$reverse".fastq.gz \
	# $TRIM_RES_DIR/"$forward"_paired_trim.fastq.gz \
	# $TRIM_RES_DIR/"$forward"_unpaired_trim.fastq.gz \
	# $TRIM_RES_DIR/"$reverse"_paired_trim.fastq.gz \
	# $TRIM_RES_DIR/"$reverse"_unpaired_trim.fastq.gz \
	# ILLUMINACLIP:$TRIM_DIR/adapters/TruSeq3-PE-2.fa:2:30:10 \
	# LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50

	# ----- AVG READ LENGTH (after trimming) ------

	# https://bioinformatics.stackexchange.com/questions/4911/calculating-read-average-length-in-a-fastq-file-with-bioawk-awk
	#awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' $TRIM_RES_DIR/$prefix"1_paired_trim".fastq.gz

	# # ----- MERGE WITH FLASH -----

	# FOR MG
	# flash $TRIM_RES_DIR/"$forward"_paired_trim.fastq.gz \
	# $TRIM_RES_DIR/"$reverse"_paired_trim.fastq.gz \
	# -m 20 -r 290 -f 300 -s 500 --output-prefix="$smp" \
	# --output-directory=$FLASH_RES_DIR 2>&1 | tee $FLASH_RES_DIR/flash.log

	# FOR AMPLICON
	# flash $TRIM_RES_DIR/"$forward"_paired_trim.fastq.gz \
	# $TRIM_RES_DIR/"$reverse"_paired_trim.fastq.gz \
	# -m 15 -r 270 -f 390 --output-prefix="$smp" \
	# --output-directory=$FLASH_RES_DIR 2>&1 | tee $FLASH_RES_DIR/flash.log
	# # Percent combined: 46.43%

	# ----- KRAKEN2 -----

	$KRAKEN2_DIR/kraken2 \
	--db $KRAKEN2_DB_DIR/$KRAKEN2_DB \
	--threads 4 \
	--output $KRAKEN_RES_DIR/"$smp"_"$KRAKEN2_DB".kraken \
	$FLASH_RES_DIR/"$smp".extendedFrags.fastq \
	--report $KRAKEN_RES_DIR/"$smp"_"$KRAKEN2_DB".kraken.report \
	--unclassified-out $KRAKEN_RES_DIR/"$smp"_"$KRAKEN2_DB"_unclassified.kraken \
	--classified-out $KRAKEN_RES_DIR/"$smp"_"$KRAKEN2_DB"_classified.kraken

	# conver kraken output to json (biom) for later importation in R
	# kraken-biom $KRAKEN_RES_DIR/"$smp"_"$KRAKEN2_DB".kraken.report \
	# --max C --min G --fmt json -o $KRAKEN_RES_DIR/"$smp"_"$KRAKEN2_DB".biom.json
	
	# scp -r clara@130.208.252.141:/project/microme/clara/skafta/res/kraken2/*.biom.json /Users/Clara/Projects/skafta/res/kraken2/
	
	# ----- KRONA -----

	# cat $KRAKEN_RES_DIR/"$smp"_"$KRAKEN2_DB".kraken | cut -f 2,3 > $KRONA_RES_DIR/"$smp"_"$KRAKEN2_DB".krona
	# $KRONA_DIR/scripts/ImportTaxonomy.pl $KRONA_RES_DIR/"$smp"_"$KRAKEN2_DB".krona -o $KRONA_RES_DIR/"$smp"_"$KRAKEN2_DB".krona.html

	# # ---- MEGAHIT -----
	# megahit -r $FLASH_RES_DIR/"$smp".extendedFrags.fastq \
	# -m 0.5 -t 12  -o  $MEGAHIT_RES_DIR/"$smp"

	# ----- ASSEMBLY ASSESSMENT -----

	# check contigs
	# perl /project/microme/clara/mime/scripts/pl/contig-stats.pl < $MEGAHIT_RES_DIR/"$smp"/final.contigs.fa

	# with quast
	# $QUAST_DIR/metaquast.py -o $QUAST_DIR $MEGAHIT_RES_DIR/"$smp"/final.contigs.fa

	# # ---- PRODIGAL -----
	# prodigal -i $MEGAHIT_RES_DIR/"$smp"/final.contigs.fa \
	# -o $PRODIGAL_RES_DIR/"$smp".genes \
	# -a $PRODIGAL_RES_DIR/"$smp".faa -p meta

done
