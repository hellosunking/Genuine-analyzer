#!/bin/bash
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :
#
set -o nounset
set -o errexit
#command || { echo "command failed"; exit 1; }

if [ $# -lt 5 ]
then
	echo
	echo "Usage: $0 <R1.fq[.gz]> <R2.fq[.gz]> <sid> <design.seq> <genome> [enzyme=spCas9|spRY] [noise='genome'.bgNoise] [thread=16]"
	echo
	exit 2
fi >>/dev/stderr

## change log v231105
## circos plot for Crispr
## add raw sequence in *.targets.loci.info
## pipe bwa to cutpoints
## support both spCas9 and spRY
## v230929
## add mm10 support
## add Crispr sequence to both hg38p14 and mm10p6
## v230905
## use BWA for alignment
## discard no-Spike-in reads (most of them contains Spike-in, but have too many mismatch/indel)
## rescue the head/tail reads, serving as SUPPORTING reads for loop/cut point
## use Smith-Waterman algorithm to calculate the distance between designed sequence and genome
## the rmdup step is now optional
## bg noise is updated

wkSHELL=`readlink -f $0`
PRG=`dirname $wkSHELL`
sid=$3
designSeq=$4
genome=$5
enzyme=${6:-spCas9}
bgNoise=${7:-$PRG/bgNoise/$genome.bgNoise}
thread=${8:-32}

[ $enzyme == "Cas9" ] && enzyme="spCas9"
[ $enzyme == "RY"   ] && enzyme="spRY"

if [ "$enzyme" != "spCas9" ] && [ "$enzyme" != "spRY" ]
then
	echo "ERROR: only spCas9 and spRY are acceptable enzymes!"
	exit 1
fi >> /dev/stderr

seqKit=illumina
#bwaIndex=$PRG/genome/$genome$enzyme
bwaIndex=$PRG/genome/$genome

## check files
if [ ! -s $bwaIndex.sa ] || [ ! -s $bwaIndex.fa ]
then
	echo "ERROR: genome (index) files are missing!"
	exit 11
fi >>/dev/stderr

if [ ! -s $bgNoise ] && [ ! -L $bgNoise ]
then
	echo "ERROR: background noise file ($bgNoise) is missing!"
	exit 12
fi >>/dev/stderr

if [ -f $sid ]
then
	echo "WARNING: there is a file with the name of your output directory $sid!"
	exit 13
fi >>/dev/stderr

if [ -d $sid ]
then
	echo "WARNING: output directory $sid exists! Files will be over-written!"
fi >>/dev/stderr

## start analysis
echo "## Parameters"
echo Reads : $1 / $2
echo Genome: $genome
echo Enzyme: $enzyme
echo Index : $bwaIndex
echo Noise : $bgNoise

## preprocessing
echo "## Step 1: Read preprocess"
$PRG/ktrim -k $seqKit -t 8 -o $sid -1 $1 -2 $2 -c | \
$PRG/flash -m 10 -M 150 -t 8 -o $sid --interleaved-input - >$sid.flash.log
if [ ! -s $sid.trim.log ]
then
	echo "Error occured during data preprocessing. QUIT!!!"
	exit 1
fi >>/dev/stderr

mkdir -p $sid
mv $sid.*fastq $sid.*log $sid.hist* $sid
cd $sid
perl $PRG/process.ext.pl $sid.extendedFrags.fastq ext &
perl $PRG/process.unc.pl $sid.notCombined_1.fastq $sid.notCombined_2.fastq unc
wait

echo "## Step 2: Alignment"
$PRG/bwa mem -v 2 -t $thread $bwaIndex unc.border.fq | perl $PRG/extract.border.no.rmdup.pl - unc.border
$PRG/bwa mem -5SP -v 2 -t $thread $bwaIndex unc.pass.R1.fq unc.pass.R2.fq | perl $PRG/extract.pass.no.rmdup.pl - unc.pass

$PRG/bwa mem -v 2 -t $thread $bwaIndex ext.border.fq | perl $PRG/extract.border.no.rmdup.pl - ext.border
$PRG/bwa mem -5SP -v 2 -t $thread $bwaIndex ext.pass.R1.fq ext.pass.R2.fq | perl $PRG/extract.pass.no.rmdup.pl - ext.pass
## Note: exit here for bg noise modeling

echo "## Step 3: Extract candidates"
cat ext*.cutpoints unc*.cutpoints | perl $PRG/get.cutpoint.pl - $sid
gzip ext*cutpoints unc*cutpoints &

echo "## Step 4: Filter candidates"
perl $PRG/filter.v4.pl $bgNoise $designSeq $sid $genome
perl $PRG/make.html.v3.pl $designSeq $sid.targets.loci.info >$sid.targets.html
perl $PRG/make.stat.pl $sid >$sid.final.stat

## circos plots
echo "## Step 5: Data visualization"
perl $PRG/prepare.circos.v2.pl $sid $designSeq $genome
[ -s $sid.circos.conf ] && circos --conf $sid.circos.conf >>circos.log

perl $PRG/prepare.circos.CRISPR.pl $sid $genome $enzyme
[ -s $sid.circos.crispr.conf ] && circos --conf $sid.circos.crispr.conf >>circos.log

## clean up
echo "## Step 6: Clean up"
rm -f *fq *fastq *sam

cd ../
wait
echo "#### Analysis finished. ####"

