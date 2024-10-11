#!/bin/bash
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :
#
set -o nounset
#set -o errexit
#command || { echo "command failed"; exit 1; }

echo "================================================================================================="
echo "This script is designed to configure and run Genuine-Analyzer on a testing dataset automatically."
echo "================================================================================================="
echo

echo "Checking testing datasets ..."
if [ -s HEK293T-spCas9-V3-r2.read1.fq.gz ]
then
	echo "INFO: file HEK293T-spCas9-V3-r2.read1.fq.gz exists! I will use it."
else
	echo "Downloading HEK293T-spCas9-V3-r2.read1.fq.gz ..."
	curl --insecure "https://zenodo.org/records/13911526/files/HEK293T-spCas9-V3-r2.read1.fq.gz?download=1" -o HEK293T-spCas9-V3-r2.read1.fq.gz
	if [ $? != 0 ]
	then
		echo "FATAL ERROR: download file failed! Please check your internet connection." >/dev/stderr
		exit 1
	fi
fi

if [ -s HEK293T-spCas9-V3-r2.read2.fq.gz ]
then
	echo "INFO: file HEK293T-spCas9-V3-r2.read2.fq.gz exists! I will use it."
else
	echo "Downloading HEK293T-spCas9-V3-r2.read2.fq.gz ..."
	curl --insecure "https://zenodo.org/records/13911526/files/HEK293T-spCas9-V3-r2.read2.fq.gz?download=1" -o HEK293T-spCas9-V3-r2.read2.fq.gz
	if [ $? != 0 ]
	then
		echo "FATAL ERROR: download file failed! Please check your internet connection." >/dev/stderr
		exit 1
	fi
fi

echo "Checking files ..."
md5sum -c testing.dataset.md5
if [ $? != 0 ]
then
	echo "FATAL ERROR: the fastq files are incorrect, please delete them and run this script again." >/dev/stderr
	exit 1
fi

echo
echo "Checking genome index ..."
if [ -s ../genome/hg38.pac ]
then
	echo "INFO: genome index has been generated before."
else
	cd ../genome/
	if [ -s hg38.p14.fa.gz ]
	then
		echo "INFO: use the genome fasta file in this directory."
	else
		echo "INFO: download genome fasta from UCSC ..."
		curl --insecure "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/p14/hg38.p14.fa.gz" -o hg38.p14.fa.gz
		if [ $? != 0 ]
		then
			echo "FATAL ERROR: download file failed! Please check your internet connection." >/dev/stderr
			exit 1
		fi
	fi

	echo "INFO: building genome index (be patient) ..."
	perl process.genome.pl hg38.p14.fa.gz CRISPR.seq.fa >hg38.with.CRISPR.fa
	../bwa index -p hg38 hg38.with.CRISPR.fa
	ln -s hg38.with.CRISPR.fa hg38.fa

	cd ../walkthrough
fi

echo
echo "Runing Genome-Analyzer (with 16 threads) ..."
guide_sequence=GGTGAGTGAGTGTGTGCGTGNGG
echo "The command line is:"
echo 
echo "sh ../Genuine.analyzer.sh HEK293T-spCas9-V3-r2.read1.fq.gz HEK293T-spCas9-V3-r2.read2.fq.gz HEK293T-V3-r2 $guide_sequence hg38"
echo

sh ../Genuine.analyzer.sh HEK293T-spCas9-V3-r2.read1.fq.gz HEK293T-spCas9-V3-r2.read2.fq.gz HEK293T-V3-r2 $guide_sequence hg38

echo
echo "Done. The outputs are stored in 'HEK293T-V3-r2' directory, and"
echo "the major results should be identical to those in 'example.output/' directory."
echo "Please refer to README.md file for interpreting the results."

