#!/bin/bash

## example codes for building index
[ -s hg38.fa.gz ] || wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
perl process.genome.pl hg38.fa.gz CRISPR.seq.fa >hg38.with.CRISPR.fa
../bwa index -p hg38 hg38.with.CRISPR.fa
ln -s hg38.with.CRISPR.fa hg38.fa

## optional for spRY
#perl process.genome.pl hg38.fa.gz spRY.fa >hg38spRY.fa
#../bwa index -p hg38spRY hg38spRY.fa

## optional for mm10
#wget https://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.fa.gz
#perl process.genome.pl mm10.fa.gz CRISPR.seq.fa >mm10.with.CRISPR.fa
#../bwa index -p mm10 mm10.with.CRISPR.fa
#ln -s mm10.with.CRISPR.fa mm10.fa

