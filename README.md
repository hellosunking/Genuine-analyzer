## Pipeline for analyzing GENUINE-seq data (Han et al.)

## 1. Installation
Use the following code to download all the files in this repo
```
git clone https://github.com/hellosunking/Genuine-analyzer.git
```
Now you will get a directory named `Genuine-analyzer`, the main script is `Genuine.analyzer.sh` under this directory.
Note that this package also contains several annotation files that are needed during analysis.

In addition, you need to install the [circos](https://circos.ca/ "circos") software for data visualization.

## 2. Build genome index
Genome index must be built before actual analysis (just do it once). You need to prepare the genome sequence in FASTA format,
then combine it with CRISPR sequence and build an index for BWA aligner. If you already had the FASTA file (e.g. for hg38),
you can link it to the `genome` directory (gzipped file is OK):
```
ln -s /path/to/hg38.fa.gz genome/
```
Or you can download it from UCSC Genome Browser:
```
wget --no-check-certificate https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz -O genome/hg38.fa.gz
```
Now change to `genome` directory, and use the `genome/process.genome.pl` script to process the FASTA files and build the index:
```
cd genome
perl process.genome.pl hg38.fa.gz CRISPR.seq.fa >hg38.with.CRISPR.fa
../bwa index -p hg38 hg38.with.CRISPR.fa
ln -s hg38.with.CRISPR.fa hg38.fa
```
Note that `CRISPR.seq.fa` is provided within this package. Now you can leave the `genome` directory and use `hg38` as reference genome.

## 3. Perform analysis
Run the `Genuine.analyzer.sh` without parameters to see the usage:
```
user@linux$ sh Genuine.analyzer.sh 
Usage: Genuine.analyzer.sh <R1.fq[.gz]> <R2.fq[.gz]> <sid> <design.seq> <genome> [enzyme=spCas9|spRY] [noise='genome'.bgNoise] [thread=16]
```
The first 2 paramaters are the raw sequencing files (in FASTQ format) generated using GENUINE-seq protocol, the 3rd parameter is the
sample id (a directory with this name would be created and all outputs would be written into this directory), the 4th parameter is the
designed guide sequence (with "N" allowed), and the 5th paramter is the genome id (e.g., hg38). The remaining paramters are optional:
6th is the enzyme, could be spCas9 or spRY (Note that they require different index files), 7th is the background noise files (we provided
the files for human and mouse), and 8th is the running thread (default is 16).

We had prepared a walkthrough script that would automatically download a testing dataset, build genome index, and run Genuine-analyzer
on this dataset. Make sure that you have `curl` installed (which is installed by default in most Linux distributions) and internet access,
then change to `walkthrough` directory and run it:
```
user@linux$ cd walkthrough; sh walkthrough.sh
```

You should expect to see outputs very similar to the following (note that this script may takes time to finish):
```
=================================================================================================
This script is designed to configure and run Genuine-Analyzer on a testing dataset automatically.
=================================================================================================

Checking testing datasets ...
Downloading HEK293T-spCas9-V3-r2.read1.fq.gz ...
Downloading HEK293T-spCas9-V3-r2.read2.fq.gz ...
Checking files ...
HEK293T-spCas9-V3-r2.read1.fq.gz: OK
HEK293T-spCas9-V3-r2.read2.fq.gz: OK

Checking genome index ...
INFO: download genome fasta from UCSC ...
INFO: building genome index (be patient) ...
[bwa_index] Pack FASTA... 18.40 sec
[bwa_index] Construct BWT for the packed sequence...
[BWTIncCreate] textLength=6176591378, availableWord=446607504
[BWTIncConstructFromPacked] 10 iterations done. 99999986 characters processed.
... [lots of outputs are skipped here]
[BWTIncConstructFromPacked] 680 iterations done. 6165967586 characters processed.
[bwt_gen] Finished constructing BWT in 686 iterations.
[bwa_index] 1974.43 seconds elapse.
[bwa_index] Update BWT... 15.58 sec
[bwa_index] Pack forward-only FASTA... 11.14 sec
[bwa_index] Construct SA from BWT and Occ... 799.81 sec
[main] Version: 0.7.17-r1188
[main] CMD: ../bwa index -p hg38 hg38.with.CRISPR.fa
[main] Real time: 2900.643 sec; CPU: 2819.367 sec

Runing Genome-Analyzer (with 16 threads) ...
The command line is:
sh ../Genuine.analyzer.sh HEK293T-spCas9-V3-r2.read1.fq.gz HEK293T-spCas9-V3-r2.read2.fq.gz HEK293T-V3-r2 GGTGAGTGAGTGTGTGCGTGNGG hg38

## Parameters
Reads : HEK293T-spCas9-V3-r2.read1.fq.gz / HEK293T-spCas9-V3-r2.read2.fq.gz
Genome: hg38
Enzyme: spCas9
Index : /mnt/github/Genuine-analyzer/genome/hg38
Noise : /mnt/github/Genuine-analyzer/bgNoise/hg38.bgNoise
# Step 1: Read preprocess
# Step 2: Alignment
[main] Version: 0.7.17-r1188
[main] CMD: /mnt/github/Genuine-analyzer/bwa mem -v 2 -t 32 /mnt/github/Genuine-analyzer/genome/hg38 unc.border.fq
[main] Real time: 48.060 sec; CPU: 1081.455 sec
[main] Version: 0.7.17-r1188
[main] CMD: /mnt/github/Genuine-analyzer/bwa mem -5SP -v 2 -t 32 /mnt/github/Genuine-analyzer/genome/hg38 unc.pass.R1.fq unc.pass.R2.fq
[main] Real time: 25.421 sec; CPU: 447.230 sec
[main] Version: 0.7.17-r1188
[main] CMD: /mnt/github/Genuine-analyzer/bwa mem -v 2 -t 32 /mnt/github/Genuine-analyzer/genome/hg38 ext.border.fq
[main] Real time: 297.799 sec; CPU: 9033.967 sec
[main] Version: 0.7.17-r1188
[main] CMD: /mnt/github/Genuine-analyzer/bwa mem -5SP -v 2 -t 32 /mnt/github/Genuine-analyzer/genome/hg38 ext.pass.R1.fq ext.pass.R2.fq
[main] Real time: 70.823 sec; CPU: 1249.871 sec
# Step 3: Extract candidates
# Step 4: Filter candidates
Loading bgNoise: /mnt/github/Genuine-analyzer/bgNoise/hg38.bgNoise
=> Valid loci: 127096, supporting reads: 39848021.
Loading genome: hg38 (/mnt/github/Genuine-analyzer/genome/hg38.fa)
Loading HEK293T-V3-r2.cutpoints
=> Valid loci: 83648, supporting reads: 13259999.
Filtering cutpoints against GGTGAGTGAGTGTGTGCGTGNGG
Loading HEK293T-V3-r2.circles
Loading HEK293T-V3-r2.translocation
## Step 5: Data visualization
## Step 6: Clean up
#### Analysis finished. ####

Done. The outputs are stored in 'HEK293T-V3-r2' directory, and
the major results should be identical to those in 'example.output/' directory.
Please refer to README.md file for interpreting the results.
```

## 4. Output files
The outputs are all written in the directory with the name of the sample id paramater.
The following files are the main results (where `XXX` is the sample id in the 3r parameter):

- `XXX.final.stat`: key statistics during the analysis
- `XXX.targets.html`: target sites in colored HTML format
- `XXX.targets.loci.info`: target sites in plain text
- `XXX.structural.variants`: structural variants in plain text (empty for none-detected)
- `XXX.links.png` and `XXX.links.svg`: visualization for structural variants
- `XXX.crispr.insertions.png` and `XXX.crispr.insertions.svg`: visualization for CRISPR insertions

Other files are intermediate, you can ignore them after a successful analysis.

