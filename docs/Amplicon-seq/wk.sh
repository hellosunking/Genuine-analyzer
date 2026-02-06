#!/bin/bash
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :
#
set -o nounset
set -o errexit
#command || { echo "command failed"; exit 1; }

for C in `ls *html`
do
	cat $C | perl -ne 's/hg38.HBV.EBV/hg38/; print' >zzz
	mv zzz $C
done 
