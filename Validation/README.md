## Validations using Amplicon-seq

We have validated a handful of off-target loci using Amplicon-seq. The results are provided in html files
analyzed using [CRISPResso2](https://github.com/pinellolab/CRISPResso2).

Here is the example codes (with key parameters) we used during the analysis:
```
R1=/path/to/Read1.fq.gz
R2=/path/to/Read2.fq.gz
amplicon=Sequence_of_the_amplicon
guideRNA=Sequence_of_the_guideRNA
bowtie2index=/path/to/bowtie2.index/hg38

CRISPResso -x $bowtie2index \
	--allele_plot_pcts_only_for_assigned_reference --write_detailed_allele_table \
	--plot_window_size 40 --max_rows_alleles_around_cut_to_plot 100 \
	--min_identity_score 50 --min_frequency_alleles_around_cut_to_plot 0.00001 \
	-r1 $R1 -r2 $R2 -a $amplicon -g $guideRNA -n $sid
```

