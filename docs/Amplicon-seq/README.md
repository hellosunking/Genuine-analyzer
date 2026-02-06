## Validations using Amplicon-seq

We have validated a handful of off-target loci using Amplicon-seq. The results are provided in files
analyzed using [CRISPResso2](https://github.com/pinellolab/CRISPResso2).

Here are the codes (with key parameters) we used during the analysis:
```
sid=sample_ID
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

File lists:

[CHANGE-V2-OT5](CRISPResso_on_CHANGE-V2-OT5.html)

[CIRCLE-V2-OT1](CRISPResso_on_CIRCLE-V2-OT1.html)

[CIRCLE-V2-OT5](CRISPResso_on_CIRCLE-V2-OT5.html)

[DISCOVER-V2-ot1](CRISPResso_on_DISCOVER-V2-ot1.html)

[DISCOVER-V2-ot2](CRISPResso_on_DISCOVER-V2-ot2.html)

[DISCOVER-V2-ot3](CRISPResso_on_DISCOVER-V2-ot3.html)

[DISCOVER-V2-ot4](CRISPResso_on_DISCOVER-V2-ot4.html)

[DISCOVER-V2-ot5](CRISPResso_on_DISCOVER-V2-ot5.html)

[Genuine-E-ot1](CRISPResso_on_Genuine-E-ot1.html)

[Genuine-E-ot2](CRISPResso_on_Genuine-E-ot2.html)

[Genuine-E-ot3](CRISPResso_on_Genuine-E-ot3.html)

[Genuine-E-ot4](CRISPResso_on_Genuine-E-ot4.html)

[Genuine-F-ot1](CRISPResso_on_Genuine-F-ot1.html)

[Genuine-F-ot2](CRISPResso_on_Genuine-F-ot2.html)

[Genuine-F-ot3](CRISPResso_on_Genuine-F-ot3.html)

[Genuine-F-ot4](CRISPResso_on_Genuine-F-ot4.html)

[Genuine-F-ot5](CRISPResso_on_Genuine-F-ot5.html)

[Genuine-H4-ot1](CRISPResso_on_Genuine-H4-ot1.html)

[Genuine-H4-ot2](CRISPResso_on_Genuine-H4-ot2.html)

[Genuine-H4-ot3](CRISPResso_on_Genuine-H4-ot3.html)

[Genuine-H4-ot4](CRISPResso_on_Genuine-H4-ot4.html)

[Genuine-H4-ot5](CRISPResso_on_Genuine-H4-ot5.html)

[Genuine-V1-ot1](CRISPResso_on_Genuine-V1-ot1.html)

[Genuine-V1-ot2](CRISPResso_on_Genuine-V1-ot2.html)

[Genuine-V1-ot3](CRISPResso_on_Genuine-V1-ot3.html)

[Genuine-V1-ot4](CRISPResso_on_Genuine-V1-ot4.html)

[Genuine-V1-ot5](CRISPResso_on_Genuine-V1-ot5.html)

[Genuine-V2-ot1](CRISPResso_on_Genuine-V2-ot1.html)

[Genuine-V2-ot2](CRISPResso_on_Genuine-V2-ot2.html)

[Genuine-V2-ot3](CRISPResso_on_Genuine-V2-ot3.html)

[Genuine-V2-ot4](CRISPResso_on_Genuine-V2-ot4.html)

[Genuine-V2-ot5](CRISPResso_on_Genuine-V2-ot5.html)

[Genuine-V3-ot1](CRISPResso_on_Genuine-V3-ot1.html)

[Genuine-V3-ot2](CRISPResso_on_Genuine-V3-ot2.html)

[Genuine-V3-ot3](CRISPResso_on_Genuine-V3-ot3.html)

[Genuine-V3-ot4](CRISPResso_on_Genuine-V3-ot4.html)

[Genuine-V3-ot5](CRISPResso_on_Genuine-V3-ot5.html)

[GUIDE-H4-ot1](CRISPResso_on_GUIDE-H4-ot1.html)

[GUIDE-H4-ot2](CRISPResso_on_GUIDE-H4-ot2.html)

[GUIDE-H4-ot3](CRISPResso_on_GUIDE-H4-ot3.html)

[GUIDE-V2-ot1](CRISPResso_on_GUIDE-V2-ot1.html)

[GUIDE-V3-ot1](CRISPResso_on_GUIDE-V3-ot1.html)

[GUIDE-V3-ot2](CRISPResso_on_GUIDE-V3-ot2.html)

[GUIDE-V3-ot3](CRISPResso_on_GUIDE-V3-ot3.html)

[GUIDE-V3-ot4](CRISPResso_on_GUIDE-V3-ot4.html)

[GUIDE-V3-ot5](CRISPResso_on_GUIDE-V3-ot5.html)

[PEAC-F-ot1](CRISPResso_on_PEAC-F-ot1.html)

[PEAC-V1-ot1](CRISPResso_on_PEAC-V1-ot1.html)

[PEAC-V1-ot2](CRISPResso_on_PEAC-V1-ot2.html)

[PEAC-V1-ot3](CRISPResso_on_PEAC-V1-ot3.html)

[PEAC-V1-ot4](CRISPResso_on_PEAC-V1-ot4.html)

[PEAC-V2-ot1](CRISPResso_on_PEAC-V2-ot1.html)

[PEAC-V2-ot2](CRISPResso_on_PEAC-V2-ot2.html)

[PEAC-V2-ot3](CRISPResso_on_PEAC-V2-ot3.html)

[PEAC-V2-ot4](CRISPResso_on_PEAC-V2-ot4.html)

[PEAC-V2-ot5](CRISPResso_on_PEAC-V2-ot5.html)

[PEAC-V3-ot1](CRISPResso_on_PEAC-V3-ot1.html)

[PEAC-V3-ot2](CRISPResso_on_PEAC-V3-ot2.html)

[PEAC-V3-ot3](CRISPResso_on_PEAC-V3-ot3.html)

[PEAC-V3-ot4](CRISPResso_on_PEAC-V3-ot4.html)

[PEAC-V3-ot5](CRISPResso_on_PEAC-V3-ot5.html)

[PEM-E-ot1](CRISPResso_on_PEM-E-ot1.html)

[PEM-E-ot2](CRISPResso_on_PEM-E-ot2.html)

[PEM-E-ot3](CRISPResso_on_PEM-E-ot3.html)

[PEM-E-ot4](CRISPResso_on_PEM-E-ot4.html)

[PEM-E-ot5](CRISPResso_on_PEM-E-ot5.html)

[V1-OT6](CRISPResso_on_V1-OT6.html)

[V1-OT7](CRISPResso_on_V1-OT7.html)

[V1-OT8](CRISPResso_on_V1-OT8.html)

[V1-OT9](CRISPResso_on_V1-OT9.html)

[V2-OT10](CRISPResso_on_V2-OT10.html)

[V2-OT11](CRISPResso_on_V2-OT11.html)

[V2-OT12](CRISPResso_on_V2-OT12.html)

[V2-OT13](CRISPResso_on_V2-OT13.html)

[V2-OT14](CRISPResso_on_V2-OT14.html)

[V2-OT15](CRISPResso_on_V2-OT15.html)

[V2-OT16](CRISPResso_on_V2-OT16.html)

[V2-OT17](CRISPResso_on_V2-OT17.html)

[V2-OT18](CRISPResso_on_V2-OT18.html)

[V2-OT19](CRISPResso_on_V2-OT19.html)

[V2-OT1](CRISPResso_on_V2-OT1.html)

[V2-OT20](CRISPResso_on_V2-OT20.html)

[V2-OT21](CRISPResso_on_V2-OT21.html)

[V2-OT22](CRISPResso_on_V2-OT22.html)

[V2-OT23](CRISPResso_on_V2-OT23.html)

[V2-OT24](CRISPResso_on_V2-OT24.html)

[V2-OT25](CRISPResso_on_V2-OT25.html)

[V2-OT26](CRISPResso_on_V2-OT26.html)

[V2-OT27](CRISPResso_on_V2-OT27.html)

[V2-OT28](CRISPResso_on_V2-OT28.html)

[V2-OT29](CRISPResso_on_V2-OT29.html)

[V2-OT2](CRISPResso_on_V2-OT2.html)

[V2-OT30](CRISPResso_on_V2-OT30.html)

[V2-OT31](CRISPResso_on_V2-OT31.html)

[V2-OT32](CRISPResso_on_V2-OT32.html)

[V2-OT33](CRISPResso_on_V2-OT33.html)

[V2-OT34](CRISPResso_on_V2-OT34.html)

[V2-OT35](CRISPResso_on_V2-OT35.html)

[V2-OT36](CRISPResso_on_V2-OT36.html)

[V2-OT37](CRISPResso_on_V2-OT37.html)

[V2-OT38](CRISPResso_on_V2-OT38.html)

[V2-OT39](CRISPResso_on_V2-OT39.html)

[V2-OT3](CRISPResso_on_V2-OT3.html)

[V2-OT4](CRISPResso_on_V2-OT4.html)

[V2-OT5](CRISPResso_on_V2-OT5.html)

[V2-OT6](CRISPResso_on_V2-OT6.html)

[V2-OT7](CRISPResso_on_V2-OT7.html)

[V2-OT8](CRISPResso_on_V2-OT8.html)

[V2-OT9](CRISPResso_on_V2-OT9.html)

[VEGFAs1-OT10](CRISPResso_on_VEGFAs1-OT10.html)

[VEGFAs1-OT11](CRISPResso_on_VEGFAs1-OT11.html)

[VEGFAs1-OT12](CRISPResso_on_VEGFAs1-OT12.html)

[VEGFAs1-OT13](CRISPResso_on_VEGFAs1-OT13.html)

[VEGFAs1-OT14](CRISPResso_on_VEGFAs1-OT14.html)

[VEGFAs1-OT15](CRISPResso_on_VEGFAs1-OT15.html)

[VEGFAs1-OT1](CRISPResso_on_VEGFAs1-OT1.html)

[VEGFAs1-OT2](CRISPResso_on_VEGFAs1-OT2.html)

[VEGFAs1-OT3](CRISPResso_on_VEGFAs1-OT3.html)

[VEGFAs1-OT4](CRISPResso_on_VEGFAs1-OT4.html)

[VEGFAs1-OT5](CRISPResso_on_VEGFAs1-OT5.html)

[VEGFAs1-OT6](CRISPResso_on_VEGFAs1-OT6.html)

[VEGFAs1-OT7](CRISPResso_on_VEGFAs1-OT7.html)

[VEGFAs1-OT8](CRISPResso_on_VEGFAs1-OT8.html)

[VEGFAs1-OT9](CRISPResso_on_VEGFAs1-OT9.html)

[VEGFA_site1_V1-OT1](CRISPResso_on_VEGFA_site1_V1-OT1.html)

[VEGFA_site1_V1-OT2](CRISPResso_on_VEGFA_site1_V1-OT2.html)

[VEGFA_site1_V1-OT3](CRISPResso_on_VEGFA_site1_V1-OT3.html)

[VEGFA_site1_V1-OT4](CRISPResso_on_VEGFA_site1_V1-OT4.html)

[VEGFA_site1_V1-OT5](CRISPResso_on_VEGFA_site1_V1-OT5.html)

[VEGFA_site1_V1-OT6](CRISPResso_on_VEGFA_site1_V1-OT6.html)

[VEGFA_site1_V1-OT7](CRISPResso_on_VEGFA_site1_V1-OT7.html)

[VEGFA_site1_V1-OT8](CRISPResso_on_VEGFA_site1_V1-OT8.html)

[VEGFA_site1_V1-OT9](CRISPResso_on_VEGFA_site1_V1-OT9.html)

