#!/bin/bash
# for GRO-seq and CTCF motif analysis
# prepare hg19.refGene.tss.uniq.srt.bed
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
awk -F"\t" -v OFS="\t" '{if($4=="-"){s=$6-1;e=$6;}else{s=$5-1;e=$5}print $3,s,e,$13"|"$2,$12,$4}' refGene.txt | \
	sed 's/^chr//' | bedtools sort -i - > hg19.refGene.tss.srt.bed
perl unique_tss_bed.pl hg19.refGene.tss.srt.bed | \
	bedtools sort -i - > hg19.refGene.tss.uniq.srt.bed

################### GRO-seq ###################
# bowtie2; HOMER;

################### ChIA-PET ###################
# GSE39495 
# bedtools window -w 5000 -a ChIA-PET.tag -b GROseq

################### CTCF motif ###################
# random regions by bedtools
# bedtools random -l 1 -n 1000000 -g /work/04935/shaojf/stampede2/refs/hg19.clean.genome | awk -vOFS="\t" '{print $1,$2,$3,"random1M"}' | bedtools sort -i - > random_1e6.srt.bed
# Scan genome with fimo and prepare CTCF motif coordinate files
# CTCF motif: MA0139.1 in jaspar
# URL http://jaspar.genereg.net?ID=MA0139.1&rm=present&collection=CORE
fimo --o fimo_ctcf --max-stored-scores 100000000 \
	CTCF.meme Homo_sapiens.GRCh37.dna.primary_assembly.fa
grep -v "^#" fimo_ctcf/fimo.gff | \
	awk -F"\t" -v OFS="\t" '{print $1,$4,$5,$9,$6,$7}' | \
	bedtools sort -i - > hg19.fimo.MA0139.1.srt.bed

# Pick up centers for each CTCF sites.
sed 's/qvalue= /qvalue=/' hg19.fimo.MA0139.1.srt.bed | \
	awk -vOFS="\t" '{print $1,int(($2+$3)/2-1),int(($2+$3)/2),$4,$5,$6}' > hg19.fimo.MA0139.1.mid.srt.bed
bedtools closest -D a -a hg19.refGene.tss.uniq.srt.bed -b hg19.fimo.MA0139.1.mid.srt.bed | \
	awk '$(NF-1)!="."' > hg19.refGene.tss.uniq.srt.bed.MA0139.1.mid

