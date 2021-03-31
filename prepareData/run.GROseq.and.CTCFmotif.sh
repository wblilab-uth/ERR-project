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
bedtools window -w 5000 -a ChIA-PET.tag -b GROseq

################### CTCF motif ###################
# random regions by bedtools
# bedtools random -l 1 -n 1000000 -g /work/04935/shaojf/stampede2/refs/hg19.clean.genome | awk -vOFS="\t" '{print $1,$2,$3,"random1M"}' | bedtools sort -i - > random_1e6.srt.bed
# Scan genome with fimo and prepare CTCF motif coordinate files
# CTCF motif: MA0139.1 in jaspar
# URL http://jaspar.genereg.net?ID=MA0139.1&rm=present&collection=CORE
fimo --o fimo_ctcf --max-stored-scores 100000000 CTCF.meme Homo_sapiens.GRCh37.dna.primary_assembly.fa
grep -v "^#" fimo_ctcf/fimo.gff | \
	awk -F"\t" -v OFS="\t" '{print $1,$4,$5,$9,$6,$7}' | \
	bedtools sort -i - > hg19.fimo.MA0139.1.srt.bed

# Pick up centers for each CTCF sites.
sed 's/qvalue= /qvalue=/' hg19.fimo.MA0139.1.srt.bed | awk -vOFS="\t" '{print $1,int(($2+$3)/2-1),int(($2+$3)/2),$4,$5,$6}' > hg19.fimo.MA0139.1.mid.srt.bed
bedtools closest -D a -a hg19.refGene.tss.uniq.srt.bed -b hg19.fimo.MA0139.1.mid.srt.bed | \
	awk '$(NF-1)!="."' > hg19.refGene.tss.uniq.srt.bed.MA0139.1.mid


# RAD21.5000.MCF-7.all.pairs.final
###////////////////////////###
##############################
# for motif of ChIA-PET paired EP
#####
awk -vOFS="\t" '{print $1,int(($2+$3)/2-1),int(($2+$3)/2),$4,$5,$6}' LWB.CTCF.bed | \
	bedtools sort -i - > LWB.CTCF.mid.srt.bed
motif=hg19.fimo.MA0139.1.mid.srt.bed
peaks=LWB.CTCF.mid.srt.bed
fac=CTCF
bedtools closest -D a -a $motif -b $peaks | awk '$(NF-6)!="."' > $motif.$peaks
awk '$NF> -100 && $NF <100' $motif.$peaks | cut -f 1-6 > $motif.$peaks.100
bedtools intersect -wao -a eRNA.all.srt.bed -b $motif.$peaks.100 | \
	awk '$(NF-6)!="."' > overlap.eRNA.all.$fac
bedtools intersect -wao -a gene.all.srt.bed -b $motif.$peaks.100 | \
	awk '$(NF-6)!="."' > overlap.gene.all.$fac
perl prepareData/add_any_2files_together.pl overlap.eRNA.all.$fac RAD21.5000.MCF-7.all.pairs 3 5 > RAD21.5000.MCF-7.all.pairs.$fac.eRNA
perl prepareData/add_any_2files_together.pl overlap.gene.all.$fac RAD21.5000.MCF-7.all.pairs.$fac.eRNA 3 9 > RAD21.5000.MCF-7.all.pairs.$fac.eRNA.gene
cut -f 1-10,18-19,22,31-32,35 RAD21.5000.MCF-7.all.pairs.$fac.eRNA.gene > RAD21.5000.MCF-7.all.pairs.$fac.eRNA.gene.sim
grep -v "/" RAD21.5000.MCF-7.all.pairs.$fac.eRNA.gene.sim > RAD21.5000.MCF-7.all.pairs.$fac.eRNA.gene.sim.paired
awk '{print $3":"$4":"$5":"$6":"$7":"$8":"$9":"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16}' RAD21.5000.MCF-7.all.pairs.$fac.eRNA.gene.sim | sort | uniq > RAD21.5000.MCF-7.all.pairs.$fac
#####
sed 's/chr//' hg19.Snyder.opt.IDR.MCF-7.FOXA1.bed | awk -vOFS="\t" '{print $1,int(($2+$3)/2-1),int(($2+$3)/2), $4,$5,$6}' | bedtools sort -i - > hg19.Snyder.opt.IDR.MCF-7.FOXA1.mid.srt.bed
motif=hg19.fimo.FOXA1.mid.srt.bed
peaks=hg19.Snyder.opt.IDR.MCF-7.FOXA1.mid.srt.bed
fac=FOXA1
bedtools closest -D a -a $motif -b $peaks | awk '$(NF-6)!="."' > $motif.$peaks
awk '$NF> -100 && $NF <100' $motif.$peaks | cut -f 1-6 > $motif.$peaks.100
bedtools intersect -wao -a eRNA.all.srt.bed -b $motif.$peaks.100 | awk '$(NF-6)!="."' > overlap.eRNA.all.$fac
bedtools intersect -wao -a gene.all.srt.bed -b $motif.$peaks.100 | awk '$(NF-6)!="."' > overlap.gene.all.$fac
perl ~/myScripts/add_any_2files_together.pl overlap.eRNA.all.$fac RAD21.5000.MCF-7.all.pairs 3 5 > RAD21.5000.MCF-7.all.pairs.$fac.eRNA
perl ~/myScripts/add_any_2files_together.pl overlap.gene.all.$fac RAD21.5000.MCF-7.all.pairs.$fac.eRNA 3 9 > RAD21.5000.MCF-7.all.pairs.$fac.eRNA.gene
cut -f 1-10,18-19,22,31-32,35 RAD21.5000.MCF-7.all.pairs.$fac.eRNA.gene > RAD21.5000.MCF-7.all.pairs.$fac.eRNA.gene.sim
grep -v "/" RAD21.5000.MCF-7.all.pairs.$fac.eRNA.gene.sim > RAD21.5000.MCF-7.all.pairs.$fac.eRNA.gene.sim.paired
awk '{print $3":"$4":"$5":"$6":"$7":"$8":"$9":"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16}' RAD21.5000.MCF-7.all.pairs.$fac.eRNA.gene.sim | sort | uniq > RAD21.5000.MCF-7.all.pairs.$fac

perl prepareData/add_any_2files_together.pl RAD21.5000.MCF-7.all.pairs.CTCF RAD21.5000.MCF-7.all.pairs.ESR1 0 0 | cut -f 1-7,9- > RAD21.5000.MCF-7.all.pairs.CTCF.ESR1
perl prepareData/add_any_2files_together.pl RAD21.5000.MCF-7.all.pairs.CTCF.ESR1 RAD21.5000.MCF-7.all.pairs.FOXA1 0 0 | cut -f 1-7,9- > RAD21.5000.MCF-7.all.pairs.CTCF.ESR1.FOXA1
sort RAD21.5000.MCF-7.all.pairs.CTCF.ESR1.FOXA1 | uniq | sed 's/:/\t/g' > RAD21.5000.MCF-7.all.pairs.final

# RAD21.5000.MCF-7.all.pairs.CTCF.FOXA1
######
perl prepareData/add_any_2files_together.pl RAD21.5000.MCF-7.all.pairs.CTCF RAD21.5000.MCF-7.all.pairs.FOXA1 0 0 | cut -f 1-7,9- | sort | uniq | sed 's/:/\t/g' > RAD21.5000.MCF-7.all.pairs.CTCF.FOXA1

