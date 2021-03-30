#!/bin/bash
####
# download ICGC data from ICGC data portal
# combined the WGS mutation 
for f in ICGC.release28/simple_somatic_mutation.open.*.tsv.gz
do
    gunzip -c $f | awk -F "\t" '$34=="WGS" || NR==1' | cut -f 1-4,9-17,26,29-30 >> ICGC.release28.simple_somatic_mutation.open.WGS.tsv
done

awk -F"\t" -v OFS="\t" '$1!="icgc_mutation_id"{print $5,$6,$7,$3"|"$4"|"$1"|"$10"|"$12"|"$13}' \
    ICGC.release28.simple_somatic_mutation.open.WGS.tsv | uniq | \
    sed 's? ??g' > ICGC.release28.simple_somatic_mutation.open.WGS.tsv.bed
bedtools sort -i ICGC.release28.simple_somatic_mutation.open.WGS.tsv.bed | uniq > ICGC.release28.simple_somatic_mutation.open.WGS.tsv.srt.bed

# combined the WGS deletion 
for f in ICGC.release28/structural_somatic_mutation.*.tsv.gz
do
    gunzip -c $f | tail -n +2 | cut -f 1-4,7-8,10,12-13,17-18,22-23 >> ICGC.release28.structural_somatic_mutation.tsv
done

grep -i deletion ICGC.release28.structural_somatic_mutation.tsv | \
    cut -f 1-3,6,8-9,11 > ICGC.release28.structural_somatic_mutation.del.tsv
awk -v OFS="\t" '{a=$6;b=$7;if($6>$7){a=$7;b=$6;}print $5,a,b,$2"|"$1"|"$3"|"$4}' ICGC.release28.structural_somatic_mutation.del.tsv | \
    bedtools sort -i - > ICGC.release28.structural_somatic_mutation.del.bed

####
# get the neareast gene annotation for ICGC WGS mutation and deletion
bedtools closest -D b \
    -a <(awk '{print "chr"$0}' ICGC.release28.simple_somatic_mutation.open.WGS.tsv.srt.bed) \
    -b <(awk '{print "chr"$0}' hg19.refGene.tss.uniq.srt.bed) > ICGC.release28.WGS.refGene.tss.txt
 
bedtools closest -D b \
    -a <(awk '{print "chr"$0}' ICGC.release28.structural_somatic_mutation.del.bed) \
    -b <(awk '{print "chr"$0}' hg19.refGene.tss.uniq.srt.bed) > ICGC.release28.del.refGene.tss.txt

# to get del.win1k.promoter.tsv
bedtools window -w 1000 -a ICGC.release28.structural_somatic_mutation.del.bed \
    -b hg19.refGene.tss.uniq.srt.bed > ICGC.release28.del.refGene.promoter.txt  
Rscript prepare.del.win1k.promoter.tsv.R

####
# get the 200Kb neighbors for cosmic genes
tail -n +2 Census_allMon.Feb.10.17_05_43.2020.sim.tsv | \
    cut -f 1 | grep -wf /dev/stdin hg19.refGene.tss.uniq.srt.bed | grep -v -e "\\-AS" -e "\\-OT" -e "\\-IT1" -e "\\-IT3" > allCosmic.tss.bed
# bedtools window -w 200000 -a allCosmic.tss.bed -b hg19.refGene.tss.uniq.srt.bed | \
#     sed 's?|?\t?' | awk -v OFS="\t" '{print $8,$9,$10,$4":"$11,$12,$13}' > tmp
bedtools intersect -v -a <(bedtools window -w 200000 -a allCosmic.tss.bed -b hg19.refGene.tss.uniq.srt.bed | \
    sed 's?|?\t?' | awk -v OFS="\t" '{print $8,$9,$10,$4":"$11,$12,$13}') -b allCosmic.tss.bed > allCosmic.200k.neigbors.txt
cut -f 4 allCosmic.200k.neigbors.txt | cut -f 1 -d "|" | sed 's?:?\t?' > allCosmic.200k.neigbors.gene.txt

awk '{print $4"\t"$4}' hg19.refGene.tss.uniq.srt.bed | sed 's?|?\t?' | cut -f 1,3 > hg19.refGene.tss.gene.txt
bedtools window -w 200000 -a hg19.refGene.tss.uniq.srt.bed -b hg19.refGene.tss.uniq.srt.bed > refgene.200k.neighbors.txt
paste <(cut -f 4 refgene.200k.neighbors.txt | sed 's?|?\t?' | cut -f 1) \
    <(cut -f 10 refgene.200k.neighbors.txt | sed 's?|?\t?' | cut -f 1) | awk '$1!=$2' | grep -v "MIR" > refgene.200k.neighbors.sim.txt
cut -f 1 refgene.200k.neighbors.sim.txt | sort | uniq -c | awk '{print $2"\t"$1}' > refgene.200k.neighbors.sim.count.txt
####
# download ICGC data from ICGC data portal
# get expression and deletion
for f in ICGC.release28/exp_seq*.gz
do
    ff=`basename $f | sed 's?.gz??'`
    zcat $f | tail -n +2 | cut -f 1-3,13 | sort | uniq | awk '{print "exp_seq\t"$0}' > stats/$ff &
done

for f in ICGC.release28/exp_array*.gz
do
    ff=`basename $f | sed 's?.gz??'`
    zcat $f | tail -n +2 | cut -f 1-3,13 | sort | uniq | awk '{print "exp_array\t"$0}' > stats/$ff &
done


for f in ICGC.release28/structural_somatic_mutation*.gz
do
    ff=`basename $f | sed 's?.gz??'`
    zcat $f | tail -n +2 | cut -f 1-3,39 | sort | uniq | awk '{print "SV\t"$0}' > stats/$ff &
done

cd stats
cat simple_somatic_mutation.open.*.tsv > allstats.simple_somatic_mutation.open.tsv 
cat exp_seq.*.tsv > allstats.exp_seq.tsv 
cat exp_array.*.tsv > allstats.exp_array.tsv 
cat structural_somatic_mutation.*.tsv > allstats.structural_somatic_mutation.tsv 
cut -f 2-5 allstats.simple_somatic_mutation.open.tsv | sort | uniq -c | \
    awk -vOFS="\t" '{print $2,$3,$4,$5,$1}' > allstats.simple_somatic_mutation.open.count.tsv 
cd ..

echo "Cancer Donor Mutation" | sed 's? ?\t?g' > ICGC.release28.WGS.mut.donor.count
awk -vOFS="\t" '{print $3,$2,$1}' ICGC.release28.simple_somatic_mutation.open.WGS.tsv | grep -v "project_code" | sort -u >> ICGC.release28.WGS.mut.donor.count
tail -n +2 ICGC.release28.WGS.mut.donor.count | cut -f 1,3 | sort | uniq -c | \
    awk 'BEGIN{print "Cancer\tMutation\tDonorCount"}{print $2"\t"$3"\t"$1}' > ICGC.release28.WGS.mut.donor.count.1
tail -n +2 ICGC.release28.WGS.mut.donor.count | cut -f 1,2 | sort -u | cut -f 1 | sort | uniq -c | \
    awk 'BEGIN{print "Cancer\tDonorCount"}{print $2"\t"$1}' > ICGC.release28.WGS.mut.donor.count.2
cut -f 4 ICGC.release28.structural_somatic_mutation.del.bed |sed 's?|?\t?g' |  cut -f 1-2 | \
    sort -u | cut -f 1 | sort | uniq -c | awk '{print $2"\t"$1}' > ICGC.release28.WGS.del.donor.count.2

####
# prepare counts for mutations statistics
echo "Cancer Specimen Mutation" | sed 's? ?\t?g' > ICGC.release28.WGS.mut.count
sed 's?|?\t?;s?|?\t?;s?|?\t?;' ICGC.release28.WGS.refGene.tss.txt | \
    awk -F"\t" -vOFS="\t" '{print $4,$5,$6}' | sort -u >> ICGC.release28.WGS.mut.count
    
echo "Cancer Mutation Gene Distance Count" | sed 's? ?\t?g' > ICGC.release28.WGS.refGene.tss.promoter.1k.count
awk '$6!=-1 && $11>-1000 && $11<1000' ICGC.release28.WGS.refGene.tss.txt | sed 's?|?\t?;s?|?\t?;s?|?\t?;' | \
    awk -F"\t" -vOFS="\t" '{print $4,$6,$11,$14}' | sort | uniq -c | awk -v OFS="\t" '{print $2,$3,$4,$5,$1}'| \
    sort -k 5,5nr >> ICGC.release28.WGS.refGene.tss.promoter.1k.count

# download specimen.all_projects.tsv from ICGC data portal
echo "Cancer Specimen Donor Gene Mutation Distance2TSS TSS" | sed 's? ?\t?g' > ICGC.release28.WGS.refGene.tss.promoter.1k.mut.gene.dis.tss
awk -F"\t" -v OFS="\t" '$6!=-1 && $11>-1000 && $11<1000{print $4,$8,$11,$8}' ICGC.release28.WGS.refGene.tss.txt | \
    sed 's?|?\t?;s?|?\t?;s?|?\t?' | cut -f 1-3,5,6,7 | sed 's?|?\t?' | cut -f 1-3,4,6,7 | \
    perl add_any_2files_together <(cut -f 1,5 specimen.all_projects.tsv) /dev/stdin 0 1 | \
    awk -v OFS="\t" '{print $1,$2,$8,$4,$3,$5,$6}' | sort -u >> ICGC.release28.WGS.refGene.tss.promoter.1k.mut.gene.dis.tss

####
# count promoter mutations
echo "Mutation Gene Distance Count" | sed 's? ?\t?g' > ICGC.release28.WGS.refGene.tss.promoter.1k.count.1
awk '$6!=-1 && $11>-1000 && $11<1000' ICGC.release28.WGS.refGene.tss.txt | sed 's?|?\t?;s?|?\t?;s?|?\t?;' | \
    awk -F"\t" -vOFS="\t" '{print $6,$11,$14}' | sort | uniq -c | awk -v OFS="\t" '{print $2,$3,$4,$5,$1}'| \
    sort -k 4,4nr >> ICGC.release28.WGS.refGene.tss.promoter.1k.count.1
    
echo "Cancer Mutation Gene Distance Count" | sed 's? ?\t?g' > ICGC.release28.WGS.refGene.tss.promoter.10k.count
awk '$6!=-1 && $11>-10000 && $11<1000' ICGC.release28.WGS.refGene.tss.txt | sed 's?|?\t?;s?|?\t?;s?|?\t?;' | \
    awk -F"\t" -vOFS="\t" '{print $4,$6,$11,$14}' | sort | uniq -c | awk -v OFS="\t" '{print $2,$3,$4,$5,$1}'| \
    sort -k 5,5nr >> ICGC.release28.WGS.refGene.tss.promoter.10k.count

echo "Mutation Gene Distance Count" | sed 's? ?\t?g' > ICGC.release28.WGS.refGene.tss.promoter.10k.count.1
awk '$6!=-1 && $11>-10000 && $11<1000' ICGC.release28.WGS.refGene.tss.txt | sed 's?|?\t?;s?|?\t?;s?|?\t?;' | \
    awk -F"\t" -vOFS="\t" '{print $6,$11,$14}' | sort | uniq -c | awk -v OFS="\t" '{print $2,$3,$4,$5,$1}'| \
    sort -k 4,4nr >> ICGC.release28.WGS.refGene.tss.promoter.10k.count.1

head -1 ICGC.release28.WGS.refGene.tss.promoter.10k.count.1 > ICGC.release28.WGS.refGene.tss.promoter.10k.count.11
tail -n +2 ICGC.release28.WGS.refGene.tss.promoter.10k.count.1 | sed 's?|?\t?' | \
    awk -vOFS="\t" '{print $1,$2,$4,$5}' >> ICGC.release28.WGS.refGene.tss.promoter.10k.count.11


####
# prepare input for extended.Fig.9a.plot.R
awk '$6!=-1' ICGC.release28.WGS.refGene.tss.txt | cut -f 4,11 | \
    sed 's?|?\t?;s?|?\t?;s?|?\t?;' | cut -f 1-3,5 > ICGC.release28.WGS.refGene.tss.dis.sim.txt
# Rscript extended.Fig.9a.plot.R
####
# make MAF file for oncoplot
awk '$6!=-1' ICGC.release28.WGS.refGene.tss.txt | sed 's?|?\t?;s?|?\t?;s?|?\t?;s?|?\t?;s?|?\t?;s?|?\t?;' | \
    awk -F"\t" -vOFS="\t" '{type="nonPromoter";if($17>0 && $17<100){type="TSS+100"}if($17>-500 && $17<0){type="TSS-500"}if($17>-5000 && $17<=-500){type="TSS-5K"}if($17>=100 && $17<1000){type="TSS+1K"}print $13,$1,$2,$3,$8,$9,$7,type,$4"|"$5}' \
    > ICGC.release28.WGS.refGene.tss.detail.maf
grep -v "nonPromoter" ICGC.release28.WGS.refGene.tss.detail.maf > ICGC.release28.WGS.refGene.tss.detail.promoter.maf

head -1 ICGC.release28.WGS.refGene.tss.detail.promoter.maf > ICGC.release28.WGS.refGene.tss.promoter.maf
awk '$6!=-1 && $11>-1000 && $11<1000' ICGC.release28.WGS.refGene.tss.txt | sed 's?|?\t?;s?|?\t?;s?|?\t?;s?|?\t?;s?|?\t?;s?|?\t?;' | \
    awk -F"\t" -vOFS="\t" '{type="TSS+(-1000,1000)";print $13,$1,$2,$3,$8,$9,$7,type,$4"|"$5}' \
    >> ICGC.release28.WGS.refGene.tss.promoter.maf
# Rscript extended.Fig.9efg.plot.R

####
# scan CTCF motif and find out CTCF-motif-disrupting mutation
# prepare mutation list file
# download simple_somatic_mutation.aggregated.vcf.gz from ICGC data portal
zcat simple_somatic_mutation.aggregated.vcf.gz | tail -n +15 | grep -v "^MT" | cut -f 1-5 > ICGC.release28.simple.mutation.txt
awk -vOFS="\t" '{print $1,$2-1,$2+length($4)-1,$3,$4,$5}' ICGC.release28.simple.mutation.txt > ICGC.release28.simple.mutation.bed
# run fimo
refseq=/data/shaojf/myReference/bwa-index/hg19.fa
CTCF=JASPAR_CORE_2018_vertebrates.CTCF.meme 
mylist=ICGC.release28.simple.mutation.bed
pre=ICGC.release28.simple.mutation
awk -F"\t" -vOFS="\t" '{a=$2-21;if($2-21<0){a=0;}print "chr"$1,a,$2}' $mylist | \
    bedtools getfasta -fi $refseq -bed - > $pre.left.fa
awk -F"\t" -vOFS="\t" '{print "chr"$1,$3,$3+20}' $mylist | \
    bedtools getfasta -fi $refseq -bed - > $pre.right.fa

paste $pre.left.fa <(awk -F"\t" '{print ">"$4"\n"$5}' $mylist) $pre.right.fa | sed 's/\t//g' > $pre.ref.fa
paste $pre.left.fa <(awk -F"\t" '{print ">"$4"\n"$6}' $mylist) $pre.right.fa | sed 's/\t//g' > $pre.mut.fa

paste <(sed -n '1~2p' $pre.ref.fa | awk -F">" '{print ">"$3}') <(sed -n '2~2p' $pre.ref.fa) | tr "\t" "\n" > a
paste <(sed -n '1~2p' $pre.mut.fa | awk -F">" '{print ">"$3}') <(sed -n '2~2p' $pre.mut.fa) | tr "\t" "\n" > b
mv a $pre.ref.fa 
mv b $pre.mut.fa


fimo --thresh 1e-3 -oc CTCF.$pre.ref $CTCF $pre.ref.fa &
fimo --thresh 1e-3 -oc CTCF.$pre.mut $CTCF $pre.mut.fa &
wait

grep -v "^#" CTCF.$pre.ref/fimo.tsv > $pre.ref.fimo.txt
grep -v "^#" CTCF.$pre.mut/fimo.tsv > $pre.mut.fimo.txt

####