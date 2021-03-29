# ERR-project
Associated Scripts for ERR (Enhancer Release and Retargeting) project

# before starting:
1. get data associated for all the R scripts from http://www.fawnshao.info.sh/data/project.ERR.figdata.tar.gz
   extract it with "tar zxvf project.ERR.figdata.tar.gz"
   
2. or get all ICGC data from https://dcc.icgc.org/
   and run prepare.data.sh followed by prepare.del.win1k.promoter.tsv.R to generate related input files for R script
   
# File structure in project.ERR.figdata.tar.gz
├── annotations<br>
│   ├── gencode.v19.annotation.gene.probemap<br>
│   ├── hg19.refGene.tss.count<br>
│   └── hg19.refGene.tss.uniq.srt.bed<br>
├── Cosmic<br>
│   ├── allCosmic.200k.neigbors.gene.txt<br>
│   ├── allCosmic.200k.neigbors.tss.txt<br>
│   ├── allCosmic.tss.txt<br>
│   └── Census_allMon.Feb.10.17_05_43.2020.sim.tsv<br>
├── ICGC<br>
│   ├── allstats.simple_somatic_mutation.open.tsv<br>
│   ├── del.win1k.promoter.tsv<br>
│   ├── ICGC.cancers<br>
│   ├── ICGC.release28.simple.mutation.mut.fimo.txt<br>
│   ├── ICGC.release28.simple.mutation.ref.fimo.txt<br>
│   ├── ICGC.release28.WGS.del.donor.count.1<br>
│   ├── ICGC.release28.WGS.del.donor.count.2<br>
│   ├── ICGC.release28.WGS.mut.count.10k.combined.txt<br>
│   ├── ICGC.release28.WGS.mut.donor.count.2<br>
│   ├── ICGC.release28.WGS.refGene.tss.dis.sim.motif.txt<br>
│   ├── ICGC.release28.WGS.refGene.tss.dis.sim.txt<br>
│   ├── ICGC.release28.WGS.refGene.tss.promoter.10k.count.11<br>
│   ├── ICGC.release28.WGS.refGene.tss.promoter.1k.mut.gene.dis.tss<br>
│   ├── ICGC.release28.WGS.refGene.tss.promoter.maf<br>
│   ├── targets.all.exp_array.sim.tsv<br>
│   ├── targets.all.exp_seq.sim.tsv<br>
│   └── targets.all.exp.sim.tsv<br>
