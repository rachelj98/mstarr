##############################################
#Test for enriched pairwise overlap of null (baseline), IFNA, and dex regulatory windows using Fisher's exact tests (FET)
##############################################
#Create dataframes of regulatory windows for each treatment
setwd('~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation')
n<-read.table('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwol31250_FDR2.txt')
i<-read.table('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatIFN_FDR2.txt')
d<-read.table('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatdex_FDR2.txt')
cc=0.01
n<-subset(n,(unmeth.rna_beta>0 & fdr_RNAinunmeth<cc) | (meth.rna_beta>0 & fdr_RNAinmeth<cc))
i<-subset(i,(unmeth.rna_beta>0 & fdr_RNAinunmeth<cc) | (meth.rna_beta>0 & fdr_RNAinmeth<cc))
d<-subset(d,(unmeth.rna_beta>0 & fdr_RNAinunmeth<cc) | (meth.rna_beta>0 & fdr_RNAinmeth<cc))

#Create dataframes of windows with MD regulatory activity
nmd<-subset(n,fdr_int<cc)
imd<-subset(i,fdr_int<cc)
dmd<-subset(d,fdr_int<cc)

#Create dataframes of all windows tested for regulatory activity in each treatment
ntest<-read.table('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/windows1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwol31250.txt',sep='\t')
itest<-read.table(paste0('windows1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeat','IFN','.txt',sep=''),sep='\t')
dtest<-read.table(paste0('windows1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeat','dex','.txt',sep=''),sep='\t')
rownames(ntest)<-ntest$V4
rownames(itest)<-itest$V4
rownames(dtest)<-dtest$V4

#Flag whether each tested window was regulatory or not regulatory (and MD-reg or not MD-reg)
ntest$nreg<-ifelse(rownames(ntest) %in% rownames(n),1,0)
ntest$nmdreg<-ifelse(rownames(ntest) %in% rownames(nmd),1,0)
itest$ireg<-ifelse(rownames(itest) %in% rownames(i),1,0)
itest$imdreg<-ifelse(rownames(itest) %in% rownames(imd),1,0)
dtest$dreg<-ifelse(rownames(dtest) %in% rownames(d),1,0)
dtest$dmdreg<-ifelse(rownames(dtest) %in% rownames(dmd),1,0)

#Perform FET for null vs IFNA regulatory windows
#Reduce to windows tested in both null (baseline) and IFNA:
nitest<-merge(ntest,itest,by='row.names')
rownames(nitest)<-nitest$Row.names
head(nitest)

#Create contingency table and run FET 
nnotinot<-subset(nitest,nreg==0 & ireg==0)
nnotireg<-subset(nitest,nreg==0 & ireg==1)
nreginot<-subset(nitest,nreg==1 & ireg==0)
nregireg<-subset(nitest,nreg==1 & ireg==1)
aa<-nrow(nnotinot)
bb<-nrow(nnotireg)
cc<-nrow(nreginot)
dd<-nrow(nregireg)
values = matrix(c(aa,bb,cc,dd), nrow = 2);values
this<-fisher.test(values,alternative="two.sided");this
log2(this$estimate)

#Perform FET for null vs dex regulatory windows
#Reduce to windows tested in both null (baseline) and dex:
ndtest<-merge(ntest,dtest,by='row.names')
rownames(ndtest)<-ndtest$Row.names
head(ndtest)
nnotdnot<-subset(ndtest,nreg==0 & dreg==0)
nnotdreg<-subset(ndtest,nreg==0 & dreg==1)
nregdnot<-subset(ndtest,nreg==1 & dreg==0)
nregdreg<-subset(ndtest,nreg==1 & dreg==1)
aa<-nrow(nnotdnot)
bb<-nrow(nnotdreg)
cc<-nrow(nregdnot)
dd<-nrow(nregdreg)
values = matrix(c(aa,bb,cc,dd), nrow = 2);values
this<-fisher.test(values,alternative="two.sided");this
log2(this$estimate)

#Perform FET for IFNA vs dex regulatory windows
#Reduce to windows tested in both IFNA and dex:
idtest<-merge(itest,dtest,by='row.names')
rownames(idtest)<-idtest$Row.names
head(idtest)
inotdnot<-subset(idtest,ireg==0 & dreg==0)
inotdreg<-subset(idtest,ireg==0 & dreg==1)
iregdnot<-subset(idtest,ireg==1 & dreg==0)
iregdreg<-subset(idtest,ireg==1 & dreg==1)
aa<-nrow(inotdnot)
bb<-nrow(inotdreg)
cc<-nrow(iregdnot)
dd<-nrow(iregdreg)
values = matrix(c(aa,bb,cc,dd), nrow = 2);values
this<-fisher.test(values,alternative="two.sided");this
log2(this$estimate)

##############################################
#Perform TFBS enrichment analysis in HOMER on dex-specific (or IFNA-specific) regulatory windows, and perform Fisher's exact tests (FETs) on results.
##############################################
#In unix, run the HOMER command, findMotifsGenome.pl
cd ~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation
genome=hg38.fa
c=dex #or ifn
findMotifsGenome.pl ${c}_reg_specific23aug2021.bed $genome ${c}_specific_reg23aug2021 -size given -bg windows1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeat${c}.txt -mknown ~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/known.motifs -nomotif

#In R, perform Fisher's Exact Tests on HOMER output
setwd('~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation')
h<-read.table('knownResults_dex.txt',sep='\t',header=T)
h$a<-1;h$b<-1;h$c<-1;h$d<-1;h$or<-1;h$loweror<-1;h$upperor<-1;h$p<-1
head(h)
for (p in 1:nrow(h)){
  out<-h[p,]
  bgwmotif<-out$no_of_Background_Sequences_with_Motifof_242395
  bgwomotif<- 242395 - bgwmotif #Number bg windows with motif, 242395, obtained from HOMER output file 'knownResults.txt'
  setwmotif<-out$no_of_Target_Sequences_with_Motifof_1131
  setwomotif<-1131-setwmotif #Number test windows with motif, 1131, obtained from HOMER output file 'knownResults.txt'
  a <- bgwomotif
  b <- bgwmotif
  c <- setwomotif
  d <- setwmotif
  values = matrix(c(a,b,c,d), nrow = 2);values
  this<-fisher.test(values,alternative="two.sided");this
  h[p,'a']=a
  h[p,'b']=b
  h[p,'c']=c
  h[p,'d']=d
  h[p,'p']=this$p.value
  h[p,'or']=this$estimate
  h[p,'loweror']=this$conf.int[1]
  h[p,'upperor']=this$conf.int[2]
} #close p

##############################################
#Full mRNA-Seq pipeline and modeling
##############################################
################
#Submit slurm job 'trimandmap.sh' to trim and map endogenous RNA
################
cd ~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation
mkdir mapped
for f in L31416 L31417 L31418 L31419 L31420 L31421 L31422 L31423 L31424 L31425 L31426 L31427 L31428 L31429 L31430 L31431 L31432 L31433; do 
cat trimandmap.sh | sed -e s/FILEINFO/${f}/g > trimandmap_${f}.sh;chmod 777 trimandmap_${f}.sh;sbatch trimandmap_${f}.sh;done

#!/bin/bash
#SBATCH -n 12 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machinel
#SBATCH -t 0-96:00 # Runtime in D-HH:MM
#SBATCH --mem-per-cpu=4G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o trim_FILEINFO.out # File to which STDOUT will be written
#SBATCH -e trim_FILEINFO.err # File to which STDERR will be written
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=youremail@gmail.com # Email to which notifications will be sent

module load cutadapt
module load TrimGalore
module load STAR/2.5.0c-fasrc01
genome_path=hg38_STAR #Directory containing STAR-indexed hg38 reference and associated files (Genome, hg38.fa, SA, SAindex)

#Trim low quality basepairs (phred score < 20) and adapters
trim_galore -q 20 -o ../trimmed --fastqc -a AGATCGGAAGAGC --stringency 2 --length 25 --paired FILEINFO_R1.fastq.gz FILEINFO_R2.fastq.gz
#Map reads using STAR
STAR --genomeDir $genome_path/ --outFileNamePrefix mapped/mapped_FILEINFO --runThreadN 12 --readFilesCommand zcat --readFilesIn FILEINFO_R1_val_1.fq FILEINFO_R2_val_2.fq

################
#Submit slurm job 'spliceindex.sh' below to create index of splice locations for STAR 2nd pass mapping:
################
cd ~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation
mkdir secondpass_hg38_15jan2021
sbatch spliceindex.sh

#!/bin/bash
#SBATCH -n 24 # Number of cores
#SBATCH -N 2 # Ensure that all cores are on one machine
#SBATCH -t 0-96:00 # Runtime in D-HH:MM
#SBATCH --mem=96G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o spliceindex.out # File to which STDOUT will be written
#SBATCH -e spliceindex.err # File to which STDERR will be written
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=youremail@gmail.com # Email to which notifications will be sent

module load STAR/2.5.0c-fasrc01
path_out=~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation/secondpass_hg38_15jan2021
splicelocationsdir=~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation/secondpass_hg38_15jan2021/mapped
genome_path=hg38_STAR
STAR --runMode genomeGenerate --outFileNamePrefix $path_out/second --runThreadN 24 --genomeDir $path_out/ --limitGenomeGenerateRAM 96000000000 --genomeFastaFiles \
$genome_path/hg38.fa --sjdbFileChrStartEnd $splicelocationsdir/mapped_L31416SJ.out.tab $splicelocationsdir/mapped_L31417SJ.out.tab \
$splicelocationsdir/mapped_L31418SJ.out.tab $splicelocationsdir/mapped_L31419SJ.out.tab $splicelocationsdir/mapped_L31420SJ.out.tab \
$splicelocationsdir/mapped_L31421SJ.out.tab $splicelocationsdir/mapped_L31422SJ.out.tab $splicelocationsdir/mapped_L31423SJ.out.tab \
$splicelocationsdir/mapped_L31424SJ.out.tab $splicelocationsdir/mapped_L31425SJ.out.tab $splicelocationsdir/mapped_L31426SJ.out.tab \
$splicelocationsdir/mapped_L31427SJ.out.tab $splicelocationsdir/mapped_L31428SJ.out.tab $splicelocationsdir/mapped_L31429SJ.out.tab \
$splicelocationsdir/mapped_L31430SJ.out.tab $splicelocationsdir/mapped_L31431SJ.out.tab $splicelocationsdir/mapped_L31432SJ.out.tab \
$splicelocationsdir/mapped_L31433SJ.out.tab

################
#Submit slurm job 'map2.sh' to perform second pass mapping, using the splice index created from first pass mapping.
################
mkdir ~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation/mapped/secondpass

for f in L31416 L31417 L31418 L31419 L31420 L31421 L31422 L31423 L31424 L31425 L31426 L31427 L31428 L31429 L31430 L31431 L31432 L31433; 
do cat map2.sh | sed -e s/FILEINFO/$f/g > map2_$f.sh; sbatch map2_$f.sh; done

#!/bin/bash
#SBATCH -n 12 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-96:00 # Runtime in D-HH:MM
#SBATCH --mem-per-cpu=4G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o map2_FILEINFO.out # File to which STDOUT will be written
#SBATCH -e map2_FILEINFO.err # File to which STDERR will be written
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=youremail@gmail.com # Email to which notifications will be sent

module load STAR/2.5.0c-fasrc01
genome_path=~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation/secondpass_hg38_15jan2021
path_out=~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation/mapped/secondpass
STAR --genomeDir $genome_path/ --outFileNamePrefix $path_out/mapped2_FILEINFO --runThreadN 12 --readFilesCommand zcat --readFilesIn ../trimmed/FILEINFO_R1_val_1.fq ../trimmed/FILEINFO_R2_val_2.fq

################
#Submit slurm job 'htseq.sh' below to count reads uniquely mapped to genes
################
for f in L31416 L31417 L31418 L31419 L31420 L31421 L31422 L31423 L31424 L31425 L31426 L31427 L31428 L31429 L31430 L31431 L31432 L31433;
do cat htseq.sh | sed -e s/FILEINFO/$f/g > htseq_$f.sh; sbatch htseq_$f.sh; done

#!/bin/bash
#SBATCH -n 2 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-96:00 # Runtime in D-HH:MM
#SBATCH --mem=8G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o htseq_FILEINFO.out # File to which STDOUT will be written
#SBATCH -e htseq_FILEINFO.err # File to which STDERR will be written
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=youremail@gmail.com # Email to which notifications will be sent

f=FILEINFO
module load samtools
module load HTSeq
cd ~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation/mapped/secondpass
pathgtf=~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation
samtools view -H mapped2_${f}Aligned.out.sam > head_${f}.sam
awk '($5=="255") {print $0;}' mapped2_${f}Aligned.out.sam > select_${f}.sam
cat head_${f}.sam select_${f}.sam > ${f}_uniquelymapped.sam 
samtools view -bS ${f}_uniquelymapped.sam | samtools sort -n - -o ${f}_uniquesortedbychr.bam
samtools view ${f}_uniquesortedbychr.bam | htseq-count --mode=union --stranded=no --idattr=gene_id - ${pathgtf}/Homo_sapiens.GRCh38.102_wchr.gtf > ${f}_HTSeq-counts_uniq.txt

################
#Create htseq count file containing all samples
################
setwd('~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation')
sequence=c('L31416','L31417','L31418','L31419','L31420','L31421','L31422','L31423','L31424','L31425',
                 'L31426','L31427','L31428','L31429','L31430','L31431','L31432','L31433')
for (p in sequence) {
  perm<-paste0(p,"_HTSeq-counts_uniq.txt",sep="")
  permfile<-read.table(perm,header=F,sep="\t")
  
  if(p=="L31416")
  {
    h <- permfile
  } else {
    h <- cbind(h,permfile)
  }
} #close p

rownames(h)<-h$V1
h2<-h[,seq(2,ncol(h),by=2)]
colnames(h2)<-sequence
head(h2)
write.table(h2,file='allmstarr_htseq_winfo.txt',sep='\t')
tail(h2)
dim(h2)
#Remove extra info:
h3 <- h2[-c(60676:60680),]
colSums(h3)/colSums(h2)
h2['__no_feature',]/colSums(h2)
write.table(h3,file='allmstarr_htseq.txt',sep='\t')

################
#Calculate average transcript length per gene (includes exons only)
################
bedtools nuc -fi hg38.fa -bed Homo_sapiens.GRCh38.102_wchr.gtf > outfile.gtf.metrics
# Subset to exons:
awk '($3=="exon")' outfile.gtf.metrics > outfile.gtf.metrics_exon
#Count length of each transcript (for which we're only including exons). Specify input file. 
python ~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation/calculateTranscriptLength.py outfile.gtf.metrics_exon 
# Take the average transcript length for each gene. Specify input file.
python ~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation/calculateMeanGeneLength.py outfile.gtf.metrics_transcriptLengths.txt 
#NOTE: Last gene is printed twice, so remove last line.

################
#Prepare traits and rawcounts files
################
library(coxme)
library(edgeR)
library(limma)
library(spaa)
library("MASS")
library(ggplot2)
library(EMMREML)
library(stats)
library(cobs)
library(matrixcalc)
library(dplyr)
library(qvalue)

setwd('~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation')
traits<-read.table("traits_mstarr_endorna.txt",sep="\t",header=T)
head(traits)

rawcounts=read.table("allmstarr_htseq.txt",sep="\t",header=T)
rawcounts<-rawcounts[,as.character(traits$library)]
#Calculate TPM
gene_length<-read.table("hg38.102_wchr_gtf.metrics_transcriptLengths_meanGeneLength.txt",header=T)
gene_length<-gene_length[,c("gene_id","mean_gene_length")]
counts2<-merge(rawcounts,gene_length,by.x="row.names",by.y="gene_id")
rownames(counts2)<-counts2[,1];counts2<-counts2[,-c(1)]
total<-colSums(counts2[,1:nrow(traits)])

tpm<-counts2[,1:nrow(traits)]*100*10^6/counts2$mean_gene_length

#Set final row equal to T
tpm[nrow(tpm)+1,]<-apply(X = counts2[,1:nrow(traits)],2,function(x) sum(x*100/counts2$mean_gene_length))
for (i in 1:nrow(traits)){
  tpm[,i]<-tpm[,i]/tpm[nrow(tpm),i]
}
tpm<-tpm[-nrow(tpm),]

#Require at least 3 of 6 samples in any condition to have TPM>=X
counts2<-counts2[,-ncol(counts2)] #remove last column, which is mean_gene_length

null<-subset(traits,treatment=='null')
countsnull<-tpm[,null$library]
countsnull2<-countsnull[which(apply(countsnull,1,function(a){quantile(a,0.5)})>=3),]
dim(countsnull2)

dex<-subset(traits,treatment=='dex')
countsdex<-tpm[,dex$library]
countsdex2<-countsdex[which(apply(countsdex,1,function(a){quantile(a,0.5)})>=3),]
dim(countsdex2)

IFN<-subset(traits,treatment=='IFN')
countsIFN<-tpm[,IFN$library]
countsIFN2<-countsIFN[which(apply(countsIFN2,1,function(a){quantile(a,0.5)})>=3),]
dim(countsIFN2)
keep<-c(rownames(countsnull2),rownames(countsdex2),rownames(countsIFN2))
head(keep)
length(keep)
keep<-unique(keep) #17197
rawcounts2<-counts2[keep,]
dim(rawcounts2)

#Reduce to protein-coding genes
symbols<-read.table('genesymbols_GRCh38.p13_18jan2021.txt',sep='\t',header=T)
raw<-merge(rawcounts2,symbols,by.x='row.names',by.y='Gene.stable.ID')
dim(rawcounts2)
dim(raw)
raw2<-subset(raw,Gene.type=='protein_coding')
dim(raw2)
head(raw2)
rownames(raw2)<-raw2$Row.names
raw2<-raw2[,-c(1)]
raw2<-raw2[,1:18]
head(raw2)
rawcounts2<-raw2
dim(raw2)
write.table(rawcounts2,file='rawcounts_protcoding_3sampstpm3_anycondition.txt',sep='\t')

################
#Voom normalize and run model and 100 permutations
################
library(coxme)
library(edgeR)
library(limma)
library(spaa)
library("MASS")
library(ggplot2)
library(EMMREML)
library(stats)
library(cobs)
library(matrixcalc)
library(dplyr)
library(qvalue)

setwd('~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation')
traits<-read.table("traits_mstarr_endorna.txt",sep="\t",header=T)
traits<-subset(traits,treatment!='IFN')

rawcounts2<-read.table('rawcounts_protcoding_3sampstpm3_anycondition.txt',sep='\t',header=T)
rawcounts2<-rawcounts2[,traits$library]
traits$meth_condition<-as.factor(traits$meth_condition)
traits$meth_condition<-relevel(traits$meth_condition,'sham') #sham is reference
traits$treatment<-as.factor(traits$treatment)
traits$treatment<-relevel(traits$treatment,'null') 

perms<-as.data.frame(matrix(nrow=12,ncol=100))

for (p in 1:100) {
  traits$permtreatment<-sample(traits$treatment,replace=F)
  perms[,p]<-traits$permtreatment
  design<-model.matrix(~traits$meth+traits$permtreatment+traits$mapped_no_feature)
  design
  
  dge <- DGEList(counts=rawcounts2)
  dge <- calcNormFactors(dge)
  v <- voom(dge,design,plot=T)
  fit <-lmFit(v,design)
  fit2 <- eBayes(fit)
  se.coef <- sqrt(fit2$s2.post) * fit2$stdev.unscaled
  fit2 <- as.data.frame(cbind(fit2$p.value[,2:4],fit2$coefficient[,2:4],se.coef[,2:4]))
  #colnames(fit2)<-c('meth_p','IFN_p','nofeat_p','meth_beta','IFN_beta','nofeat_beta','meth_se','IFN_se','nofeat_se')
  colnames(fit2)<-c('meth_p','dex_p','nofeat_p','meth_beta','dex_beta','nofeat_beta','meth_se','dex_se','nofeat_se')
  write.table(fit2,paste0('model_meth_permdex_nofeat_',p,'.txt'),sep='\t')
  print(p)
} #close p

head(perms)
write.table(perms,file='perms100dex20aug2021.txt',sep='\t')

# fit2 <- as.data.frame(cbind(fit2$p.value[,2:5],fit2$coefficient[,2:5],se.coef[,2:5]))
# colnames(fit2)<-c('meth_p','dex_p','ifn_p','nofeat_p','meth_beta','dex_beta','ifn_beta','nofeat_beta','meth_se','dex_se','ifn_se','nofeat_se')
head(fit2)

#Model the actual observed data. Here, we want to make sure there's no effect on the data of 1) mSTARR-seq methylation or 2) efficiency of mRNA selection during library prep.
design<-model.matrix(~traits$meth+traits$treatment+traits$mapped_no_feature)
design
dge <- DGEList(counts=rawcounts2)
dge <- calcNormFactors(dge)
v <- voom(dge,design,plot=T)
fit <-lmFit(v,design)
fit2 <- eBayes(fit)
se.coef <- sqrt(fit2$s2.post) * fit2$stdev.unscaled
fit2 <- as.data.frame(cbind(fit2$p.value[,2:4],fit2$coefficient[,2:4],se.coef[,2:4]))
#colnames(fit2)<-c('meth_p','IFN_p','nofeat_p','meth_beta','IFN_beta','nofeat_beta','meth_se','IFN_se','nofeat_se')
colnames(fit2)<-c('meth_p','dex_p','nofeat_p','meth_beta','dex_beta','nofeat_beta','meth_se','dex_se','nofeat_se')
write.table(fit2,paste0('model_meth_dex_nofeat.txt'),sep='\t')
################
#Correct for multiple testing. You don't need to change any code for this section; it just defines the perm.fdr function.
################
library(cobs)
perm.fdr=function(input_df,perm_df,Pvals_col_name,name){
  
  pvals_index=which(colnames(input_df)==Pvals_col_name)
  ro<-input_df[order(input_df[,pvals_index]),]
  p_obs <- data.frame(pvalue=ro[,pvals_index])
  p_vector<-matrix(as.matrix(perm_df),ncol=1)
  p_vector=data.frame(p_vector[order(p_vector)])
  
  F<-p_obs[,1]
  F_o<-p_obs[,1]
  pi_hat<-p_obs[,1]
  
  j=1
  observed<-length(p_obs[,1])
  randoms<-length(p_vector[,1])
  
  for(i in 1:observed)
  {
    repeat
    {
      if((p_vector[j,1]<p_obs[i,1])&j<randoms){j<-j+1}else{break}
    }
    F[i]=i/observed
    F_o[i]=(j-1)/randoms
    if(F_o[i]<1){pi_hat[i]=(1-F[i])/(1-F_o[i])}else{pi_hat[i]=1}
  }
  tabla <-data.frame(pi_hat,pval=p_obs[,1])
  
  tabla[1,]=c(1,0)
  last_percentile_average=mean(tabla$pi_hat[as.integer(min((length(tabla[,1])*0.99),(nrow(tabla)-1)):length(tabla[,1]))])
  tabla[nrow(tabla),]=c(last_percentile_average,1)
  constraint_matrix=as.matrix(data.frame(c(0,2),c(0,1),c(1,0)))
  f_hat<-suppressWarnings(cobs(tabla$pval,tabla$pi_hat,constraint="convex",pointwise=constraint_matrix,maxiter=1000,print.warn=FALSE,print.mesg=FALSE))
  
  f_hat_serie=f_hat$fitted
  pi_o=f_hat_serie[length(f_hat_serie)]
  pi_o=min(pi_o,1)
  pi_o=max(pi_o,0)
  
  Fdr_ST_perm=pi_o*F_o/F
  
  for(i in 1:length(p_obs[,1]))
  {
    Fdr_ST_perm[i]=pi_o*F_o[i]/F[i]
    if(i>1)
    {
      for(j in 1:(i-1))
      {
        if(Fdr_ST_perm[i-j]>Fdr_ST_perm[i]){Fdr_ST_perm[i-j]=Fdr_ST_perm[i]}else{break}
      }
    }
    if(Fdr_ST_perm[i]>1)  Fdr_ST_perm[i]=1
  }
  
  fdrs_df <-data.frame(ro,q_ST_perm=Fdr_ST_perm)
  rownames(fdrs_df)=rownames(ro)
  colnames(fdrs_df)[ncol(fdrs_df)]=paste0("fdr_",name)
  
  return(fdrs_df)
}

################
#Apply function perm.fdr to calculate the empirical FDR from the 100 permutations.
################
for (p in 1:100) {
  perm<-paste0("model_meth_permdex_nofeat_",p,'.txt',sep="")
  permfile<-read.table(perm,header=T,sep="\t")
  
  if(p=="1")
  {
    shuffled_treatment <-data.frame(x=permfile[,"dex_p"])
    rownames(shuffled_treatment)=rownames(permfile)
  } else {
    shuffled_treatment <- cbind(shuffled_treatment,x=permfile[,"dex_p"])
  }
}

real<-read.table("model_meth_dex_nofeat.txt",header=T,sep="\t")
res_full=perm.fdr(data.frame(real),shuffled_treatment,"dex_p","dex")
res_full=res_full[order(rownames(res_full)),]
write.table(res_full,file="model_meth_dex_nofeat_FDR.txt",sep="\t")
##############################################
#Ask, If genes are nearest to the IFNA-specific regulatory windows (as opposed to genes nearest IFNA regulatory windows *and shared* in the baseline and/or dex treatment), do they show a higher GE response to IFNA?
#Note: Repeat this analysis, switching out IFNA with dex.
##############################################
#In Unix, obtain closest genes to all tested windows in IFN mstarr dataset
cd ~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation
bedtools closest -b human_genes20mar2020_GRCh38.p13.bed -a windows1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatIFN.txt -d -t first > ifnwindowsdistancetonearestTSS.txt 

#In R, obtain the enhancers from the 3 treatments (baseline, IFNA, dex)
setwd('~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation')
real<-read.table('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwol31250_FDR2.txt',header=T,sep='\t')
real_IFN<-read.table('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatIFN_FDR2.txt',header=T,sep='\t')
real_dex<-read.table('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatdex_FDR2.txt',header=T,sep='\t')

cc<-0.01
m<-subset(real,(unmeth.rna_beta>0 & fdr_RNAinunmeth<cc) | (meth.rna_beta>0 & fdr_RNAinmeth<cc))
m_IFN<-subset(real_IFN,(unmeth.rna_beta>0 & fdr_RNAinunmeth<cc) | (meth.rna_beta>0 & fdr_RNAinmeth<cc))
m_dex<-subset(real_dex,(unmeth.rna_beta>0 & fdr_RNAinunmeth<cc) | (meth.rna_beta>0 & fdr_RNAinmeth<cc))

#Read in dataframe d. d is all tested windows (for the treatment of interest), and their closest genes:
d<-read.table('ifnwindowsdistancetonearestTSS.txt',sep='\t')
head(d)
rownames(d)<-d$V4
d<-d[,c('V4','V9','V10')]

#add nearest genes to regulatory windows
m_IFNd<-merge(m_IFN,d,by='row.names')
head(m_IFNd)

#reduce regulatory windows to those that had nearest genes that were *tested* in the endogenouse gene expression dataset
model<-read.table('model_meth_ifn_nofeat_methFDR_ifnFDR.txt',sep='\t',header=T)
nonifnge<-subset(model, fdr_IFN_p > 0.01)
ifnge<-subset(model, fdr_IFN_p < 0.01)
modelsign<-ifnge
m_IFNd<-m_IFNd[m_IFNd$V9 %in% rownames(model),]
head(m_IFNd)

#Add column indicating whether gene significantly responded to IFNA stimulation or not
m_IFNd$geresponse<-ifelse(m_IFNd$V9 %in% rownames(ifnge),1,0)
rownames(m_IFNd)<-m_IFNd$Row.names

#Separate regulatory windows according to whether they are IFNA-specific or shared across conditions
mnull_mdex<-merge(m,m_dex,all=T,by='row.names') #Creates list of all regulatory windows in baseline and/or dex treatments
dim(mnull_mdex)
names=sort(unique(mnull_mdex$Row.names))
length(names)
IFN_spec <- m_IFNd[!rownames(m_IFNd) %in% names, ]
IFNreg_notspec<- m_IFNd[!rownames(m_IFNd) %in% rownames(IFN_spec), ]

#Reduce results to unique genes 
IFN_spec<-unique(IFN_spec[,c('V9','geresponse')])
IFNreg_notspec<-unique(IFNreg_notspec[,c('V9','geresponse')])

#Obtain each gene's nearest window
names<-unique(m_IFNd$V9)
result<-as.data.frame(matrix(nrow=length(names),ncol=ncol(m_IFNd)))

for (i in 1:length(names)) {
  out<-subset(m_IFNd,V9==names[i])
  out<-out[which(out$V10==min(out$V10)),]
  result[i,]<-as.matrix(out[sample(nrow(out),1),]) #Take closest. If there are multiple frags with same distance, randomly take 1, noting that this will lead to minor stochasticity in result. 
}
colnames(result)<-colnames(out)

#Only keep window-gene pairs within 100kb distance of eachother
result$V10<-as.numeric(result$V10)
result<-subset(result,V10<=100000)
dim(result)
#1425 one-to-one reg_window-geneGE closest pairs. 1350 when reducing to less than 100,000kb.

#Merge info on regulatory activity and matched nearest gene expression
colnames(model)<-paste0('GE_',colnames(model))
head(model)
result2<-merge(result,model,by.x='V9',by.y='row.names')

result2$ifn_spec<-ifelse(result2$Row.names %in% rownames(IFN_spec),'blue','black')
result2$ifn_spec<-as.factor(result2$ifn_spec)
spec<-subset(result2,ifn_spec=='blue')
nonspec<-subset(result2,ifn_spec=='black')
tt<-t.test(spec$GE_IFN_beta,nonspec$GE_IFN_beta,paired=F)
tt

##############################################
#1. Reduce IFNA-specific enhancers to those with empirically supported binding by ISGF3, the transcription factor that binds to ISRE motifs.
#2. Link this enhancer subset to target genes. 
#3. Test whether the gene expression response to IFNA stimulation is different between the gene targets of IFNA-specific enhancers vs gene targets of IFNA enhancers shared with baseline and/or dex conditions. 
##############################################
#1.Narrow in on joint intersection of STAT1 *AND* STAT2 first (The 2 of 3 components of ISGF3 for which CHIP data are available).
cd ~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation/
bedtools intersect -a ENCFF478XGE.bed -b ENCFF394KTR.bed > STAT1andSTAT2.bed
bedtools intersect -a ifn_reg_specific23aug2021.bed -b STAT1andSTAT2.bed -wa > ifn_reg_specific_STAT1andSTAT2.bed

#2. Using the EnhancerAtlas linkage (http://www.enhanceratlas.org/data/AllEPs/hs/K562_EP.txt accessed Sept 10, 2021), relate this enhancer subset to their linked genes (used liftover at http://genome.ucsc.edu/cgi-bin/hgLiftOver).
cd ~/Desktop/tunglab/capture-nrsa/isre/
bedtools intersect -a ifn_reg_specific_STAT1andSTAT2.bed -b k562gene_enhancer_hg38.bed -wb > ifn_reg_specific_STAT1andSTAT2_linkedgene.txt
bedtools intersect -a ifn_allreg_windows_hg38.bed -b ~/Desktop/tunglab/capture-nrsa/isre/k562gene_enhancer_hg38.bed -wb > ifn_allreg_windows_hg38linkedgene.txt

#3. Test whether the gene expression response to IFNA stimulation is different between 
#the gene targets of IFNA-specific enhancers vs gene targets of IFNA enhancers shared under baseline and/or dex conditions. 
setwd('~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation')
ge<-read.table('model_meth_ifn_nofeat_methFDR_ifnFDR.txt')
gg<-read.table('ifn_reg_specific_STAT1andSTAT2_linkedgene.txt')
gg<-gg[,c('V8','V4')]
gg<-gg[!duplicated(gg$V8),]
ge2<-merge(ge,gg,by.x='row.names',by.y='V8')
genot <- ge[!rownames(ge) %in% gg$V8, ]

#boxplot
rownames(ge2)<-ge2[,1]
ge2<-ge2[,-c(1)]
ge2<-ge2[,-c(12)]
ge2$isretarget <- 1
genot$isretarget<-0

ggg<-rbind(ge2,genot)
ggg$isretarget<-as.factor(ggg$isretarget)
#write.table(ggg,file='endogenous_ge_IFNAeffect.txt',sep='\t')
head(ge2)
head(genot)
dim(ge2);dim(genot)

tt<-t.test(ge2$IFN_beta,genot$IFN_beta,paired=F) #IFNA-specific (and ISRE containing) vs IFNA-shared
tt
##############################################
#Comparison of prop MD in IFNA-specific vs baseline-specific enhancers (2-sided binomial test)
##############################################
setwd('~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation')
n<-read.table('null_reg_specific23aug2021.bed')
d<-read.table('dex_reg_specific23aug2021.bed')
i<-read.table('ifn_reg_specific23aug2021.bed')

#Read in nested model results for each condition.
nn<-read.table('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwol31250_FDR2.txt')
dd<-read.table('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatdex_FDR2.txt')
ii<-read.table('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatifn_FDR2.txt')

#Subset nested model results to regulatory windows.
cc=0.01
nn<-subset(nn,(unmeth.rna_beta>0 & fdr_RNAinunmeth<cc) | (meth.rna_beta>0 & fdr_RNAinmeth<cc))
ii<-subset(ii,(unmeth.rna_beta>0 & fdr_RNAinunmeth<cc) | (meth.rna_beta>0 & fdr_RNAinmeth<cc))
dd<-subset(dd,(unmeth.rna_beta>0 & fdr_RNAinunmeth<cc) | (meth.rna_beta>0 & fdr_RNAinmeth<cc))

#Create dataframes of windows with MD regulatory activity
nmd<-subset(nn,fdr_int<cc)
imd<-subset(ii,fdr_int<cc)
dmd<-subset(dd,fdr_int<cc)

#Prop MD in ALL ifna reg vs ALL null reg: 
nrow(imd)
nrow(ii)
nrow(imd)/nrow(ii)
nrow(nmd)
nrow(nn)
nrow(nmd)/nrow(nn)
xx<-binom.test(3463,4632, p = 0.4751411, alternative = c("two.sided"), conf.level = 0.95)
xx$p.value

#Subset nested model results to those that are condition-specific
nn2<-nn[rownames(nn) %in% n$V4, ]
dd2<-dd[rownames(dd) %in% d$V4, ]
ii2<-ii[rownames(ii) %in% i$V4, ]

#Create dataframes of windows with MD regulatory activity
nmd2<-subset(nn2,fdr_int<cc)
imd2<-subset(ii2,fdr_int<cc)
dmd2<-subset(dd2,fdr_int<cc)

#Prop MD in IFNA-specific reg vs null-specific reg:
nrow(imd2)
nrow(ii2)
nrow(imd2)/nrow(ii2)
nrow(dmd2)/nrow(dd2)
nrow(nmd2)/nrow(nn2)
xx2<-binom.test(1341,1614, p = 0.4640235, alternative = c("two.sided"), conf.level = 0.95)
xx2$p.value

##############################################
#Assess whether regulatory regions that harbor TFBMs for TFs central to the interferon response (ISRE, IRF1, IRF2, IRF3, IRF4, IRF8; N = 663 windows) show methylation dependence
##############################################
# #Pull coordinates of TFBMs for ISRE, IRF1, IRF2, IRF3, IRF4, IRF8 from HOMER database bed file.
# #Note, homer.KnownMotifs.hg38.191020.bed is large (28Gb) and therefore not provided here; the output ifn_reg_specific_w_ISRE_IRF_nonredundant.txt is provided.
# grep -w "ISRE(IRF)" homer.KnownMotifs.hg38.191020.bed > ISRE_homer_hg38.bed
# grep -w "T1ISRE(IRF)" homer.KnownMotifs.hg38.191020.bed > T1ISRE_homer_hg38.bed
# grep -w "IRF1(IRF)" homer.KnownMotifs.hg38.191020.bed > IRF1_homer_hg38.bed
# grep -w "IRF2(IRF)" homer.KnownMotifs.hg38.191020.bed > IRF2_homer_hg38.bed
# grep -w "IRF3(IRF)" homer.KnownMotifs.hg38.191020.bed > IRF3_homer_hg38.bed
# grep -w "IRF4(IRF)" homer.KnownMotifs.hg38.191020.bed > IRF4_homer_hg38.bed
# grep -w "IRF8(IRF)" homer.KnownMotifs.hg38.191020.bed > IRF8_homer_hg38.bed
# cat ISRE_homer_hg38.bed T1ISRE_homer_hg38.bed IRF1_homer_hg38.bed IRF2_homer_hg38.bed IRF3_homer_hg38.bed IRF4_homer_hg38.bed IRF8_homer_hg38.bed > ISRE_IRF_homer_hg38.bed
# 
# bedtools intersect -a ISRE_IRF_homer_hg38.bed -b ifn_reg_specific23aug2021.bed -f 1 -wb > ifn_reg_specific_w_ISRE_IRF.bed 
# 
# #In R, remove redundancy in regions.
# setwd('~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation')
# ii<-read.table('~/Desktop/tunglab/capture-nrsa/mSTARR_data/counts/counts_dupped/motifs/ifn_spec_enhancers_isre_irf/ifn_reg_specific_w_ISRE_IRF.txt',header=F)
# head(ii)
# ii<-ii$V10
# ii<-unique(ii) 
# write.table(ii,'ifn_reg_specific_w_ISRE_IRF_nonredundant.txt')

#Read in IFN nested model regults, and IFNA-specific regulatory regions harboring transcription factor binding motifs (TFBMs) for ISRE, IRF1, IRF2, IRF3, IRF4, or IRF8.
library('dplyr')
setwd('~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation')
i<-read.table('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatIFN_FDR2.txt',header=T)
head(i)
ii<-read.table('ifn_reg_specific_w_ISRE_IRF_nonredundant.txt')

#Number IFNA-specific regulatory regions harboring TFBMs for ISRE, IRF1, IRF2, IRF3, IRF4, or IRF8.
nrow(ii)

#Proportion regulatory windows with more regulatory activity in unmethylated condition, for IFNA-specicic regulatory regions with ISRE or IRF1, IRF2, IRF3, IRF4, IRF8.
i2 <- i[rownames(i) %in% ii$V1,] #IFNA-specicic regulatory regions with ISRE, IRF1, IRF2, IRF3, IRF4, IRF8
i2_int<-subset(i2,fdr_int<0.01)
nrow(i2_int) #562 with significant MD reg activity 
nrow(i2_int)/nrow(i2) #84.76621% IFN-reg-spec windows with ISRE/IRFs show MD
length(which(i2_int$int_beta<0)) #561/562 MD IFN-reg-spec windows show more regulatory activity in the unmethylated than methylated condition.

##############################################
#For 1005 IFNA regulatory windows containing ISRE motifs, test for methylation dependence of those windows in the null condition and in the IFNA condition using paired t-test of MD betas.
##############################################
#In R, create file of all regions tested for regulatory activity in the baseline (null) condition and in the IFNA condition.
setwd('~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation')
n<-read.table('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwoL31250.txt',sep='\t')
i<-read.table('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatIFN.txt',sep='\t')
head(n)
head(i)
ni<-merge(n,i,by='row.names')
rownames(ni)<-ni$Row.names

#Save and manually create bed file
write.table(ni,file='ni_nested_13sept2021.bed.txt',sep='\t',quote=F,row.names=F,col.names=F)

#In unix, obtain null+IFNA tested regions that have ISRE element in K562s
cd ~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation
bedtools intersect -a ni_nested_13sept2021.bed.txt -b ENCFF478XGE.bed -b ENCFF394KTR.bed -wa > ni_alltested_STAT1orSTAT2.bed

#In R, separate ni by whether region has ISRE element
setwd('~/Desktop/github_mstarr/3_Env_perturbation_reveals_cryptic_regulatory_elements_and_cryptic_effects_of_DNA_methylation')
n<-read.table('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwoL31250.txt',sep='\t')
i<-read.table('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatIFN.txt',sep='\t')
head(n)
head(i)
ni<-merge(n,i,by='row.names')
rownames(ni)<-ni$Row.names

isre<-read.table('ni_alltested_STAT1orSTAT2.bed',sep='\t')
head(isre)
dim(isre)
ni_isre <- ni[rownames(ni) %in% isre$V4, ]
ni_notisre <- ni[!rownames(ni) %in% isre$V4, ]
#Check that they all add up:
dim(ni_isre)
dim(ni_notisre)
dim(ni)

#Within tested regions with ISRE, is there a MD effect in *null*?
tt<-t.test(ni_isre$unmeth.rna_beta.x,ni_isre$meth.rna_beta.x,paired=T)
tt

#Within tested regions with ISRE, is there a MD effect in *IFNA*?
tt2<-t.test(ni_isre$unmeth.rna_beta.y,ni_isre$meth.rna_beta.y,paired=T)
tt2
tt2$p.value







                                              