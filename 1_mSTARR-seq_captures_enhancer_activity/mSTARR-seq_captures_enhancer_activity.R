##############################################
#SECTION: mSTARR-seq_captures_enhancer_activity
##############################################

##############################################
#trimmap.sh Trim and map mSTARR-seq RNA and mSTARR-Seq DNA reads to hg38. Create and submit multiple jobs (our computer cluster uses slurm) by replacing FILEINFO with sample ids from ids.txt, using the command:
for f in `cat ids.txt`; do cat trimmap.sh | sed -e s/FILEINFO/$f/g > trimmap_$f.sh; sbatch trimmap_$f.sh; done
##############################################
#!/bin/bash
#SBATCH -n 8 # Number of cores
#SBATCH -t 0-96:00 # Runtime in D-HH:MM
#SBATCH --mem=32G #Memory for job
#SBATCH -o trimmap_FILEINFO.out # File to which STDOUT will be written
#SBATCH -e trimmap_FILEINFO.err # File to which STDERR will be written
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=myemail@gmail.com # Email to which notifications will be sent

module load TrimGalore
module load cutadapt
module load bwa
module load samtools
module load bedtools2

inputdir=directory_containing_fastq
cd ${inputdir}
trim_galore -q 20 -o ${inputdir} -a AGATCGGAAGAGC --stringency 2 --length 25 --paired ${inputdir}/FILEINFO_R1.fastq.gz ${inputdir}/FILEINFO_R2.fastq.gz
bwa mem -t 8 hg38.fa FILEINFO_R1_val_1.fq.gz FILEINFO_R2_val_2.fq.gz > FILEINFO.sam

samtools view -bSq 10 -f 0x2 FILEINFO.sam | samtools sort -n - | bedtools bamtobed -bedpe | awk '{OFS="\t"; print $1,$2,$6,$7}' > FILEINFO_pe.bed
awk '($3-$2) >= 0 && ($3-$2) <= 2000' FILEINFO_pe.bed > FILEINFO_pe_filt.bed #Remove lines in which fragment is negative value or greater than 2000bp (rare, but occasionally happens)
bedtools coverage -a ~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/hg38_600bp_windows.bed -b FILEINFO_pe_filt.bed > counts_FILEINFO_600bpwin.txt 
awk '{OFS="\t"; print $4}' counts_FILEINFO_600bpwin.txt > counts2_FILEINFO_600bpwin.txt

##############################################
#Subset to unique fragments for comparison of library diversity to other datasets
for f in `cat ids.txt`; do cat uniqueit.sh | sed -e s/FILEINFO/$f/g > uniqueit_$f.sh; sbatch uniqueit_$f.sh; done
##############################################
#!/bin/bash
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-96:00 # Runtime in D-HH:MM
#SBATCH --mem=4G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o uniqueit_FILEINFO.out # File to which STDOUT will be written
#SBATCH -e uniqueit_FILEINFO.err # File to which STDERR will be written
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=myemail@gmail.com # Email to which notifications will be sent

inputdir=directory_containing_bedfiles
cd ${inputdir}
sort -k1,1 -k2,2n -k3,3n -u FILEINFO_pe_filt.bed > FILEINFO_uniq_pe_filt.bed

#After completing uniqueit.sh on all bed files:
wc -l *_uniq_pe_filt.bed > linecounts.txt 

##############################################
#In R, create counts file for mSTARR-seq RNA and mSTARR-Seq DNA reads
##############################################
b<-read.table('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/hg38_600bp_windows.txt',header=F,sep='\t')
samps<-read.table('ids.txt',sep='\t')
samps<-samps$V1
for (p in 1:length(samps)) {
  myfile<-paste0('counts2_',samps[p],'_600bpwin.txt')
  a<-read.table(myfile,header=F,sep='\t')
  colnames(a)<-samps[p]
  b<-cbind(b,a)
}
rownames(b)<-b[,1]
b<-b[,-1]
write.table(b,file='allcounts_dupped_600bpwin.txt',sep='\t',quote=F)

##############################################
#Reduce to 600 basepair windows meeting DNA and RNA thresholds. 
##############################################
#Count DNA. Keep windows where at least half of the mSTARR-seq DNA samples had at least 1 read covering the window within the meth *and* sham conditions.
b<-read.delim('allcounts_dupped_600bpwin.txt')
info<-read.delim('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/coverage_report_Novaseq.txt')
info<-subset(info,libraryid!='L31250' & libraryid!='L31286') #Remove outlier RNA sample, and corresponding DNA sample
info<-subset(info,treatment=='null') #Reduce to treatment of interest ('null','dex', or 'IFN')
info<-subset(info,sampletype=='dna')
info<-droplevels(info)
dna<-b[,as.character(info$libraryid)]
dim(counts)
dim(dna)

meth<-dna[,which(info$meth=='meth')]
sham<-dna[,which(info$meth=='sham')]
methzeros<-length(which(info$meth=='meth'))-3
meth$zero_counts<-apply(meth,1,function(a) length(which(a==0)))
meth<-subset(meth,zero_counts<=methzeros)
meth<-meth[,-c(ncol(meth))]
sham$zero_counts<-apply(sham,1,function(a) length(which(a==0)))
sham<-subset(sham,zero_counts<=3)
sham<-sham[,-c(ncol(sham))]
dim(sham)
dim(meth) 
head(sham)
head(meth)

write.table(sham,file='shamdna_3ormoresamps_dupped600bpwin_nullwol31250.txt',sep='\t')
write.table(meth,file='methdna_3ormoresamps_dupped600bpwin_nullwol31250.txt',sep='\t')
dna<-merge(meth,sham,by="row.names")
rownames(dna)<-dna$Row.names
dna<-dna[,-c(1)]
head(dna)
write.table(dna,file='dnacounts_3repswith1ormore_dupped600bpwin_nullwol31250.txt',sep='\t',quote=F)  

##############################################
#FET of windows enriched for targeted sites of the genome (log2OR=3.971) baseline
##############################################
#Obtain all genome windows passing DNA filter (ie, background dataset)
cd ~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity
bedtools coverage -b windowsdnacounts_3repswith1ormore_dupped600bpwin_nullwol31250.txt -a hg38_600bp_windows.bed > genome600bp_overlapped_nullwol31250_6dna_highrepeat.txt

#Obtain target-covering genomic windows. Targets are specified in NRSA_final3_regions.bed:
bedtools coverage -b NRSA_final3_regions.bed -a hg38_600bp_windows.bed > genome600bp_overlapped_NRSA_final3_regions.txt

#In R, create dataframes of windows that contain ('b2'), or do not contain ('b3'), the target space:
setwd('~/Desktop/tunglab/capture-nrsa/mSTARR_data/counts/counts_dupped/')
b<-read.table('genome600bp_overlapped_NRSA_final3_regions.txt',header=F,sep='\t')
head(b)
b2<-subset(b,V4!=0)
dim(b)
dim(b2)
nrow(b2)/nrow(b) #Target space (0.1435303 of all 600bp genomic windows contain a target site)
b3<-subset(b,V4==0)
b2<-b2[,1:3]
head(b2)
b3<-b3[,1:3]
write.table(b2,file='genome600bp_overlapped_NRSA_final3_regions.bed',sep='\t',quote=F,row.names=F,col.names=F)
write.table(b3,file='genome600bp_NOToverlapped_NRSA_final3_regions.bed',sep='\t',quote=F,row.names=F,col.names=F)

#Now, in Unix, obtain target-covering genomic windows *passing DNA filter*, or NON-target-covering genomic windows passing DNA filter:
bedtools coverage -b windowsdnacounts_3repswith1ormore_dupped600bpwin_nullwol31250.txt -a genome600bp_overlapped_NRSA_final3_regions.bed > targetwindows_overlapped_nullwol31250_6dna_highrepeat.txt
bedtools coverage -b windowsdnacounts_3repswith1ormore_dupped600bpwin_nullwol31250.txt -a genome600bp_NOToverlapped_NRSA_final3_regions.bed > NOTtargetwindows_overlapped_nullwol31250_6dna_highrepeat.txt

#Perform Fisher's exact test in R:
nottarget<-read.table('NOTtargetwindows_overlapped_nullwol31250_6dna_highrepeat.txt',sep='\t')
target<-read.table('targetwindows_overlapped_nullwol31250_6dna_highrepeat.txt',sep='\t')
dim(nottarget)
dim(target)

a<-nrow(subset(nottarget,V4==0)) # number not-target windows not passing DNA
c<-nrow(subset(nottarget,V4!=0)) #number not-target windows passing DNA
b<-nrow(subset(target,V4==0)) #number target-containing windows not passing DNA
d<-nrow(subset(target,V4!=0)) #number target-containing windows passing DNA
values = matrix(c(a,c,b,d), nrow = 2);values #create contingency table
this<-fisher.test(values,alternative="two.sided")
this
log2(this$estimate)

##############################################
#Reduce to windows where at least 3 replicate samples produced at least 1 RNA-seq read in either the methylated condition or in the unmethylated condition. Reduce to windows with high repeatability (n = 216,091).
##############################################
#Count RNA. Keep windows where at least half of the mSTARR-seq RNA samples had at least 1 read covering the window within the meth *or* sham conditions.
b<-read.delim('allcounts_dupped_600bpwin.txt')
info<-read.delim('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/coverage_report_Novaseq.txt')
info<-subset(info,libraryid!='L31250' & libraryid!='L31286')
info<-subset(info,treatment=='null')
info<-subset(info,sampletype=='rna')
info<-droplevels(info)

counts<-b[,as.character(info$libraryid)]
meth<-counts[,which(info$meth=='meth')]
sham<-counts[,which(info$meth=='sham')]
methzeros<-length(which(info$meth=='meth'))-3
meth$zero_counts<-apply(meth,1,function(a) length(which(a==0)))
meth<-subset(meth,zero_counts<=methzeros)
meth<-meth[,-c(ncol(meth))]
head(meth)
dim(meth)
shamzeros<-length(which(info$meth=='sham'))-3
sham$zero_counts<-apply(sham,1,function(a) length(which(a==0)))
sham<-subset(sham,zero_counts<=shamzeros)
sham<-sham[,-c(ncol(sham))]
write.table(sham,file='shamrna_3ormoresamps_dupped600bpwin_nullwol31250.txt',sep='\t')
write.table(meth,file='methrna_3ormoresamps_dupped600bpwin_nullwol31250.txt.txt',sep='\t')

#Require 3 nonzero RNA samps in *either* meth or unmeth:
l<-c(rownames(meth),rownames(sham)) 
l2<-unique(l) 
counts2<-counts[l2,]
write.table(counts2,file='countsrnaw3ormoresampseither_dupped600bpwin_nullwol31250.txt',sep='\t')
counts<-counts2
dim(counts)

#Combine RNA and DNA:
dna<-read.delim('dnacounts_3repswith1ormore_dupped600bpwin_nullwol31250.txt',sep='\t')  
all<-merge(counts,dna,by='row.names',)
dim(all) 
write.table(all,file='rnadnacounts_3rna6dna_dupped600bpwin_nullwol31250.txt',sep='\t',quote=F)  

#Normalize DNA counts with limma voom
tt='null'
all<-read.table('rnadnacounts_3repswith1ormore_dupped600bpwin_nullwol31250.txt',sep='\t')  
rownames(all)<-all$Row.names;all<-all[,-c(1)];head(all)

info<-read.delim('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/coverage_report_Novaseq.txt')
info<-subset(info,libraryid!='L31250' & libraryid!='L31286')
info<-subset(info,treatment=='null') 
info<-droplevels(info)
all<-all[,info$libraryid]

dna<-all[,which(info$sampletype=='dna')]
require(limma); require(edgeR); require(statmod);library(data.table)
dge <- DGEList(counts=dna)
dge <- calcNormFactors(dge)
infodna<-subset(info,sampletype=='dna')
infodna$meth

design <- model.matrix(~c(rep(1,length(which(infodna$meth=='meth'))),rep(0,length(which(infodna$meth=='sham')))))
v <- voom(dge,design,plot=T) 
write.table(exp,file=paste0('26apr21_exp_3rna6dna_dupped600bpwin_nullwol31250.txt'),sep='\t')
fwrite(exp,file=paste0('26apr21_exptable_3rna6dna_dupped600bpwin_nullwol31250.txt')) #writes comma separated table without rownames, to be compatible with fread

#Filter for high DNA repeatability.
library(data.table)
library(qvalue)
tt='nullwol31250'
k562<-fread(paste0('26apr21_exptable_3rna6dna_dupped600bpwin_',tt,'.txt'))
head(infodna)

pairs<-as.data.frame(t(combn(which(infodna$meth=='sham'), 2)))
out<-as.data.frame(matrix(ncol=dim(pairs)[1],nrow=dim(k562)[1]))
for (i in c(1:dim(pairs)[1])) {
  rep1<-pairs[i,1]
  rep2<-pairs[i,2]
  tmp<-as.data.frame(k562[,rep1,with=FALSE] - k562[,rep2,with=FALSE])
  out[,i]<-tmp }
k562_unmeth<-out

pairs<-as.data.frame(t(combn(which(infodna$meth=='meth'), 2))) 
out<-as.data.frame(matrix(ncol=dim(pairs)[1],nrow=dim(k562)[1]))
for (i in c(1:dim(pairs)[1])) {
  rep1<-pairs[i,1]
  rep2<-pairs[i,2]
  tmp<-as.data.frame(k562[,rep1,with=FALSE] - k562[,rep2,with=FALSE])
  out[,i]<-tmp }
k562_meth<-out

# For every site, record the count difference between every pair of DNA samples 
nn<-nrow(infodna)
pairs<-as.data.frame(t(combn(c(1:nn), 2)))
out<-as.data.frame(matrix(ncol=dim(pairs)[1],nrow=dim(k562)[1]))
for (i in c(1:dim(pairs)[1])) {
  rep1<-pairs[i,1]
  rep2<-pairs[i,2]
  tmp<-as.data.frame(k562[,rep1,with=FALSE] - k562[,rep2,with=FALSE])
  out[,i]<-tmp }
k562_both<-out

x5<-quantile(as.matrix(k562_both),seq(0,1,0.025))[3] #obtain 5%
x6<-quantile(as.matrix(k562_both),seq(0,1,0.025))[39] #obtain 95%
cut1<-mean(c(x5))
cut2<-mean(c(x6))

#For each site, count how many pairs of samples are outside the central 90th percentile of distribution of ALL pairwise differences
k562_both_c1<-apply(k562_both,1,function(x) length(which(x>cut2 |x<cut1)))

#Print the proportion of sites that will be removed if >25% pairs fall outside the quantiles above (ie keeping cental 90% dist)
filt2<-dim(k562_both)[2]*0.25 #25% of number of columns
length(which(k562_both_c1>filt2)) #number of windows that are removed
length(which(k562_both_c1>filt2))/length(k562_both_c1) #proportion of windows that are removed
filt3<-which(k562_both_c1>filt2)
length(k562_both_c1)-length(which(k562_both_c1>filt2)) #number of windows that remain
# k562s
all<-read.delim(paste0('rnadnacounts_3rna6dna_dupped600bpwin_',tt,'.txt'))
rownames(all)<-all$Row.names;all<-all[,-c(1)]
allkeep<-all[c(-filt3),]
write.table(allkeep,file=paste0('rnadnacounts_3rna_6dna_highdnarepeat_dups_600bpwin_',tt,'.txt'),sep='\t')

##############################################
#Testing for regulatory activity and MD reg activity, with permutations to estimate empirical false discovery rate (FDR).
##############################################
#Create PERMUTATIONS of RNA/DNA label, within meth and within unmeth, keeping paired structure of RNA/DNA.
require(limma); require(edgeR); require(statmod); library(qvalue)

# #Prep 'info' dataframe:
zz<-'nullwol31250' #IFN or dex or nothing
tt<-'null' #IFN or dex or null
info<-read.delim('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/coverage_report_Novaseq.txt')
info<-subset(info,treatment==tt) #Keep treatment of interest
info<-subset(info,libraryid!='L31250' & libraryid!='L31286')
info2<-info

#For each set of replicate pairs that contain a rna/dna, flip coin to see if first sample in pair is dna or rna.
library(combinat)
library(gtools)
xx=c("rna","dna") #possible values for coin flip
coinflips<-permutations(n=2,r=nrow(info2)/2,v=xx,repeats.allowed=T) #n=number of values to choose from (ie 2 for heads/tails); r=number of times you choose (ie # flips)

#Permute within the meth samples
perm<-as.data.frame(matrix(nrow=nrow(info2),ncol=100)) #ncol designates number of permutations

rr<-sample(1:nrow(coinflips),100)
for (i in 1:100){
  first<-c(coinflips[rr[i],])
  second<-ifelse(first=="rna","dna","rna")
  perm[,i]<-c(first,second)
}

write.table(perm,file=paste0(zz,'_perm_sampletypewithinmethcondition28jun2021_100.txt',sep=''),sep='\t')

#Now run nested model to just get nested betas. Permute 100x by RNA (maintaining paired stucture) to identify windows with signifant regulatory activity.
#Then separate permuted results by whether permuted window has positive RNA, and take subset to have consistent # windows.
#Then run FDR estimation script. 
#Separately, run 100x by methylation status to identify windows with methylation effects. 

#Prepare voom normalized counts file
require(limma); require(edgeR); require(statmod); library(qvalue)
tt<-'nullwol31250' #IFN or dex or nullwol31250
zz<-'null' #IFN or dex or null
all<-read.table(paste0('rnadnacounts_3rna_6dna_highdnarepeat_dups_600bpwin_',tt,'.txt',sep=''),sep='\t')
info<-read.delim('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/coverage_report_Novaseq.txt')
info<-subset(info,libraryid!='L31250' & libraryid!='L31286')
info<-subset(info,treatment==zz) #Keep treatment of interest
info<-droplevels(info)
info$sampletype<-as.factor(info$sampletype)
all<-all[,info$libraryid]
info1<-info$meth
info2<-info$sampletype
design1<-model.matrix(~info2)
design1
all[is.na(all)]<-0
dge <- DGEList(counts=(all))
dge1 <- calcNormFactors(dge)
v <- voomWithQualityWeights(dge1,design=design1,plot=T)

#Run permutations:
perm<-read.table(paste0(tt,'_perm_sampletypewithinmethcondition28jun2021_100.txt',sep=''),sep='\t')
for (i in 1:100) {
  info1<-as.factor(info$meth)
  info2<-perm[,i]
  design2<-model.matrix(~info1 + info1:info2)
  design2
  fit <-lmFit(v,design2)
  fit_dna <- eBayes(fit)
  se.coef <- sqrt(fit_dna$s2.post) * fit_dna$stdev.unscaled
  fit_dna <- as.data.frame(cbind(fit_dna$p.value[,2:4],fit_dna$coefficient[,2:4],se.coef[,2:4]))
  colnames(fit_dna)<-c('meth_p','meth.rna_p','unmeth.rna_p','meth_beta','meth.rna_beta','unmeth.rna_beta','meth_se','meth.rna_se','unmeth.rna_se')
  rownames(v$E)<-rownames(all)
  rownames(fit_dna)<-rownames(all)
  write.table(fit_dna,paste0('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeat',tt,'_perm',i,'.txt',sep=''),sep='\t')
  print(i)
}

#Run model of observed (ie, non-permuted) data
info1<-as.factor(info$meth)
info2<-as.factor(info$sampletype)
design2<-model.matrix(~info1 + info1:info2)
design2
fit <-lmFit(v,design2)
fit_dna <- eBayes(fit)
se.coef <- sqrt(fit_dna$s2.post) * fit_dna$stdev.unscaled
fit_dna <- as.data.frame(cbind(fit_dna$p.value[,2:4],fit_dna$coefficient[,2:4],se.coef[,2:4]))
colnames(fit_dna)<-c('meth_p','meth.rna_p','unmeth.rna_p','meth_beta','meth.rna_beta','unmeth.rna_beta','meth_se','meth.rna_se','unmeth.rna_se')
rownames(v$E)<-rownames(all)
rownames(fit_dna)<-rownames(all)
write.table(fit_dna,paste0('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeat',tt,'.txt',sep=''),sep='\t')

#To account for multiple hypothesis testing, we calculate here the empirically derived false discovery rate.
#To do so, we'll create a dataframe *shuffled_pvals_rnainunmeth* and *shuffled_pvals_rnainmeth* with pvals from permutations (one column=1 permutation)
#We will subsample each permfile to have the same number of rows in the p-value dataframes.

#Read in p-values from permutations:
tt<-'nullwol31250'
real<-read.table(paste0('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeat',tt,'.txt',sep=''),sep='\t')
posrna<-subset(real,meth.rna_beta>0 | unmeth.rna_beta>0)

for (p in 1:100) {
  perm<-paste0('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeat',tt,'_perm',p,'.txt',sep='')
  permfile<-read.table(perm,header=T,sep="\t")
  permfile<-permfile[rownames(posrna),] #Reduce to windows that had positive RNA beta in the observed data
  permfile<-permfile[sample(nrow(permfile), nrow(posrna),replace=F), ] 
  write.table(permfile,file=paste0('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeat',tt,'_perm',p,'subsample.txt',sep=''),sep='\t')
  permfile<-permfile[,c('unmeth.rna_p','meth.rna_p')]
  
  if(p=="1")
  {
    shuffled_pvals <-data.frame(x=permfile[,c("unmeth.rna_p","meth.rna_p")])
    rownames(shuffled_pvals)=rownames(permfile)
  } else {
    shuffled_pvals <- cbind(shuffled_pvals,permfile)
  }
  print(p)
} #close p  

dim(shuffled_pvals)
shuffled_pvals_rnainunmeth<-shuffled_pvals[,seq(1,199,by=2)]
shuffled_pvals_rnainmeth<-shuffled_pvals[,seq(2,200,by=2)]

head(shuffled_pvals_rnainunmeth)
head(shuffled_pvals_rnainmeth)

real<-real[rownames(posrna),]
length(which(real$meth.rna_beta>0))
length(which(real$unmeth.rna_beta>0))

#Correct for Multiple testing and filter results to report
#perm.fdr is the function we use to calculate the empirical FDR by comparing the observed p-values to the permuted p-values.
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

#perm.fdr=function(input_df,permuted_df,input_column_name_observed_pvalues,output_column_name_prefix)

#Apply perm.fdr to calculate FDR for effect of sampletype (RNA or DNA) in the methylated samples:
res_full=perm.fdr(data.frame(real),shuffled_pvals_rnainmeth,"meth.rna_p","RNAinmeth")
res_full=res_full[order(rownames(res_full)),]
head(res_full)

#Apply perm.fdr to calculate FDR for effect of sampletype (RNA or DNA) in the unmethylated (ie, sham) samples:
res_full=perm.fdr(data.frame(res_full),shuffled_pvals_rnainunmeth,"unmeth.rna_p","RNAinunmeth")
res_full=res_full[order(rownames(res_full)),]

head(res_full)
write.table(res_full,paste0('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeat',tt,'_FDR.txt',sep=''),sep="\t")

#Now create dataframe of permuted methylation condition. Create permutation file, permuting meth/unmeth within dna samples and matching to rna samples given the paired structure of the data.
#Prep info for nrsa:
zz<-'null' #or IFN or null or dex
tt<-'nullwol31250' #IFN or dex or nullwol31250
info<-read.delim('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/coverage_report_Novaseq.txt')
info<-subset(info,treatment==zz) #Keep treatment of interest
info<-subset(info,libraryid!='L31250' & libraryid!='L31286')
info<-droplevels(info)
info$sampletype<-as.factor(info$sampletype)
infodna<-subset(info,sampletype=='dna')
cols<-as.data.frame(matrix(nrow=nrow(info),ncol=100)) #designate # perms

for (i in 1:100) {
  random<-sample(infodna$meth,nrow(infodna),replace=F)
  cols[,i]<-c(random,random)
}
write.table(cols,file=paste0('1july2021permutedmeth100',tt,'.txt'),sep='\t')

#Run permutations of meth
require(limma); require(edgeR); require(statmod); library(qvalue)

zz<-'null' #or IFN or null or dex
tt<-'nullwol31250' #IFN or dex or nullwol31250
all<-read.table(paste0('rnadnacounts_3rna_6dna_highdnarepeat_dups_600bpwin_',tt,'.txt',sep=''),sep='\t')
info<-read.delim('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/coverage_report_Novaseq.txt')
info<-subset(info,treatment==zz) #Keep treatment of interest
info<-subset(info,libraryid!='L31250' & libraryid!='L31286')
info<-droplevels(info)
info$sampletype<-as.factor(info$sampletype)
all<-all[,info$libraryid]
info2<-info$sampletype
design1<-model.matrix(~info2)
design1
all[is.na(all)]<-0
dge <- DGEList(counts=(all))
dge1 <- calcNormFactors(dge)
v <- voomWithQualityWeights(dge1,design=design1,plot=T)

cols<-read.table(paste0('1july2021permutedmeth100',tt,'.txt'),sep='\t')
for (i in 1:100) {
  info1<-as.factor(cols[,i])
  info2<-info$sampletype
  design2<-model.matrix(~info1 + info1:info2)
  design2
  fit <-lmFit(v,design2)
  fit_dna <- eBayes(fit)
  se.coef <- sqrt(fit_dna$s2.post) * fit_dna$stdev.unscaled
  fit_dna <- as.data.frame(cbind(fit_dna$p.value[,2:4],fit_dna$coefficient[,2:4],se.coef[,2:4]))
  colnames(fit_dna)<-c('meth_p','meth.rna_p','unmeth.rna_p','meth_beta','meth.rna_beta','unmeth.rna_beta','meth_se','meth.rna_se','unmeth.rna_se')
  rownames(v$E)<-rownames(all)
  rownames(fit_dna)<-rownames(all)
  write.table(fit_dna,paste0('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeat',tt,'_permmeth',i,'.txt',sep=''),sep='\t')
  print(i)
}

#Now apply perm.fdr code to permutations of meth in nested model, subsetting to the window subsamples that 
#had been used for generating the empirical null of regulatory windows.
tt<-'nullwol31250' #IFN or dex or nullwol31250
df<-read.table(paste0('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeat',tt,'_permmeth1.txt'),sep='\t')

for (p in 1:100) {
  perm<-paste0('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeat',tt,'_permmeth',p,'.txt',sep='')
  permfile<-read.table(perm,header=T,sep="\t")
  permfile$int_beta<-permfile$meth.rna_beta-permfile$unmeth.rna_beta
  permfile$int_se<-sqrt(permfile$meth.rna_se^2 + permfile$unmeth.rna_se^2)
  permfile$int_t=permfile$int_beta/permfile$int_se
  permfile$int_p=2*pt(-abs(permfile$int_t),df=20) 
  rownames(permfile)<-rownames(df)
  window<-paste0('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeat',tt,'_perm',p,'subsample.txt',sep='')
  windowfile<-read.table(window,header=T,sep="\t")
  
  permfile<-permfile[rownames(windowfile),]
  
  if(p=="1")
  {
    shuffled_pvals_rnaxmeth <-data.frame(x=permfile[,"int_p"])
    rownames(shuffled_pvals_rnaxmeth)=rownames(permfile)
  } else {
    shuffled_pvals_rnaxmeth <- cbind(shuffled_pvals_rnaxmeth,x=permfile[,"int_p"])
  }
  print(p)
} #close p

#Read in results of model performed on observed data
real<-read.table(paste0('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeat',tt,'_FDR.txt'),header=T,sep="\t")

#Manually calculate p-values for methylation*sampletype interaction term.
real$int_beta<-real$meth.rna_beta-real$unmeth.rna_beta
real$int_se<-sqrt(real$meth.rna_se^2 + real$unmeth.rna_se^2)
real$int_t=real$int_beta/real$int_se
real$int_p=2*pt(-abs(real$int_t),df=20) 

#Apply perm.fdr to calculate FDR for interaction term methylation*sampletyp:
#perm.fdr=function(input_df,permuted_df,input_column_name_observed_pvalues,output_column_name_prefix)
res_full=perm.fdr(data.frame(real),shuffled_pvals_rnaxmeth,"int_p","int")
res_full=res_full[order(rownames(res_full)),]
write.table(res_full,paste0("1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeat",tt,"_FDR2.txt"),sep="\t")

##############################################
#Reg activity enrichment in ENCODE's 15 ChromHMM genomic region types
##############################################
#In unix terminal. For each ChromHMM region type, create file that shows the overlap of ChromHMM region with tested windows
#option -f requires that at least half of a is covered by b. -u says to only report a once if it meets overlap criteria at least once.
d=nullwol31250
for t in 10_Txn_Elongation 11_Weak_Txn 12_Repressed 13_Heterochromlo 14_RepetitiveCNV 15_RepetitiveCNV 1_Active_Promoter 2_Weak_Promoter 3_Poised_Promoter 4_Strong_Enhancer 5_Strong_Enhancer 6_Weak_Enhancer 7_Weak_Enhancer 8_Insulator 9_Txn_Transition;do
bedtools intersect -a ~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/windows1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeat${d}.txt -b ~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/references/K562_Hg38_${t}.bed -f 0.5 -u > testedwindows600_METH_3rna_6dna_nested${d}overlappedby_${t}.txt 
done

#In R, perform Fisher's exact tests. 'All' is background set of all tested windows; 'reg' is the set of regulatory windows
setwd('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/')
tt<-'nullwol31250'
real<-read.table(paste0('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeat',tt,'.txt',sep=''),sep='\t')
result0<-read.table(paste0('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeat',tt,'_FDR.txt'),sep='\t')

cc=0.01
result<-subset(result0,(unmeth.rna_beta>0 & fdr_RNAinunmeth<cc)|(meth.rna_beta>0 & fdr_RNAinmeth<cc))
final<-as.data.frame(matrix(nrow=15,ncol=4))
colnames(final)<-c('OR','lower','upper','p')
rownames(final)<-c("10_Txn_Elongation","11_Weak_Txn","12_Repressed","13_Heterochromlo","14_RepetitiveCNV","15_RepetitiveCNV","1_Active_Promoter","2_Weak_Promoter","3_Poised_Promoter","4_Strong_Enhancer","5_Strong_Enhancer","6_Weak_Enhancer","7_Weak_Enhancer","8_Insulator","9_Txn_Transition")
final$color<-rep('gray',15)
final[c(10,11,12,13),5]<-'darkblue'  
regionlist<-c("10_Txn_Elongation","11_Weak_Txn","12_Repressed","13_Heterochromlo","14_RepetitiveCNV","15_RepetitiveCNV","1_Active_Promoter","2_Weak_Promoter","3_Poised_Promoter","4_Strong_Enhancer","5_Strong_Enhancer","6_Weak_Enhancer","7_Weak_Enhancer","8_Insulator","9_Txn_Transition")
for (t in 1:15) {
  region<-regionlist[t]
  all<-read.table(paste0('testedwindows600_METH_3rna_6dna_nested',tt,'overlappedby_',region,'.txt'))
  reg_in_region<-all[all$V4 %in% rownames(result),]
  notreg_in_region<-all[!all$V4 %in% rownames(result),]
  d<-nrow(reg_in_region);d
  c<-nrow(result)-d;c
  b<-nrow(notreg_in_region);b
  a<-nrow(real)-b-c-d #total windows tested - b-c-d
  values = matrix(c(a,b,c,d), nrow = 2) #create contingency table
  values
  this<-fisher.test(values,alternative="two.sided")
  this
  
  final[t,1]<-this$estimate
  final[t,2]<-this$conf.int[1]
  final[t,3]<-this$conf.int[2]
  final[t,4]<-this$p.value
  print(t)
} #close t

final

#Plot the results as barplot, ordered by OR
ccol<-c(rep('#f3ebf0',10),'#a1b7c4','#f3ebf0',rep('#a1b7c4',3))
cc=0.01
par(mfrow=c(1,1),mar=c(12,4,1,1))
final$log2or<-log2(final$OR)
final2<-final[order(final$log2or),]
barplot(final2$log2or,col=as.character(ccol),names.arg=rownames(final2),las=2,ylab="log2(OR)",ylim=c(-3,3.7),main="")
abline(h=1,col='#753b19',lty='dashed');abline(h=2,col='#753b19',lty='dashed');abline(h=3,col='#753b19',lty='dashed')

##############################################
#Perform Fisher's Exact Test for overlap of regulatory regions between this study's baseline dataset and Lea et al 2018 dataset, 
#and look at correlation of the two-studies' effects of methylation on regulatory activity
##############################################
#Apply the pipeline above to the Lea et al 2018 dataset (obtained from NCBI’s SRA accession number SRP120556) 
#to trim and map reads to hg38, create counts file, and perform nested model and estimate empirical FDR.

#Create dataframe of 600 bp windows tested in both datasets
setwd('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity')
n<-read.table('windows1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwol31250.txt')
k<-read.table('windows18nov2020_nested_model_METH_3rna_6dna_dup_highdnarepeatk562.txt')
rownames(n)<-n$V4
rownames(k)<-k$V4
nk<-merge(n,k,by='row.names')
dim(n)
dim(k)
dim(nk)
rownames(nk)<-nk$Row.names

#Create dataframe of regulatory regions identified in this study that were tested in both studies
cc=0.01
nreg<-read.table('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwol31250_FDR2.txt')
nreg<-nreg[rownames(nreg) %in% rownames(nk),]
nreg<-subset(nreg,(unmeth.rna_beta>0 & fdr_RNAinunmeth<cc) | (meth.rna_beta>0 & fdr_RNAinmeth<cc))
dim(nreg)

#Create dataframe of regulatory regions identified in the Lea et al 2018 study that were tested in both studies
kreg<-read.table('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatk562_FDR2.txt')
kreg<-kreg[rownames(kreg) %in% rownames(nk),]
kreg<-subset(kreg,(unmeth.rna_beta>0 & fdr_RNAinunmeth<cc) | (meth.rna_beta>0 & fdr_RNAinmeth<cc))
dim(kreg)

#Create dataframes for regulatory regions shared in both studies, specific to one study, are non-regulatory in both studies 
nregkreg<-merge(nreg,kreg,by='row.names')
knot<-nk[!rownames(nk) %in% rownames(kreg),]
nrow(nregkreg)
nnot<-nk[!rownames(nk) %in% rownames(nreg),]
kregnnot<-kreg[!rownames(kreg) %in% rownames(nreg),]
nregknot<-nreg[!rownames(nreg) %in% rownames(kreg),]
nnotknot<-merge(nnot,knot,by='row.names')

#Perform Fisher's Exact Test to test for overlap of regulatory regions between the two studies
aa<-nrow(nnotknot)
bb<-nrow(kregnnot)
cc<-nrow(nregknot)
dd<-nrow(nregkreg)
values = matrix(c(aa,bb,cc,dd), nrow = 2)
values<-as.data.frame(values)
rownames(values)<-c('elife_not_reg','elife_reg')
colnames(values)<-c('NRSA_not_reg','NRSA_reg')
values
this<-fisher.test(values,alternative="two.sided");this
log2(this$estimate)

#Perform FET for overlap of methylation-dependent regulatory windows
cc=0.01
head(nregkreg)
nnotmdknotmd<-subset(nregkreg, fdr_int.x>cc & fdr_int.y>cc)
kmdnnot<-subset(nregkreg, fdr_int.x>cc & fdr_int.y<cc)
nmdknot<-subset(nregkreg, fdr_int.x<cc & fdr_int.y>cc)
nmdkmd<-subset(nregkreg, fdr_int.x<cc & fdr_int.y<cc)

aa<-nrow(nnotmdknotmd)
bb<-nrow(kmdnnot)
cc<-nrow(nmdknot)
dd<-nrow(nmdkmd)
values = matrix(c(aa,bb,cc,dd), nrow = 2)
values<-as.data.frame(values)
rownames(values)<-c('elife_not_MD','elife_MD')
colnames(values)<-c('NRSA_not_MD','NRSA_MD')
values
this<-fisher.test(values,alternative="two.sided");this

plot(nregkreg$int_beta.x,nregkreg$int_beta.y,pch=19,cex.lab=1.4,cex.axis=1.4,xlab='Int beta NRSA',ylab='Int beta eLife',xlim=c(-12.1,5.5),ylim=c(-12.1,5.5))
l<-lm(nregkreg$int_beta.y~nregkreg$int_beta.x)
abline(h=0,lty='dashed')
abline(v=0,lty='dashed')
abline(l)
cor(nregkreg$int_beta.x,nregkreg$int_beta.y)
xx<-summary(l)
xx$r.squared
xx$coefficients

##############################################
#Bisulfite sequencing full pipeline
##############################################
#Create directory for trimmed reads
mkdir trimmed

#Submit the job below to slurm:
for f in `cat ids.txt`; do cat methratioredo2.sh | sed -e s/FILE/$f/g > methratioredo2_$f.sh; sbatch methratioredo2_$f.sh; done #methratioredo2.sh created28jan2021

#!/bin/bash
#SBATCH -n 1 # Number of cores
#SBATCH -t 0-96:00 # Runtime in D-HH:MM
#SBATCH --mem=32G 
#SBATCH -o redo_FILE.out # File to which STDOUT will be written
#SBATCH -e redo_FILE.err # File to which STDERR will be written
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=youremail@gmail.com # Email to which notifications will be sent

module load TrimGalore
module load cutadapt
module load BSMAP
module load python/2.7.6-fasrc01
module load samtools/0.1.19-fasrc01

#Trim adapters from reads
trim_galore -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o trimmed --gzip --rrbs --length 15 --stringency 4 --paired FILE_R1.fastq.gz FILE_R2.fastq.gz

#Map reads to reference
#Specify path to genome of interest. Here, a combined human + mstarr plasmid + Lambda Phage genome so reads could be mapped and bisulfite conversion rate could be estimated in one step.
path_genome=hg38_wphage_wplasmid_wmSTARR_primerR2.fa #this was methratioredo2.sh
#Path to desired sam output file
path_sam=FILE_redo2.sam
#Map. Option -r 0 reports unique hits only. -v sets allowed # of mismatches. If between 0 and 1, it is interpreted as a % of the read length (i.e. here 10% of read length is allowed to be mismatched)
bsmap -a trimmed/FILE_R1_val_1.fq.gz -b trimmed/FILE_R2_val_2.fq.gz -p 8 -d $path_genome -o $path_sam -v 0.1 -r 0

#Call methratios
#--context=CG will only return mratio calls at CpG dinucleotides
#--combine-CpG will combine mratio call information from both strands
path_out=FILE.methratio_redo2.txt
python methratio.py -o $path_out -d $path_genome --combine-CpG --context=CG $path_sam
#end of methratioredo2.sh script submitted to slurm


#Pull plasmid CpGs, which will have the highest coverage and will therefore be best sites for estimating overall methylation levels
for f in `cat ids.txt`;do
grep -w plasmid_w_adapter ${f}.methratio_redo2.txt > ${f}.methratio_redo2_plasmid.txt
cat methratio_header.txt ${f}.methratio_redo2_plasmid.txt > ${f}.methratio_redo2_plasmid2.txt
done

#In R, compare methylation levels between methylated and unmethylated samples (baseline K562 samples):
setwd('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity')
library(beeswarm)
info<-read.table('libinfobs.txt',header=T)
info$meth_condition<-as.factor(info$meth_condition)
summary(info$meth_condition)
ids<-c("L31374","L31375","L31376","L31377","L31378","L31379","L31380","L31381","L31382","L31383",
       "L31384","L31385","L31386","L31387","L31388","L31389","L31390","L31391","L31392","L31393","L31394",
       "L31395","L31396","L31397","L31398","L31399","L31400","L31401","L31402","L31404","L31405",
       "L31406","L31407","L31411","L31434")
#Leave out L31403, L31408, and L31409 because, due to low sequencing coverage they are each missing some of the 18 plasmid sites

rownames(info)<-info$sample
out<-as.data.frame(matrix(nrow=18,ncol=length(ids)))
colnames(out)<-ids

for (p in 1:length(ids)) {
  print(p)
  filename<-paste0(ids[p],'.methratio_redo2_plasmid2.txt')
  myfile<-read.table(filename,header=T)
  myfile$pos2<-paste0(myfile$chr,'_',myfile$pos)
  myfile$meth_level<-myfile$C_count/myfile$eff_CT_count
  out[,p]<-myfile$meth_level
} #close p
rownames(out)<-myfile$pos2

#Get average and median methylation level for methylated and unmethylated samples at a representative site for plotting
out2<-out['plasmid_w_adapter_2294',]
info2<-info[ids,]
info2meth<-subset(info2,meth_condition=='meth')
info2unmeth<-subset(info2,meth_condition=='sham')
out2meth<-out2[info2meth$sample]
out2unmeth<-out2[info2unmeth$sample]

#Summarize methylation levels of methylated samples
out3meth<-t(out2meth[1:17])
summary(out3meth) #only includes post-transformation
#mean: 0.8848
#median: 0.8965

#Summarize methylation levels of unmethylated samples
out3unmeth<-t(out2unmeth[1:15])
summary(out3unmeth) #only includes post-transformation
#mean: 0.066006
#median: 0.010438

t.test(out3unmeth,out3meth,paired=F)

#Create PCA plot of overall mSTARR-Seq reads for the Dex-treated RNA samples, color-coded by methylation status
setwd('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity')
library(beeswarm)
all<-read.table(paste0('rnadnacounts_3rna_6dna_highdnarepeat_dups_600bpwin_dex.txt',sep=''),sep='\t')
info<-read.delim('coverage_report_Novaseq.txt')
info<-subset(info,treatment=='dex') #Keep treatment of interest
info<-subset(info,sampletype=='rna')
info<-droplevels(info)
info$sampletype<-as.factor(info$sampletype)
all<-all[,info$libraryid]
ids<-info$libraryid

par(mfrow=c(1,1))
pca<-prcomp(t(cov(all)),center=T,scale=T)
plot(pca$x[,1],pca$x[,2],pch=ifelse(info$sampletype=='rna',19,21),bty='n',col=as.factor(info$meth),main='PCA of cov matrix')
text(pca$x[,1]+0.5,pca$x[,2],labels=ids,cex=0.7,cex.axis=1.4,cex.lab=1.4)

##############################################
#FET for overlap of regulatory regions baseline dataset from this study vs Johnson et al 2018 dataset ((at t=0 unstimulated A549s))
##############################################
#Create bedfile of tested windows that were not regulatory 
setwd('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity')
tested<-read.table('windows1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwol31250.txt',header=F)
reg<-read.table('regfdr0.01_1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwol31250.bed',header=F)
rownames(tested)<-tested$V4
rownames(reg)<-reg$V4
notreg<-subset(tested, !rownames(tested) %in% rownames(reg))
head(notreg)
write.table(notreg,file='notregfdr0.01_1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwol31250.bed',sep='\t',row.names=F,col.names=F,quote=F)

#Creat files of overlap of regulatory activity
cd ~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity
bedtools intersect -a regfdr0.01_1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwol31250.bed -b GDJohnson_NatCom_peaks_t00_union.bed -wo > regfdr0.01nullwol31250_natcompeakst00.txt
bedtools intersect -a windows1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwol31250.txt -b GDJohnson_NatCom_peaks_t00_union.bed -wo > testednullwol31250_natcompeakst00.txt
bedtools intersect -a notregfdr0.01_1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwol31250.bed -b GDJohnson_NatCom_peaks_t00_union.bed -wo > notregfdr0.01nullwol31250_natcompeakst00.txt

#Perform Fisher's Exact Test for overlap of regulatory space between the two datasets
setwd('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity')
tested<-read.table('windows1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwol31250.txt',header=F)
reg<-read.table('regfdr0.01_1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwol31250.bed',header=F)
rownames(tested)<-tested$V4
rownames(reg)<-reg$V4
notreg<-subset(tested, !rownames(tested) %in% rownames(reg))
nrow(notreg)+nrow(reg)
nrow(tested)

#This study's tested space null (bp):
nrow(tested)*600
129654600

#This study's regulatory space (bp):
regthisstudy<-nrow(reg)*600
regthisstudy
2232600

#My reg and Johnson reg
regboth<-read.table('~/Desktop/tunglab/capture-nrsa/johnson_2018/regfdr0.01nullwol31250_natcompeakst00.txt')
regboth<-sum(regboth$V9)
102976

#My reg and Johnson not reg
regthisstudy_notregjohnson<-regthisstudy-regboth
regthisstudy_notregjohnson
2129624

#Johnson reg space in this study's tested
johnsonreg<-read.table('~/Desktop/tunglab/capture-nrsa/johnson_2018/testednullwol31250_natcompeakst00.txt')
johnsonreg<-sum(johnsonreg$V9)
johnsonreg
1219688

#Johnson reg, not reg this study
johnsonreg_notregthisstudy<-johnsonreg-regboth
johnsonreg_notregthisstudy
1116712

#Not reg this study
notregthisstudy<-nrow(notreg)*600
notregthisstudy
127422000

#Not reg this study and Johnson not reg
notregthisstudy_johnstonnotreg<-notregthisstudy - johnsonreg_notregthisstudy
notregthisstudy_johnstonnotreg

values = matrix(c(notregthisstudy_johnstonnotreg,johnsonreg_notregthisstudy,regthisstudy_notregjohnson,regboth), nrow = 2);values
values
this<-fisher.test(values,alternative="two.sided");this
log2(this$estimate)

##############################################
#Correlation of MD effects in baseline vs HepG2
##############################################
setwd('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity')

#Create a dataframe of windows tested for MD regulatory activity in both the baseline K562 and HepG2 datasets
ktest<-read.table('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwol31250.txt',sep='\t')
htest<-read.table('18nov2020_nested_model_METH_3rna_6dna_dup_highdnarepeathepg2.txt',sep='\t')
hktest<-merge(htest,ktest,by='row.names')
rownames(hktest)<-hktest$Row.names

#Pull MD regulatory results from windows tested in both the baseline K562 and HepG2 datasets
h<-read.table('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeathepg2_FDR2.txt')
k<-read.table('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwol31250_FDR2.txt')
h <- h[rownames(h) %in% rownames(hktest),]
k <- k[rownames(k) %in% rownames(hktest),]
cc=0.01
hreg<-subset(h,(unmeth.rna_beta>0 & fdr_RNAinunmeth<cc) | (meth.rna_beta>0 & fdr_RNAinmeth<cc))
kreg<-subset(k,(unmeth.rna_beta>0 & fdr_RNAinunmeth<cc) | (meth.rna_beta>0 & fdr_RNAinmeth<cc))
hkreg<-merge(hreg,kreg,by='row.names')

#Measure correlation of MD regulatory effects in K562 and HepG2 regulatory windows
l<-lm(hkreg$int_beta.y~hkreg$int_beta.x)
b<-summary(l)
b
b$coefficients

##############################################
#HOMER for regulatory windows showing higher regulatory activity in the methylated relative to unmethylated condition (baseline K562s)
##############################################
setwd('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity')
realf<-read.table('1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwol31250_FDR2.txt',sep='\t')
real2<-subset(realf, (meth.rna_beta>0 & fdr_RNAinmeth<0.01) | (unmeth.rna_beta>0 & fdr_RNAinunmeth<0.01))
real3<-subset(real2,fdr_int<0.01)
real4<-subset(real3,int_beta>0) #Unexpected direction-- regulatory windows with more regulatory activity in the *methylated* relative to unmethylated condition
dim(real4) #24 windows

#Create text files of background and test sets to read into HOMER for TFBS motif enrichment analysis
real2[,1]<-rownames(real2)
real2<-real2[,1:2]
write.table(real2,file='enhancers_3721.txt',sep='\t',quote = F,row.names=F,col.names=F)

real4[,1]<-rownames(real4)
real4<-real4[,1:2]
write.table(real4,file='mdenhancers_activeinmeth24.txt',sep='\t',quote = F,row.names=F,col.names=F)

#Run HOMER from terminal to perform TFBS motif enrichment
cd ~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity
genome=hg38.fa
findMotifsGenome.pl mdenhancers_activeinmeth24.txt $genome activeinmeth -size given -bg enhancers_3721.txt -mknown known.motifs -nomotif
mv knownResults.txt knownResultsactiveinmeth.txt

#In R, perform FET for HOMER output ('moreactiveinmeth')
h<-read.table('knownResultsactiveinmeth.txt',sep='\t',header=T)
h$or<-1;h$loweror<-1;h$upperor<-1;h$p<-1
head(h)
for (p in 1:nrow(h)){
  out<-h[p,]
  a <- out$a
  b <- out$b
  c <- out$c
  d <- out$d
  values = matrix(c(a,b,c,d), nrow = 2);values
  this<-fisher.test(values,alternative="two.sided");this
  h[p,'p']=this$p.value
  h[p,'or']=this$estimate
  h[p,'loweror']=this$conf.int[1]
  h[p,'upperor']=this$conf.int[2]
} #close p
write.table(h,file='knownResultsactiveinmeth_fet.txt',sep='\t',quote=F)

