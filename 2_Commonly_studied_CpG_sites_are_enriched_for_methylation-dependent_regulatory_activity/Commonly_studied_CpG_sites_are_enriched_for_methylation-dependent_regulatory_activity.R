##############################################
#Test for enrichment of windows EPIC vs NR3C1 vs control loci passing DNA filter, that also pass RNA filter
##############################################
#For each 600bp window passing DNA filter, mark the number of target CpG sites overlapping it (separated by target set- EPIC, MSP1 RRBS digest, NR3C1, or control)
cd ~/Desktop/github_mstarr/2_Commonly_studied_CpG_sites_are_enriched_for_methylation-dependent_regulatory_activity
bedtools coverage -a dnacounts_3repswith1ormore_dupped600bpwin_nullwol31250.bed -b NRSA_final3_regions_cg.bed > epic_overlapped_dnacounts_dna6reps1ormore_dupped600bpwin_nullwol31250_2dec22.txt
bedtools coverage -a dnacounts_3repswith1ormore_dupped600bpwin_nullwol31250.bed -b NRSA_final3_regions_control.bed > control_overlapped_dnacounts_dna6reps1ormore_dupped600bpwin_nullwol31250_2dec22.txt
bedtools coverage -a dnacounts_3repswith1ormore_dupped600bpwin_nullwol31250.bed -b NRSA_final3_regions_insilico.bed > msp1_overlapped_dnacounts_dna6reps1ormore_dupped600bpwin_nullwol31250_2dec22.txt
bedtools coverage -a dnacounts_3repswith1ormore_dupped600bpwin_nullwol31250.bed -b NRSA_final3_regions_NR3C1.bed > NR3C1_overlapped_dnacounts_dna6reps1ormore_dupped600bpwin_nullwol31250_2dec22.txt

#In R, perform Fisher's Exact Tests
setwd('~/Desktop/github_mstarr/2_Commonly_studied_CpG_sites_are_enriched_for_methylation-dependent_regulatory_activity')
regions<-c('epic','msp1','NR3C1','control')
result<-as.data.frame(matrix(nrow=length(regions),ncol=8))
my<-read.delim('dnacounts_3repswith1ormore_dupped600bpwin_nullwol31250.bed',sep='\t')

#Read in dataframe ('m') that specifies all windows that passed the RNA filter
m<-read.delim('../1_mSTARR-seq_captures_enhancer_activity/1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwol31250.txt',sep='\t')

#Perform FETs
rownames(result)<-regions;colnames(result)<-c('OR','OR_lower_95CI','OR_upper_95CI','p-value','bg_notpass','focal_notpass','bgpass','focalpass')
for (i in 1:length(regions)) {
  f=regions[i]
  d<-'dnacounts_dna6reps1ormore_dupped600bpwin_nullwol31250_2dec22'
  all<-read.delim(paste0(f,'_overlapped_',d,'.txt'),header=F)
  all<-subset(all,V5>0)
  bg<-read.delim('control_overlapped_dnacounts_dna6reps1ormore_dupped600bpwin_nullwol31250.txt',header=F)
  bg<-subset(bg,V5>0)
  rownames(all)<-paste(all$V1,all$V2,all$V3,sep='_')
  rownames(bg)<-paste(bg$V1,bg$V2,bg$V3,sep='_')
  allnotsurvived<-all[!rownames(all) %in% rownames(m), ]
  allsurvived<-all[rownames(all) %in% rownames(m), ]
  bgnotsurvived<-bg[!rownames(bg) %in% rownames(m), ]
  bgsurvived<-bg[rownames(bg) %in% rownames(m), ]
  aa<-nrow(bgnotsurvived)
  bb<-nrow(allnotsurvived)
  cc<-nrow(bgsurvived)
  dd<-nrow(allsurvived)
  values = matrix(c(aa,bb,cc,dd), nrow = 2);values
  this<-fisher.test(values,alternative="two.sided");this
  # rownames(values)<-c('control',f);colnames(values)<-c('notregulatory','regulatory')
  # write.table(as.data.frame(values),file=paste0('fet_contingencytable_regactivity_',f,'.txt'),sep='\t',row.names=T,col.names=T,quote=F)
  result[i,1]<-this$estimate
  result[i,2]<-this$conf.int[1]
  result[i,3]<-this$conf.int[2]
  result[i,4]<-this$p.value
  result[i,5]<-aa
  result[i,6]<-bb
  result[i,7]<-cc
  result[i,8]<-dd
  print(i)
}
result$log2OR<-log2(result$OR)
result$log2OR_lower_95CI<-log2(result$OR_lower_95CI)
result$log2OR_upper_95CI<-log2(result$OR_upper_95CI)
result
##############################################
#EPIC, NR3C1, RRBS vs control passing RNA filter showing reg activity
##############################################
#Create file that specifies, for each target type, the windows it was covered by in tested dataset (baseline K562s)
#Target types: cg=EPIC, insilico=RRBS, control=control, NR3C1=NR3C1
cd ~/Desktop/github_mstarr/2_Commonly_studied_CpG_sites_are_enriched_for_methylation-dependent_regulatory_activity
c=nullwol31250
for f in control cg insilico NR3C1; do
bedtools coverage -b ~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/windows1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwol31250.txt -a NRSA_final3_regions_${f}.bed > ${f}_overlapped_${c}.txt;done

#Perform FETs
setwd('~/Desktop/github_mstarr/2_Commonly_studied_CpG_sites_are_enriched_for_methylation-dependent_regulatory_activity')
cc=0.01
regions<-c('cg','insilico','NR3C1','control')
my<-read.delim('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwol31250_FDR2.txt',sep='\t')
m<-subset(my,(unmeth.rna_beta>0 & fdr_RNAinunmeth<cc) | (meth.rna_beta>0 & fdr_RNAinmeth<cc))

result<-as.data.frame(matrix(nrow=4,ncol=8))
rownames(result)<-regions;colnames(result)<-c('OR','OR_lower_95CI','OR_upper_95CI','p-value','a','b','c','d')

for (i in 1:4) {
  f=regions[i]
  d<-'nullwol31250'
  all<-read.delim(paste0(f,'_overlapped_',d,'.txt'),header=F)
  all<-subset(all,V5>0) #Subsets to tested windows that contained region CpG
  bg<-read.delim(paste0('control_overlapped_',d,'.txt'),header=F)
  bg<-subset(bg,V5>0)
  rownames(all)<-paste(all$V1,all$V2,all$V3,sep='_')
  rownames(bg)<-paste(bg$V1,bg$V2,bg$V3,sep='_')
  allnotreg<-all[!rownames(all) %in% rownames(m), ]
  allreg<-all[rownames(all) %in% rownames(m), ]
  bgnotreg<-bg[!rownames(bg) %in% rownames(m), ]
  bgreg<-bg[rownames(bg) %in% rownames(m), ]
  aa<-nrow(bgnotreg)
  bb<-nrow(allnotreg)
  cc<-nrow(bgreg)
  dd<-nrow(allreg)
  values = matrix(c(aa,bb,cc,dd), nrow = 2);values
  this<-fisher.test(values,alternative="two.sided");this
  result[i,1]<-this$estimate
  result[i,2]<-this$conf.int[1]
  result[i,3]<-this$conf.int[2]
  result[i,4]<-this$p.value
  result[i,5]<-aa
  result[i,6]<-bb
  result[i,7]<-cc
  result[i,8]<-dd
  print(i)
}
result$log2OR<-log2(result$OR)
result$log2OR_lower_95CI<-log2(result$OR_lower_95CI)
result$log2OR_upper_95CI<-log2(result$OR_upper_95CI)
result
write.table(result,file=paste0("fet_captureregactivity_",d,".txt"),sep='\t')

##############################################
#EPIC, NR3C1, RRBS vs control regulatory windows showing MD activity
##############################################
#Create file that specifies, for each target type, the windows it was covered by in regulatory windows dataset 'nullwol31250'
cd ~/Desktop/github_mstarr/2_Commonly_studied_CpG_sites_are_enriched_for_methylation-dependent_regulatory_activity
c=regfdr0.01nullwol31250
for f in control cg insilico NR3C1; do
bedtools coverage -b ~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/regfdr0.01_1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwol31250.bed -a NRSA_final3_regions_${f}.bed > ${f}_overlapped_${c}.txt;done

#In R, create dataframe of windows with MD regulatory activity
setwd('~/Desktop/github_mstarr/2_Commonly_studied_CpG_sites_are_enriched_for_methylation-dependent_regulatory_activity')
regions<-c('cg','insilico','NR3C1','control')
result<-as.data.frame(matrix(nrow=4,ncol=8))
rownames(result)<-regions;colnames(result)<-c('OR','OR_lower_95CI','OR_upper_95CI','p-value','a','b','c','d')
my<-read.delim('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/1july2021_nested_model_METH_3rna_6dna_dup_highdnarepeatnullwol31250_FDR2.txt',sep='\t')
cc=0.01
m<-subset(my,(unmeth.rna_beta>0 & fdr_RNAinunmeth<cc) | (meth.rna_beta>0 & fdr_RNAinmeth<cc))
m<-subset(m,fdr_int<cc)

#Perform FETs
for (i in 1:4) {
  f=regions[i]
  d<-'regfdr0.01nullwol31250'
  all<-read.delim(paste0(f,'_overlapped_',d,'.txt'),header=F)
  all<-subset(all,V5>0)
  bg<-read.delim(paste0('control_overlapped_',d,'.txt'),header=F)
  bg<-subset(bg,V5>0)
  rownames(all)<-paste(all$V1,all$V2,all$V3,sep='_')
  rownames(bg)<-paste(bg$V1,bg$V2,bg$V3,sep='_')
  allnotreg<-all[!rownames(all) %in% rownames(m), ]
  allreg<-all[rownames(all) %in% rownames(m), ]
  bgnotreg<-bg[!rownames(bg) %in% rownames(m), ]
  bgreg<-bg[rownames(bg) %in% rownames(m), ]
  aa<-nrow(bgnotreg)
  bb<-nrow(allnotreg)
  cc<-nrow(bgreg)
  dd<-nrow(allreg)
  values = matrix(c(aa,bb,cc,dd), nrow = 2);values
  this<-fisher.test(values,alternative="two.sided");this
  result[i,1]<-this$estimate
  result[i,2]<-this$conf.int[1]
  result[i,3]<-this$conf.int[2]
  result[i,4]<-this$p.value
  result[i,5]<-aa
  result[i,6]<-bb
  result[i,7]<-cc
  result[i,8]<-dd
}
result$log2OR<-log2(result$OR)
result$log2OR_lower_95CI<-log2(result$OR_lower_95CI)
result$log2OR_upper_95CI<-log2(result$OR_upper_95CI)
result
