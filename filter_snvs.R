#------------------------------------------
# Script to filter variant calls - ALL/AML switch case
# Tim Coorens - August 2022 
#------------------------------------------
options(stringsAsFactors = F)

#----------------
# FUNCTIONS
#----------------

exact.binomial=function(gender,NV,NR,cutoff=-5){
  # Function to filter out germline variants based on unmatched
  # variant calls of multiple samples from same individual (aggregate coverage
  # ideally >150 or so, but will work with less). NV is matrix of reads supporting 
  # variants and NR the matrix with total depth (samples as columns, mutations rows, 
  # with rownames as chr_pos_ref_alt or equivalent). Returns a logical vector, 
  # TRUE if mutation is likely to be germline.
  
  XY_chromosomal = grepl("X|Y",rownames(NR))
  autosomal = !XY_chromosomal
  
  if(gender=="female"){
    NV_vec = rowSums(NV)
    NR_vec = rowSums(NR)
    pval = rep(1,length(NV_vec))
    for (n in 1:length(NV_vec)){
      if(NR_vec[n]>0){
        pval[n] = binom.test(x=NV_vec[n],
                             n=NR_vec[n],
                             p=0.5,alt='less')$p.value
      }
      if (n%%1000==0){
        print(n)
      }
    }
  }
  # For male, split test in autosomal and XY chromosomal part
  if(gender=="male"){
    pval=rep(1,nrow(NV))
    NV_vec = rowSums(NV)[autosomal]
    NR_vec = rowSums(NR)[autosomal]
    pval_auto = rep(1,sum(autosomal))
    pval_XY = rep(1,sum(XY_chromosomal))
    
    for (n in 1:sum(autosomal)){
      if(NR_vec[n]>0){
        pval_auto[n] = binom.test(x=NV_vec[n],
                                  n=NR_vec[n],
                                  p=0.5,alt='less')$p.value
      }
      if (n%%1000==0){
        print(n)
      }
    }
    NV_vec = rowSums(NV)[XY_chromosomal]
    NR_vec = rowSums(NR)[XY_chromosomal]
    for (n in 1:sum(XY_chromosomal)){
      if(NR_vec[n]>0){
        pval_XY[n] = binom.test(x=NV_vec[n],
                                n=NR_vec[n],
                                p=0.95,alt='less')$p.value
      }
      if (n%%1000==0){
        print(n)
      }
    }
    pval[autosomal]=pval_auto
    pval[XY_chromosomal]=pval_XY
  }
  qval = p.adjust(pval,method="BH")
  germline = log10(qval)>cutoff
  return(germline)
}

require(VGAM)
estimateRho_gridml = function(NV_vec,NR_vec) {
  #Make sure depth is non-zero
  NV_vec=NV_vec[NR_vec>0]
  NR_vec=NR_vec[NR_vec>0]
  
  # Function to estimate maximum likelihood value of rho for beta-binomial
  rhovec = 10^seq(-6,-0.05,by=0.05) # rho will be bounded within 1e-6 and 0.89
  mu=sum(NV_vec)/sum(NR_vec)
  ll = sapply(rhovec, function(rhoj) sum(dbetabinom(x=NV_vec, size=NR_vec, rho=rhoj, prob=mu, log=T)))
  return(rhovec[ll==max(ll)][1])
}

beta.binom.filter = function(NR,NV,cutoff=0.1, binom.pval=F,pval.cutoff=0.05){
  # Function to apply beta-binomial filter for artefacts. Works best on sets of
  # clonal samples (ideally >10 or so). As before, takes NV and NR as input. 
  # Optionally calculates pvalue of likelihood beta-binomial with estimated rho
  # fits better than binomial. This was supposed to protect against low-depth variants,
  # but use with caution. Returns logical vector with good variants = TRUE
  
  rho_est = pval = rep(NA,nrow(NR))
  for (k in 1:nrow(NR)){
    rho_est[k]=estimateRho_gridml(NV_vec = as.numeric(NV[k,]),
                                  NR_vec = as.numeric(NR[k,]))
    if (k%%1000==0){
      print(k)
    }
  }
  flt_rho=rho_est>=cutoff
  return(rho_est)
}

library("GenomicRanges")
library("Rsamtools")
library("MASS")

plot_spectrum = function(bed,save,add_to_title="",genomeFile="/nfs/cancer_ref01/Homo_sapiens/37/genome.fa"){
  mutations=as.data.frame(bed)
  colnames(mutations) = c("chr","pos","ref","mut")
  mutations$pos=as.numeric(mutations$pos)
  mutations = mutations[(mutations$ref %in% c("A","C","G","T")) & (mutations$mut %in% c("A","C","G","T")) & mutations$chr %in% c(1:22,"X","Y"),]
  mutations$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(mutations$chr, IRanges(as.numeric(mutations$pos)-1, 
                                                                                     as.numeric(mutations$pos)+1))))
  # 2. Annotating the mutation from the pyrimidine base
  ntcomp = c(T="A",G="C",C="G",A="T")
  mutations$sub = paste(mutations$ref,mutations$mut,sep=">")
  mutations$trinuc_ref_py = mutations$trinuc_ref
  for (j in 1:nrow(mutations)) {
    if (mutations$ref[j] %in% c("A","G")) { # Purine base
      mutations$sub[j] = paste(ntcomp[mutations$ref[j]],ntcomp[mutations$mut[j]],sep=">")
      mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
    }
  }
  
  # 3. Counting subs
  freqs = table(paste(mutations$sub,paste(substr(mutations$trinuc_ref_py,1,1),substr(mutations$trinuc_ref_py,3,3),sep="-"),sep=","))
  sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
  ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
  full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
  freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec
  
  xstr = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
  
  if(!is.null(save)){pdf(save,width=12,height=4)}
  if(is.null(save)){dev.new(width=12,height=4)}
  colvec = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
  y = freqs_full; maxy = max(y)
  h = barplot(y, las=2, col=colvec, border=NA, ylim=c(0,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab="Number mutations", main=paste0("Number of mutations: ",sum(freqs_full), add_to_title))
  for (j in 1:length(sub_vec)) {
    xpos = h[c((j-1)*16+1,j*16)]
    rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colvec[j*16])
    text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
  }    
  if(!is.null(save)){ dev.copy(pdf,save,width=12,height=4); dev.off(); dev.off()}
  
}

#----------------
# READ INPUT
#---------------- 

data = read.table("PD53609.snp.tsv", comment.char="",header=T)
data_allQ=read.table("PD53609.snp.allQ.tsv", comment.char="",header=T)

Muts = paste(data$Chrom,data$Pos,data$Ref,data$Alt,sep="_")
Genotype = data[,grepl("VAF",colnames(data))&colnames(data)!="PDv38is_wgs_VAF"]
NR = data[,grepl("DEP",colnames(data))&colnames(data)!="PDv38is_wgs_DEP"]
NR_allQ = data_allQ[,grepl("DEP",colnames(data_allQ))&colnames(data_allQ)!="PDv38is_wgs_DEP"]
NV = data[,grepl("MTR",colnames(data))&colnames(data)!="PDv38is_wgs_MTR"]
NV_allQ = data_allQ[,grepl("MTR",colnames(data_allQ))&colnames(data_allQ)!="PDv38is_wgs_MTR"]
rownames(Genotype)=rownames(NV)=rownames(NR)=rownames(NV_allQ)=rownames(NR_allQ)=Muts
samples = colnames(Genotype)=colnames(NR)=colnames(NV)=gsub("_VAF","",colnames(Genotype))
colnames(NV_allQ)=colnames(NR_allQ)=gsub("_DEP","",colnames(NR_allQ))

gender='female'

#Select samples with high degree
good_samples=paste0('PD53609',letters[3:10])

#----------------
# FILTER
#---------------- 
#Depth
Depth_filter = rowMeans(NR[,good_samples])>30&rowMeans(NR[,good_samples])<150
#Quality ratio
qual_ratio=rowSums(NV)/rowSums(NV_allQ)

indels=read.table("all_muts_sorted_unique_v2.bed")
colnames(indels)=c("Chr","Start","Ref","Alt")
indels$End=indels$Start+nchar(indels$Ref)-1

close_to_indel=rep(NA,nrow(data))
intv=1
for (n in 1:nrow(data)){
  close_to_indel[n]=any(data$Chrom[n]==indels$Chr&
                          data$Pos[n]>=(indels$Start-intv)&
                          data$Pos[n]<=(indels$End+intv))
  if(n%%1000==0) print(n)
}


Genotype_flt = Genotype[Depth_filter&!close_to_indel&qual_ratio>0.75,good_samples]
NR_flt = NR[Depth_filter&!close_to_indel&qual_ratio>0.75,good_samples]
NV_flt = NV[Depth_filter&!close_to_indel&qual_ratio>0.75,good_samples]
Muts=Muts[Depth_filter&!close_to_indel&qual_ratio>0.75]

CN_normal=paste0("PD53609",c('c','e','i','j','f'))
germline=exact.binomial(gender=gender,NV=NV_flt[,CN_normal],NR=NR_flt[,CN_normal],cutoff = -5) #determine which variants are germline


Genotype_flt2=Genotype_flt[!germline2&rowSums(NV_flt>4)>0,]
NV_flt2=NV_flt[!germline2&rowSums(NV_flt>4)>0,]
NR_flt2=NR_flt[!germline2&rowSums(NV_flt>4)>0,]

NR_flt_nonzero=NR_flt2
NR_flt_nonzero[NR_flt_nonzero==0]=1
shared_muts=rowSums(NV_flt2>0)>1

rho_est = beta.binom.filter(NR=NR_flt_nonzero[shared_muts,],NV=NV_flt2[shared_muts,])
flt_rho=rho_est<0.01
rho_filtered_out = rownames(NR_flt_nonzero[shared_muts,])[flt_rho]
write.table(rho_filtered_out,"snv_bbinom_filtered_out.txt")
write.table(rho_est,"rho_est.txt")

NR_flt3= NR_flt2[!rownames(NR_flt2)%in%rho_filtered_out,]
NV_flt3= NV_flt2[!rownames(NV_flt2)%in%rho_filtered_out,]
Genotype_flt3=Genotype_flt2[!rownames(Genotype_flt2)%in%rho_filtered_out,]

write.table(NR_flt3,"snp_NR_filtered_all.txt")
write.table(NV_flt3,"snp_NV_filtered_all.txt")