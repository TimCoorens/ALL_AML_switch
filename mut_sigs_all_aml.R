#------------------------------------------
# Script mutational signature analysis - ALL/AML switch case
# Tim Coorens - August 2022 
#-----------------------------------------
options(stringsAsFactors = F)
library(sigfit)
muts_df=read.table("muts_df_2022_07_26.txt") 
muts_df=muts_df[muts_df$Branch!="Noise",]
muts_df$Branch[muts_df$Branch=="Primary, major clone - neutral tail leading to relapse root"]="Primary, major clone"  

Muts_all=as.data.frame(muts_df[,c("Chr","Pos","Ref","Alt","Branch")])
colnames(Muts_all)=c("Chr","Pos","Ref","Alt","SampleID")
trinuc_mut_mat=mutlist_to_96_contexts(Muts_all,version=38)


library(RColorBrewer)
ref=read.csv("/lustre/scratch117/casm/team268/tc16/COSMIC_Mutational_Signatures_v3.1.csv")
rownames(ref)=paste0(substr(ref$Subtype,1,1),"[",ref$Type,"]",substr(ref$Subtype,3,3))
ref=ref[,-c(1,2)]
sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
full_vec = paste0(rep(c("A","C","G","T"),each=4),"[",rep(sub_vec,each=16),"]",rep(c("A","C","G","T"),times=4))
ref=ref[full_vec,]
ref=apply(ref,2,as.numeric)
ref[is.na(ref)|ref==0]=0.00001
ref=t(ref)/colSums(ref)
data("cosmic_signatures_v3")

extract=extract_signatures(trinuc_mut_mat,nsignatures = 2:5)
signatures <- retrieve_pars(extract[[3]],
                            par = "signatures")

extract=extract_signatures(trinuc_mut_mat,nsignatures = 3,iter = 20000,warmup = 10000)
exposures <- retrieve_pars(extract, "exposures")
signatures <- retrieve_pars(extract,"signatures")

saveRDS(exposures,"sigfit_blina_exposures.Rdata")
saveRDS(signatures,"sigfit_blina_signatures.Rdata")



exposures <- retrieve_pars(extract, "exposures")
signatures <- retrieve_pars(extract,"signatures")
match_signatures(signatures, ref)

par(mfrow = c(3, 1))
sigfit::plot_spectrum(signatures)

colnames(ref)=colnames(cosmic_signatures_v3)

sig_mat_ord=exposures$mean[c("Primary","Primary, major clone","Primary, minor clone in d and i",
                             "Primary, minor clone in i","Primary, minor clone in d","Relapse root",
                             "Shared ALL-relapses","All-relapse 1 private","ALL-relapse 2",
                             "AML Branch","Private to f" ,"Shared g and h","Private to g","Private to h"),paste0("Signature ",c("A","C","B"))]
select=rowSums(trinuc_mut_mat[rownames(sig_mat_ord),])>100

pdf("sig_barplot_2022_07_26.pdf")
barplot(t(sig_mat_ord[select,]),col=c("grey80","grey50","grey20"),las=2)
dev.off()

