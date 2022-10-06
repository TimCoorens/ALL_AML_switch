#------------------------------------------
# Plot heatmap - ALL/AML switch case
# Tim Coorens - August 2022 
#------------------------------------------


library(pheatmap)
muts_df=read.table("~/Desktop/Blina/muts_df_2022_07_26.txt") 
muts_df=muts_df[muts_df$Branch!="Noise",]
muts_df$Branch[muts_df$Branch=="Primary, major clone - neutral tail leading to relapse root"]="Primary, major clone"  

annot_branch <- data.frame(Branch = muts_df$Branch)
rownames(annot_branch) <- rownames(muts_df)
muts_df$Branch=factor(muts_df$Branch,levels=c("Primary","Primary, minor clone in d and i",
                                              "Primary, minor clone in i","Primary, minor clone in d",
                                              "Primary, major clone","Relapse root",
                                              "Shared ALL-relapses","All-relapse 1 private","ALL-relapse 2",
                                              "AML Branch","Private to f" ,"Shared g and h","Private to g","Private to h"))
muts_df=muts_df[order(muts_df$Branch),]
sample_annot=data.frame(Sample=c("ALL - relapse 1", "Primary ALL","ALL - relapse 1","Bone marrow","AML","AML","Primary ALL","ALL - relapse 2"))
rownames(sample_annot)=colnames(muts_df[,5:12])

primary_cols=brewer.pal(9,"Purples")[c(6,8)]
subclone_cols=brewer.pal(9,"Greens")[c(4,6,8)]
ALL_cols=brewer.pal(9,"Blues")[c(4,6,8)]
AML_cols=brewer.pal(9,"Reds")[c(3,4,6,7,8)]

ann_colors = list(
  Sample = c("Primary ALL" = "black", "ALL - relapse 1" = "royalblue3","ALL - relapse 2"="skyblue1", "Bone marrow"="salmon2","AML" = "firebrick"),
  Branch=c("Primary"="black", "Primary, minor clone in d and i" = subclone_cols[3],"Primary, minor clone in i"= subclone_cols[2],
           "Primary, minor clone in d"= subclone_cols[1],"Primary, major clone" = "goldenrod1","Relapse root"="darkorchid3",
           "Shared ALL-relapses"=ALL_cols[3],"All-relapse 1 private"=ALL_cols[2],"ALL-relapse 2"=ALL_cols[1],
           "AML Branch"=AML_cols[5],"Private to f"="grey60" ,"Shared g and h"=AML_cols[3],"Private to g"=AML_cols[2],"Private to h"=AML_cols[1])
)

library(RColorBrewer)
pdf("~/Desktop/Blina/big_heatmap.pdf",width=12,height=8)
pheatmap(t(muts_df[,paste0("PD53609",rev(c("d","i","c","e","j","f","g","h")))]),
         annotation_col=annot_branch,
         annotation_row=sample_annot,
         show_colnames=F, 
         treeheight_row=0, 
         treeheight_col=0,
         cluster_rows=F, 
         cluster_cols=F,
         color = c("white",brewer.pal(9,"Greys")[2],colorRampPalette(brewer.pal(9,"Greys")[3:9])(99)),
         annotation_colors = ann_colors)
dev.off()