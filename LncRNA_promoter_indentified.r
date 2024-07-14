"""
Identifying promoter region of lncRNA
"""

rm(list=ls())
# 距离lncRNA起始位点最近的TSS作为该lncRNA的TSS
setwd("E:/LUAD/lncRNA_promoter_methylation/data")
lnc_info <- as.matrix(read.table("lnc_info.txt",sep="\t"))
cna_lnc <- as.matrix(read.table("LUAD_lncRNA_DSS_matrix.txt",sep="\t",header=TRUE))
TSS_info <- as.matrix(read.table("Tss_all.txt",sep="\t"))

cnt_lnc <- nrow(lnc_info)
lnc_tss <- matrix(0,cnt_lnc,6) 
lnc_tss[,-6] <- lnc_info[,c(5,1,2,3,4)]
colnames(lnc_tss) <- c("lncRNA","chr","start","end","strand","TSS")
for(i in 1:cnt_lnc){
	chr <- lnc_tss[i,2]
	chr_index <- which(TSS_info[,1]==chr)
	chr_tss <- TSS_info[chr_index,]
	lnc_index <- which(chr_tss[,2]<=lnc_tss[i,3])
	if(length(lnc_index)>0){
	lnc_tss[i,6] <- chr_tss[lnc_index[length(lnc_index)],2]
	}
}
setwd("E:/LUAD/lncRNA_promoter_methylation/result")
write.table(lnc_tss,file="lnc_tss.txt", sep="\t", quote=F,row.names=F,col.names=T)


