### ----------------------------------------------------------------------------
# 在R中用SAM包做差异表达分析
### ----------------------------------------------------------------------------

library(impute)
library(matrixStats)
library(samr)
setwd("F:/workspace/LUAD_diff/data")
basal_tumor_lncRNA <- as.matrix(read.table("LUAD_gene_exp.txt",sep="\t"))
normal_lncRNA_exp <- as.matrix(read.table("normal_gene_exp.txt",sep="\t"))

setwd("F:/workspace/LUAD_diff/result")
basal_lncRNA <- intersect(rownames(basal_tumor_lncRNA),rownames(normal_lncRNA_exp))
basal_index <- vector()
for(i in 1:length(basal_lncRNA)){
	index <- vector()
	index <- which(rownames(basal_tumor_lncRNA)==basal_lncRNA[i])
	basal_index <- append(index,basal_index)
}
basal_lncRNA1 <- basal_tumor_lncRNA[basal_index,]
basal_lncRNA1 <- basal_lncRNA1[sort(rownames(basal_lncRNA1)),]

basal_normal_index <- vector()
for(i in 1:length(basal_lncRNA)){
index <- vector()
index <- which(rownames(normal_lncRNA_exp)==basal_lncRNA[i])
basal_normal_index <- append(index,basal_normal_index)
}
basal_normal_lncRNA <- normal_lncRNA_exp[basal_normal_index,]
basal_normal_lncRNA <- basal_normal_lncRNA[sort(rownames(basal_normal_lncRNA)),] 

basal_lncRNA_exp <- cbind(basal_lncRNA1,basal_normal_lncRNA)
basal_sample_t_n <- c(rep(2,ncol(basal_lncRNA1)),rep(1,ncol(basal_normal_lncRNA)))

# 开始做差异表达分析
basal_samfit <- SAM(basal_lncRNA_exp,basal_sample_t_n,fdr.output=0.05,nperms=1000,resp.type="Two class unpaired",logged2 = TRUE,geneid=as.character(rownames(basal_lncRNA_exp)),genenames=as.character(rownames(basal_lncRNA_exp)))
basal_diffExp_lncRNA_up1 <- basal_samfit$siggenes.table$genes.up
basal_diffExp_lncRNA_down1 <- basal_samfit$siggenes.table$genes.lo
# basal_diffExp_lncRNA_up <- basal_diffExp_lncRNA_up1[as.numeric(basal_diffExp_lncRNA_up1[,6])>4,]
# basal_diffExp_lncRNA_down <- basal_diffExp_lncRNA_down1[as.numeric(basal_diffExp_lncRNA_down1[,6])<(1/4),]

write.table(basal_diffExp_lncRNA_up1,file="LUAD_diffExp_gene_up.txt", sep="\t", quote=F,row.names=F,col.names=T)
write.table(basal_diffExp_lncRNA_down1,file="LUAD_diffExp_gene_down.txt", sep="\t", quote=F,row.names=F,col.names=T)



