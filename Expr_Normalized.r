"""
TCGA 表达谱数据的标准化处理
"""

setwd("E:/LUAD/6gene表达/2表达数据标准化/data")
library(impute)

# 处理basal_exp表达数据
basal_exp1 <- as.matrix(read.table("LUAD_tumor_gene.txt", sep="\t"))
normal_exp1 <- as.matrix(read.table("LUAD_normal_gene.txt", sep="\t"))
gene_info <- as.matrix(read.table("gene_information.txt", sep="\t",header=F))

basal_exp <- basal_exp1
basal_exp <- matrix(as.numeric(basal_exp),nrow=nrow(basal_exp))
colnames(basal_exp) <- colnames(basal_exp1)
rownames(basal_exp) <- rownames(basal_exp1)

# 统计每个基因的表达值是0的数目30%
row_index <- NULL
for(i in 1:nrow(basal_exp)){
	a <- table(as.numeric(basal_exp[i,])==0)
	if(length(a)==2&&a[[2]]>ncol(basal_exp)*0.3){
		row_index <- append(row_index,i)
	}
	else if(length(a)==1&&unique(as.numeric(basal_exp[i,])==0)==TRUE){
		row_index <- append(row_index,i)
	}
}
row_index <- as.numeric(row_index)
if(length(row_index)>0){
	basal_exp <- basal_exp[-row_index,]
}
setwd("E:/LUAD/6gene表达/2表达数据标准化/result")
# 使用2代替基因表达值为0的样本
# basal_exp[basal_exp==0] <- 0.05
basal_exp <- basal_exp + 2
normal_exp <- normal_exp1[rownames(basal_exp),]
normal_exp <- normal_exp + 2
basal_exp <- log2(basal_exp)
normal_exp <- log2(normal_exp)

write.table(basal_exp,"LUAD_gene_exp.txt", row.names=T, col.names=T, sep="\t", quote=F)
write.table(normal_exp,"normal_gene_exp.txt", row.names=T, col.names=T, sep="\t", quote=F)

basal_diff_gene <- basal_exp
basal_seq <- seq(from=1,to=(2*nrow(basal_diff_gene)),by=2)
basal_diff_gene1 <- unlist(strsplit(rownames(basal_diff_gene),"\\|"))[basal_seq]

gene_sum <- table(basal_diff_gene1)
gene_index <- which(gene_sum>=2)
gene_name <- names(gene_sum[gene_index])

cnt_index <- vector()
gene_matrix <- matrix(0,length(gene_name),ncol(basal_diff_gene))
for(i in 1:length(gene_name)){
	index <- vector()
	index <- which(basal_diff_gene1==gene_name[i])
	gene_matrix[i,] <- apply(basal_diff_gene[index,],2,mean)
	cnt_index <- append(cnt_index,index)
}
rownames(gene_matrix) <- gene_name
basal_diff_gene2 <- basal_diff_gene[-cnt_index,]
basal_seq1 <- seq(from=1,to=(2*nrow(basal_diff_gene2)),by=2)
rownames(basal_diff_gene2) <- unlist(strsplit(rownames(basal_diff_gene2),"\\|"))[basal_seq1]
basal_diff_gene3 <- rbind(basal_diff_gene2,gene_matrix)
write.table(basal_diff_gene3,"LUAD_gene_exp1.txt", sep="\t", quote=F,row.names=T,col.names=T)

normal_diff_normal_gene <- normal_exp
normal_seq <- seq(from=1,to=(2*nrow(normal_diff_normal_gene)),by=2)
normal_diff_normal_gene1 <- unlist(strsplit(rownames(normal_diff_normal_gene),"\\|"))[normal_seq]

normal_gene_sum <- table(normal_diff_normal_gene1)
normal_gene_index <- which(normal_gene_sum>=2)
normal_gene_name <- names(normal_gene_sum[normal_gene_index])

cnt_index <- vector()
normal_gene_matrix <- matrix(0,length(normal_gene_name),ncol(normal_diff_normal_gene))
for(i in 1:length(normal_gene_name)){
	index <- vector()
	index <- which(normal_diff_normal_gene1==normal_gene_name[i])
	normal_gene_matrix[i,] <- apply(normal_diff_normal_gene[index,],2,mean)
	cnt_index <- append(cnt_index,index)
}
rownames(normal_gene_matrix) <- normal_gene_name
normal_diff_normal_gene2 <- normal_diff_normal_gene[-cnt_index,]
normal_seq1 <- seq(from=1,to=(2*nrow(normal_diff_normal_gene2)),by=2)
rownames(normal_diff_normal_gene2) <- unlist(strsplit(rownames(normal_diff_normal_gene2),"\\|"))[normal_seq1]
normal_diff_normal_gene3 <- rbind(normal_diff_normal_gene2,normal_gene_matrix)
write.table(normal_diff_normal_gene3,"normal_gene_exp1.txt", sep="\t", quote=F,row.names=T,col.names=T)

pcg <- intersect(rownames(basal_diff_gene3),gene_info[,7])
pcg_exp <- basal_diff_gene3[pcg,]
normal_pcg_exp <- normal_diff_normal_gene3[pcg,]

write.table(pcg_exp,"LUAD_pcg_exp.txt", row.names=T, col.names=T, sep="\t", quote=F)
write.table(normal_pcg_exp,"normal_pcg_gene_exp.txt", row.names=T, col.names=T, sep="\t", quote=F)


