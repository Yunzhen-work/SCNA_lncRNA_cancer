### ------------------------------------------------------------------------------------------------
# 寻找lncRNA的临近基因，并计算lncRNA与临近基因的相关性
### ------------------------------------------------------------------------------------------------
# 设置工作目录，并把matrix依次读入
setwd("F:/workspace/lncRNA_process/data")
lnc_info <- as.matrix(read.table("lnc_info.txt", sep="\t",header=F))
gene_info <- as.matrix(read.table("gene_information.txt", sep="\t",header=F))
lnc_gene <- as.matrix(read.table("LUAD_SCNA_lncRNA_diff_pcg(0.6).txt", sep="\t",header=F))
amp_up <- as.matrix(read.table("amp_up_lncRNA.txt", sep="\t",header=F))
del_down_lncRNA <- as.matrix(read.table("del_down_lncRNA.txt", sep="\t",header=F))
des_lncRNA <- as.matrix(read.table("LUAD_lncRNA_DSS_matrix(n=6).txt", sep="\t",header=T))

# 设置输出目录
setwd("F:/workspace/lncRNA_process/output_result")
lnc_near_gene <- matrix("gene",2,2)
lnc_near_gene[,1] <- "lncRNA"
cnt_lnc_index <- vector()
lncRNA <- unique(lnc_gene[,1])
lncRNA <- des_lncRNA[,1]

# 开始识别lncRNA的临近基因
for(i in 1:length(lncRNA)){
	index <- vector()
	index <- which(lnc_info[,1]==lncRNA[i])
	cnt_lnc_index <- append(cnt_lnc_index,index)
	chr_index <- vector()
	chr_index <- which(gene_info[,1]==lnc_info[index,2])
	gene_info1 <- gene_info[chr_index,]
	stain_index <- vector()
	stain_index <- which(gene_info1[,4]==lnc_info[index,5])
	gene_info2 <- gene_info1[stain_index,]
	
# 临近基因的识别条件：lncRNA上下游300kb内
	lnc_upstream <- as.numeric(lnc_info[index,3]) - 300000   ###########300kb
	lnc_downstream <- as.numeric(lnc_info[index,4]) + 300000 
	upstream_start_index <- which(as.numeric(gene_info2[,2])>=lnc_upstream)
	downstream_start_index <- which(as.numeric(gene_info2[,2])<=lnc_downstream)
	start_index <- intersect(upstream_start_index,downstream_start_index)
	
	if(length(start_index)>0){
	lncRNA_gene1 <- matrix(lncRNA[i],length(start_index),2)
	lncRNA_gene1[,2] <- gene_info2[start_index,7]
	lnc_near_gene <- rbind(lnc_near_gene,lncRNA_gene1)
	}
	
	upstream_end_index <- which(as.numeric(gene_info2[,3])>=lnc_upstream)
	downstream_end_index <- which(as.numeric(gene_info2[,3])<=lnc_downstream)
	end_index <- intersect(upstream_end_index,downstream_end_index)
	
	if(length(end_index)>0){
	lncRNA_gene1 <- matrix(lncRNA[i],length(end_index),2)
	lncRNA_gene1[,2] <- gene_info2[end_index,7]
	lnc_near_gene <- rbind(lnc_near_gene,lncRNA_gene1)
	}
}
des_lnc_info <- lnc_info[cnt_lnc_index,c(2,3,4,1)]
lnc_near_gene <- unique(lnc_near_gene)

# 写出结果
write.table(lnc_near_gene, "des_lncRNA_nearby_diff_gene.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(des_lnc_info, "des_lnc_info.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

lnc_gene_pcc <- matrix("lncRNA",2,5)
lnc_gene_pcc[,2] <- "gene"
# 使用intersect交叉合并两个matrix
lncRNA1 <- intersect(unique(lnc_near_gene[,1]),lnc_gene[,1])
for(i in 1:length(lncRNA1)){
	index1 <- vector()
	index1 <- which(lnc_gene[,1]==lncRNA1[i])
	index2 <- vector()
	index2 <- which(lnc_near_gene[,1]==lncRNA1[i])
	sam_gene <- vector()
	sam_gene <- intersect(unique(lnc_near_gene[index2,2]),lnc_gene[index1,2])
	lnc_gene1 <- lnc_gene[index1,]
	if(length(sam_gene)>0){
		if(length(index1)==1){
		lnc_gene_pcc <- rbind(lnc_gene_pcc,lnc_gene1)
		}else{
			for(j in 1:length(sam_gene)){
				index3 <- vector()
				index3 <- which(lnc_gene1[,2]==sam_gene[j])
				lnc_gene_pcc <- rbind(lnc_gene_pcc,lnc_gene1[index3,])
			}
        }		
	}
}

lnc_gene_pcc <- unique(lnc_gene_pcc)
lnc_gene_pcc[1,c(3,4,5)] <- c("pcc","p","FDR_BH")
colnames(lnc_gene_pcc) <- lnc_gene_pcc[1,]
lnc_gene_pcc <- lnc_gene_pcc[-1,]
write.table(lnc_gene_pcc, "des_lncRNA_nearby_diff_gene_pcc(0.6).txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)




