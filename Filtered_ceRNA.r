### -------------------------------------------------------------------------------
# Filtered lncRNA-miRNA-mRNA ceRNA relationship
### -------------------------------------------------------------------------------
setwd("F:/workspace/Filtered_ceRNA/data")

lnc_mir1 <- as.matrix(read.csv("miRNA-LncRNA.xls", sep="\t",header=T))
mir_gene1 <- as.matrix(read.csv("MiRNA-Target.xls", sep="\t",header=T))
lnc_info <- as.matrix(read.table("lnc_info.txt", sep="\t",header=F))
lncRNA_mir1 <- as.matrix(read.table("LUAD_lncRNA_mir_pcc.txt", sep="\t",header=F))
mir_pcg1 <- as.matrix(read.table("LUAD_mir_diff_pcg_pcc.txt", sep="\t",header=F))

lncRNA_mir1 <- lncRNA_mir1[-1,]
mir_pcg1 <- mir_pcg1[-1,]
lnc_mir2 <- unique(lnc_mir1[,c(3,1)])
mir_gene2 <- unique(mir_gene1[,c(1,2)])
lncRNA_mir1[,1] <- substr(lncRNA_mir1[,1],1,15)

# 将有靶点关系的lncRNA-miRNA对挑选出来，并将我们自己project中的lncRNA对应上去
lncRNA_symbol <- intersect(unique(lnc_mir2[,1]),unique(lnc_info[,7]))
cnt_lncRNA1 <- vector()
for(i in 1:length(lncRNA_symbol)){
	index <- vector()
	index <- which(lnc_mir2[,1]==lncRNA_symbol[i])
	cnt_lncRNA1 <- append(cnt_lncRNA1,index)
	index1 <- which(lnc_info[,7]==lncRNA_symbol[i])
	if(length(index1)>1){
	lnc_mir2[index,1] <- lnc_info[index1,5][1]
	}else{
	lnc_mir2[index,1] <- lnc_info[index1,5]
	}
}
lnc_mir3 <- lnc_mir2[cnt_lncRNA1,]

lnc <- intersect(unique(lnc_mir3[,1]),unique(lncRNA_mir1[,1]))
cnt_lnc <- vector()
lnc_mir_result <- vector()
for(i in 1:length(lnc)){
	index1 <- vector()
	index1 <- which(lncRNA_mir1[,1]==lnc[i])
	index2 <- vector()
	index2 <- which(lnc_mir3[,1]==lnc[i])
	mir <- intersect(lnc_mir3[index2,2],lncRNA_mir1[index1,2])
	cnt_index <- vector()
	if(length(mir)>0){
		lncRNA_mir2 <- lncRNA_mir1[index1,]
		for(j in 1:length(mir)){
			index <- vector()
			index <- which(lncRNA_mir2[,2]==mir[j])
			cnt_index <- append(cnt_index,index)
		}
		lnc_mir_result <- rbind(lnc_mir_result,lncRNA_mir2[cnt_index,])
	}else if(length(mir)==1){
			index <- vector()
			index <- which(lncRNA_mir2[,2]==mir)
			lnc_mir_result <- rbind(lnc_mir_result,lncRNA_mir2[index,])		
		}
}

index1 <- which(as.numeric(lnc_mir_result[,3])<0)
index2 <- which(as.numeric(lnc_mir_result[,5])<0.05)
index3 <- intersect(index1,index2)
length(index3)
lnc_mir_result1 <- lnc_mir_result[index3,]

# 通过miRNA将lncRNA-miRNA对与mRNA-miRNA关系对链接起来
pcg <- intersect(unique(mir_gene2[,2]),unique(mir_pcg1[,2]))
cnt_lnc <- vector()
mir_pcg_result <- vector()
for(i in 1:length(pcg)){
	index1 <- vector()
	index1 <- which(mir_pcg1[,2]==pcg[i])
	index2 <- vector()
	index2 <- which(mir_gene2[,2]==pcg[i])
	mir <- intersect(mir_gene2[index2,1],mir_pcg1[index1,1])
	cnt_index <- vector()
	if(length(mir)>0){
		lncRNA_mir2 <- mir_pcg1[index1,]
		for(j in 1:length(mir)){
			index <- vector()
			index <- which(lncRNA_mir2[,1]==mir[j])
			cnt_index <- append(cnt_index,index)
		}
		mir_pcg_result <- rbind(mir_pcg_result,lncRNA_mir2[cnt_index,])
	}else if(length(mir)==1){
			index <- vector()
			index <- which(lncRNA_mir2[,1]==mir)
			mir_pcg_result <- rbind(mir_pcg_result,lncRNA_mir2[index,])		
		}
}

index1 <- which(as.numeric(mir_pcg_result[,3])<0)
index2 <- which(as.numeric(mir_pcg_result[,5])<0.05)
index3 <- intersect(index1,index2)
length(index3)
mir_pcg_result1 <- mir_pcg_result[index3,]

setwd("F:/workspace/Filtered_ceRNA/result")
write.table(lnc_mir_result1, "lnc_mir_result.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(mir_pcg_result1, "mir_pcg_result.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(lnc_mir_result, "lnc_mir_result_all.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(mir_pcg_result, "mir_pcg_result_all.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)






