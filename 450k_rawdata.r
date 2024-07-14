"""
Processing raw data of methylation merge it into matrix
"""

setwd("F:/workspace/Identify_promoter_methylation/data")
probe_info1 <- as.matrix(read.table("LUAD.HumanMethylation450.19.lvl-3.TCGA-73-A9RS-01.txt",sep="\t",header=TRUE))
lnc_tss <- as.matrix(read.table("lnc_tss.txt",sep="\t"))

lnc_promotor <- lnc_tss[-1,]
colnames(lnc_promotor) <- c("lncRNA","chr","promotor_start","promotor_end","strand","TSS")
lnc_promotor[,3] <- as.numeric(lnc_promotor[,6]) - 2000
lnc_promotor[,4] <- as.numeric(lnc_promotor[,6]) + 2000
write.table(lnc_promotor,file="lnc_promotor.txt", sep="\t", quote=F,row.names=F,col.names=T)

# 获取探针位置信息
probe_info2 <- probe_info1[,c(1,4,5)]
probe_index <- which(as.numeric(probe_info2[,3])==0)
probe_info <- probe_info2[-probe_index,]
write.table(probe_info,file="probe_info.txt", sep="\t", quote=F,row.names=F,col.names=T)

lnc_probe_result <- c("lncRNA","probe","chr","probe_info")
# 探针对应的lncRNA
for(i in 1:nrow(lnc_promotor)){
	chr <- lnc_promotor[i,2]
	probe_chr <- paste("chr",probe_info[,2],sep="")
	chr_probe <- probe_info[probe_chr==chr,]
	lnc_probe_index1 <- which(as.numeric(chr_probe[,3])>=as.numeric(lnc_promotor[i,3]))
	lnc_probe_index2 <- which(as.numeric(chr_probe[,3])<=as.numeric(lnc_promotor[i,4]))
	lnc_probe_index <- intersect(lnc_probe_index1,lnc_probe_index2)
	
	if(length(lnc_probe_index)==1){
		lnc_probe_result1 <- vector()
		lnc_probe_result1 <- c(lnc_promotor[i,1],chr_probe[lnc_probe_index,])
		lnc_probe_result <- rbind(lnc_probe_result,lnc_probe_result1)
	}else if(length(lnc_probe_index)>1){
		lnc_probe_result1 <- matrix(0,length(lnc_probe_index),4)
		lnc_probe_result1[,1] <- lnc_promotor[i,1]
		lnc_probe_result1[,c(2,3,4)] <- chr_probe[lnc_probe_index,]
		lnc_probe_result <- rbind(lnc_probe_result,lnc_probe_result1)
	}	
}

colnames(lnc_probe_result) <- lnc_probe_result[1,]
lnc_probe_result1 <- lnc_probe_result[-1,]
lnc_probe_result1[,3] <- paste("chr",lnc_probe_result1[,3],sep="")

# 删除一个探针对应多个lncRNA的探针
probe_num <- table(lnc_probe_result[,2])
index3 <- which(probe_num > 1)
probe1 <- names(probe_num)[index3]
cnt_probe_index2 <- vector()
for(j in 1:length(probe1)){
	probe_index2 <- vector()
	probe_index2 <- which(lnc_probe_result1[,2]==probe1[j])
	cnt_probe_index2 <- append(cnt_probe_index2,probe_index2)
}
lnc_probe_result2 <- lnc_probe_result1[-cnt_probe_index2,]
write.table(lnc_probe_result2,file="lnc_probe_result.txt", sep="\t", quote=F,row.names=F,col.names=T)

# 获取每个样本的探针甲基化信息
sample_info <- dir()
lnc_methy <- matrix(0,nrow(lnc_probe_result2),2)
lnc_methy <- lnc_probe_result2[,c(1,2)]
sample_info1 <- substr(sample_info,46,73)
sample_info2 <- substr(sample_info1,1,15)

# 提取lncRNA启动子区域对应的探针甲基化
rownames(lnc_methy) <- lnc_methy[,2]
lnc_methy <- lnc_methy[order(rownames(lnc_methy)),]
all_probe_methy <- vector()
for(i in 1:length(sample_info)){
	sample1 <- read.table(sample_info[i],header=TRUE,skip=1,sep="\t")
	rownames(sample1) <- sample1[,1]
	sample2 <- sample1[lnc_methy[,2],2]
	#sample2 <- sample2[order(rownames(sample2)),]
	lnc_methy <- cbind(lnc_methy,sample2)      # lncRNA启动子区域对应的探针甲基化
	all_probe_methy <- cbind(all_probe_methy,sample1[,2])  # 所有探针甲基化
}
rownames(all_probe_methy) <- sample1[,1]
colnames(all_probe_methy) <- sample_info   
colnames(lnc_methy) <- c("lncRNA","probe",sample_info)

setwd("F:/workspace/Identify_promoter_methylation/result")

write.table(all_probe_methy,file="all_probe_methy.txt", sep="\t", quote=F,row.names=T,col.names=T)
write.table(lnc_methy,file="lnc_methy.txt", sep="\t", quote=F,row.names=F,col.names=T)
### --------------------------------------------------------------------------------------------------------------------




