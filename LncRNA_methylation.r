"""
Identifying methylation in lncRNA
"""

setwd("F:/workspace/LncRNA_methylation/data")

# 读入数据
all_probe_methy <- as.matrix(read.table("all_probe_methy.txt",sep="\t"))
lnc_methy <- as.matrix(read.table("lnc_methy.txt",sep="\t",header=TRUE))
file_sample <- as.matrix(read.table("FILE_SAMPLE_MAP.txt",sep="\t",header=TRUE))
file_sample[,2] <- substr(file_sample[,2],1,15)
file_sample[,1] <- gsub("-",".",file_sample[,1])
rownames(file_sample) <- file_sample[,1]
all_probe_methy <- all_probe_methy[,order(colnames(all_probe_methy))]

# 处理文件
same_sample <- intersect(colnames(all_probe_methy),rownames(file_sample))
file_sample1 <- file_sample[same_sample,]
file_sample2 <- file_sample1[order(file_sample1[,1]),]
colnames(all_probe_methy) <- as.character(file_sample2[,2])

rownames(lnc_methy) <- lnc_methy[,2]
lnc_methy <- lnc_methy[,-2]
lnc_methy <- lnc_methy[,order(colnames(lnc_methy))]

same_sample1 <- intersect(colnames(lnc_methy),rownames(file_sample))
file_sample3 <- file_sample[same_sample1,]
file_sample4 <- file_sample1[order(file_sample3[,1]),]
colnames(lnc_methy) <- c(as.character(file_sample4[,2]),"lncRNA")
lnc_methy <- cbind(lnc_methy[,ncol(lnc_methy)],lnc_methy)
lnc_methy <- lnc_methy[,-ncol(lnc_methy)]
colnames(lnc_methy) <- c("lncRNA",as.character(file_sample4[,2]))

# 删除文件中值全为NA的探针
k <- vector()
for(z in 1:nrow(all_probe_methy)){
	tf <- unique(is.na(all_probe_methy[z,]))
	len <- length(tf)
	if(len==1&&tf==TRUE){
		k <- append(k,z)
	}
}
all_probe_methy1 <- all_probe_methy[-k,]

m <- vector()
for(h in 1:nrow(lnc_methy)){
	tf <- unique(is.na(lnc_methy[h,-1]))
	len <- length(tf)
	if(len==1&&tf==TRUE){
		m <- append(m,h)
	}
}
lnc_methy1 <- lnc_methy[-m,]

normal_sample <- colnames(all_probe_methy)[as.character(substr(colnames(all_probe_methy),14,14))=="1"]
normal_sample <- as.character(normal_sample)
cancer_sample <- setdiff(colnames(all_probe_methy),normal_sample)

# 获取lncRNA甲基化
lncRNA <- unique(lnc_methy1[,1])
lnc_methy2 <- matrix(0,length(lncRNA),(ncol(lnc_methy1)-1))
rownames(lnc_methy2) <- lncRNA
colnames(lnc_methy2) <- colnames(lnc_methy1)[-1]
for(i in 1:length(lncRNA)){
	index <- vector()
	index <- which(lnc_methy1[,1]==lncRNA[i])
	result <- lnc_methy1[index,-1]
	if(length(index)>1){
	result <- matrix(as.numeric(result),nrow=nrow(result)) 
	lnc_methy2[i,] <- apply(result,2,mean)
	}else{
	lnc_methy2[i,] <- as.character(result)
	}
}

cancer_probe_methy <- all_probe_methy1[,cancer_sample]
normal_probe_methy <- all_probe_methy1[,normal_sample]

cancer_lnc_methy <- lnc_methy2[,cancer_sample]
normal_lnc_methy <- lnc_methy2[,normal_sample]

setwd("F:/workspace/LncRNA_methylation/result")
write.table(cancer_probe_methy,file="cancer_probe_methy.txt", sep="\t", quote=F,row.names=T,col.names=T)
write.table(normal_probe_methy,file="normal_probe_methy.txt", sep="\t", quote=F,row.names=T,col.names=T)
write.table(cancer_lnc_methy,file="cancer_lnc_methy.txt", sep="\t", quote=F,row.names=T,col.names=T)
write.table(normal_lnc_methy,file="normal_lnc_methy.txt", sep="\t", quote=F,row.names=T,col.names=T)






