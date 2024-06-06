### ------------------------------------------------------------------------------
# Dealing with missing methylation values (Using mean value)
### ------------------------------------------------------------------------------
setwd("F:/workspace/Process_NA_methylation/data")
normal_probe_methy <- as.matrix(read.table("normal_probe_methy.txt",sep="\t"))
cancer_probe_methy <- as.matrix(read.table("cancer_probe_methy.txt",sep="\t"))
cancer_lnc_methy <- as.matrix(read.table("cancer_lnc_methy.txt",sep="\t"))
normal_lnc_methy <- as.matrix(read.table("normal_lnc_methy.txt",sep="\t"))

na.mat1 <- is.na(normal_probe_methy)
na.num1 <- apply(na.mat1, 1, sum)
cnt_len1 <- nrow(normal_probe_methy)
for(i in 1:cnt_len1){
	if(na.num1[i]>0){
		sub <- vector()
		data.complete <- vector()
		sub <- which(is.na(normal_probe_methy[i,])==T)
		data.complete <- as.numeric(normal_probe_methy[i,-sub])
		normal_probe_methy[i,sub] <- mean(data.complete)
	}
}

na.mat2 <- is.na(cancer_probe_methy)
na.num2 <- apply(na.mat2, 1, sum)
cnt_len2 <- nrow(cancer_probe_methy)
for(i in 1:cnt_len2){
	if(na.num2[i]>0){
		sub <- vector()
		data.complete <- vector()
		sub <- which(is.na(cancer_probe_methy[i,])==T)
		data.complete <- as.numeric(cancer_probe_methy[i,-sub])
		cancer_probe_methy[i,sub] <- mean(data.complete)
	}
}


na.mat3 <- is.na(normal_lnc_methy)
na.num3 <- apply(na.mat3, 1, sum)
cnt_len3 <- nrow(normal_lnc_methy)
for(i in 1:cnt_len3){
	if(na.num3[i]>0){
		sub <- vector()
		data.complete <- vector()
		sub <- which(is.na(normal_lnc_methy[i,])==T)
		data.complete <- as.numeric(normal_lnc_methy[i,-sub])
		normal_lnc_methy[i,sub] <- mean(data.complete)
	}
}


na.mat4 <- is.na(cancer_lnc_methy)
na.num4 <- apply(na.mat4, 1, sum)
cnt_len4 <- nrow(cancer_lnc_methy)
for(i in 1:cnt_len4){
	if(na.num4[i]>0){
		sub <- vector()
		data.complete <- vector()
		sub <- which(is.na(cancer_lnc_methy[i,])==T)
		data.complete <- as.numeric(cancer_lnc_methy[i,-sub])
		cancer_lnc_methy[i,sub] <- mean(data.complete)
	}
}

setwd("F:/workspace/Process_NA_methylation/result")
write.table(cancer_probe_methy,file="cancer_probe_methy.txt", sep="\t", quote=F,row.names=T,col.names=T)
write.table(normal_probe_methy,file="normal_probe_methy.txt", sep="\t", quote=F,row.names=T,col.names=T)
write.table(cancer_lnc_methy,file="cancer_lnc_methy.txt", sep="\t", quote=F,row.names=T,col.names=T)
write.table(normal_lnc_methy,file="normal_lnc_methy.txt", sep="\t", quote=F,row.names=T,col.names=T)


# 对探针进行t_test识别差异甲基化的探针
diff_probe_methy <- matrix(1,nrow(normal_probe_methy),4)
diff_probe_methy[,1] <- as.character(rownames(normal_probe_methy))
colnames(diff_probe_methy) <- c("probe","fold_change","p_value","FDR")
x_len1 <- nrow(normal_probe_methy)
for(i in 1:x_len1){
	t_value <- NULL
	t_value <- t.test(as.numeric(cancer_probe_methy[i,]),as.numeric(normal_probe_methy[i,]),alternative = c("two.sided"),paired = FALSE)
	mean_value <- as.numeric(t_value$estimate)
	diff_probe_methy[i,2] <- (mean_value[1])/(mean_value[2])
	diff_probe_methy[i,3] <- t_value$p.value
}

# 对p值进行多重检验矫正
p1 <- vector()
p1 <- as.numeric(diff_probe_methy[,3])
diff_probe_methy[,4] <- p.adjust(p1,method="BH",n=length(p1))


index1 <- which(as.numeric(diff_probe_methy[,4]) <= 0.01)
index2 <- which(as.numeric(diff_probe_methy[,2]) <= 0.5)
index3 <- which(as.numeric(diff_probe_methy[,2]) >= 2)
high_probe_index <- intersect(index1,index3)
low_probe_index <- intersect(index1,index2)
high_probe_methy <- diff_probe_methy[high_probe_index,]
low_probe_methy <- diff_probe_methy[low_probe_index,]


write.table(diff_probe_methy, "diff_probe_methy.txt", quote=FALSE, row.names=F, col.names=T, sep="\t")
write.table(high_probe_methy, "high_probe_methy.txt", quote=FALSE, row.names=F, col.names=T, sep="\t")
write.table(low_probe_methy, "low_probe_methy.txt", quote=FALSE, row.names=F, col.names=T, sep="\t")


# 对lncRNA进行t_test识别差异甲基化的lncRNA
diff_lnc_methy <- matrix(1,nrow(normal_lnc_methy),4)
diff_lnc_methy[,1] <- as.character(rownames(normal_lnc_methy))
colnames(diff_lnc_methy) <- c("lncRNA","fold_change","p_value","FDR")
x_len2 <- nrow(normal_lnc_methy)
for(i in 1:x_len2){
	t_value <- NULL
	t_value <- t.test(as.numeric(cancer_lnc_methy[i,]),as.numeric(normal_lnc_methy[i,]),alternative = c("two.sided"),paired = FALSE)
	mean_value <- as.numeric(t_value$estimate)
	diff_lnc_methy[i,2] <- (mean_value[1])/(mean_value[2])
	diff_lnc_methy[i,3] <- t_value$p.value
}

# 对p值进行多重检验矫正
p2 <- vector()
p2 <- as.numeric(diff_lnc_methy[,3])
diff_lnc_methy[,4] <- p.adjust(p2,method="BH",n=length(p2))


index1 <- which(as.numeric(diff_lnc_methy[,4]) <= 0.01)
index2 <- which(as.numeric(diff_lnc_methy[,2]) <= 0.5)
index3 <- which(as.numeric(diff_lnc_methy[,2]) >= 2)
high_lnc_index <- intersect(index1,index3)
low_lnc_index <- intersect(index1,index2)
high_lnc_methy <- diff_lnc_methy[high_lnc_index,]
low_lnc_methy <- diff_lnc_methy[low_lnc_index,]


write.table(diff_lnc_methy, "diff_lnc_methy.txt", quote=FALSE, row.names=F, col.names=T, sep="\t")
write.table(high_lnc_methy, "high_lnc_methy.txt", quote=FALSE, row.names=F, col.names=T, sep="\t")
write.table(low_lnc_methy, "low_lnc_methy.txt", quote=FALSE, row.names=F, col.names=T, sep="\t")


