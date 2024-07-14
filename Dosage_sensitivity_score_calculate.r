### --------------------------------------------------------------------------------------------------
# 计算每个lncRNA的dosage sensitivity score
### ----------------------------------------------------------------------------------------------------

library(ggplot2)
setwd("F:/workspace/DSS/data")

# 读入数据
gene_cnv<-read.table("LUAD_continue_CNA.txt",sep="\t")
gene_cnv<-as.matrix(gene_cnv)
gene_exp<-read.table("LUAD_CNA_exp.txt",sep="\t")
gene_exp<-as.matrix(gene_exp)

# gene_exp <-read.table("gene_expr.txt",sep="\t")
# gene_exp<-as.matrix(gene_exp)

sen1 <- gene_cnv
sen2 <- gene_exp

tt<-read.table("LUAD_scatter_CNA.txt",sep="\t",header=T)
rr<-tt[rowSums(tt!=0)>=0.2*ncol(tt),]

ind<-intersect(rownames(rr),rownames(sen1))
sen1<-sen1[ind,]
sen2<-sen2[ind,]

r<-nrow(sen2)
DSS_matrix <- matrix(0,r,4)
colnames(DSS_matrix) <- c("gene_symbol","monotony_score","slope","dosage sensitive score")
DSS_matrix[,1] <- rownames(sen2)
# monotony<-logical()
n <- 6
for(i in 1:r){
fit <- loess(sen2[i,] ~ sen1[i,],span=2/3,degree=1,family="symmetric")
x<-seq(min(sen1[i,]),max(sen1[i,]),length.out=n)
px <- predict(fit, newdata=x)
# monotony[i]<-all(px == cummax(px))
z<-outer(px,px,'-')
z<-z[lower.tri(z)]
z<-sign(z)
MS <- 2*sum(z)/(n*(n-1))
slope<-numeric()
lma<-numeric()
lma<-lm(sen2[i,] ~ sen1[i,])
slope<-as.numeric(lma$coefficients[2])
DSS <- MS*abs(slope)
DSS_matrix[i,2] <- MS
DSS_matrix[i,3] <- slope
DSS_matrix[i,4] <- DSS
}
# min_value <- min(as.numeric(DSS_matrix[,4]))
# max_value <- max(as.numeric(DSS_matrix[,4]))
# DSS_matrix[,5] <- as.numeric(lapply(as.numeric(DSS_matrix[,4]),function(x) (x-min_value)/(max_value-min_value)))

setwd("F:/workspace/DSS/output_result")
write.table(DSS_matrix,"LUAD_lncRNA_DSS_matrix(n=6).txt", row.names=F, col.names=T, sep="\t", quote=F)
# DSS_matrix <- as.numeric(DSS_matrix[,4])
# hist(as.numeric(DSS_matrix[,4]))


	
	
	
	
