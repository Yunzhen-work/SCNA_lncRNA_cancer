### ----------------------------------------------------------------------------
# process gene expression data that downloaded from TCGA
### ----------------------------------------------------------------------------
gene_normal <- read.table("RPKM_normal.txt",head=T,sep="\t")
gene <- gene_normal

gene1_normal <- gene[,3:length(gene_normal[1,])]

index <- c()
j <- 1
for(i in 1:length(gene1_normal[,1]))
{
  zero <- length(gene1_normal[i,gene1_normal[i,]==0])
  sample_size <- length(gene1_normal[1,])
  if((zero/sample_size)>0.3)
  {
    index[j] <- i
    j <- j+1
  }
}

expe_gene_normal <- gene1_normal[-index,]

### 对表达值进行log2标准化处理 
for(i in 1:length(expe_gene_normal[,1]))
{
  expe_gene_normal[i,expe_gene_normal[i,]==0] <- 0.05
}
expe_gene_normal <- log2(expe_gene_normal)
expe_gene_normal1 <- expe_gene_normal

mean_normal <- apply(expe_gene_normal,1,mean)
sd_normal <- apply(expe_gene_normal,1,sd)

for(i in 1:length(expe_gene_normal[,1]))
{
  expe_gene_normal[i,expe_gene_normal[i,]>(mean_normal[i]+sd_normal[i])] <- 1
  expe_gene_normal[i,expe_gene_normal[i,]<(mean_normal[i]-sd_normal[i]) & expe_gene_normal[i,]!=1] <- -1
  expe_gene_normal[i,expe_gene_normal[i,]!=1 & expe_gene_normal[i,]!=-1] <- 0
}

expe_gene_normal <- cbind(gene[-index,1:2],expe_gene_normal)
expe_gene_normal1 <- cbind(gene[-index,1:2],expe_gene_normal1)

write.table(expe_gene_normal,"normal_exp_GE.txt",row.names=F,quote=F,sep="\t")
write.table(expe_gene_normal1,"normal_exp.txt",row.names=F,quote=F,sep="\t")

### ----------------------------------------------------------------------------