setwd("F:/LUAD/5-1miRNA/hallmark_analysis/data")
diff_non_CNA_lncRNA_gene1 <- as.matrix(read.table("significant_dhype_CNA_lncRNA_pcg(0.4).txt", sep="\t",header=T))
hallmark <- as.matrix(readLines("hallmark.txt",n=-1))
#gene_info <- as.matrix(read.table("gene_information.txt", sep="\t",header=F))
#length(unique(diff_non_CNA_lncRNA_gene1[,2]))
#length(unique(intersect(diff_non_CNA_lncRNA_gene1[,2],gene_info[,7])))

dhype_non_CNA <- matrix(2,2,5)
colnames(dhype_non_CNA) <- c("hallmark","lncRNA","gene","P_value","FDR")
non_lnc <- unique(diff_non_CNA_lncRNA_gene1[,1])
non_gene <- unique(diff_non_CNA_lncRNA_gene1[,2])

hallmark_gene <- vector()
for(i in 1:length(hallmark)){
hallmark1 <- unlist(strsplit(hallmark[i],"\t"))[-2]
hallmark_gene <- append(hallmark_gene,hallmark1[-1])
}

for(i in 1:length(hallmark)){
	hallmark1 <- unlist(strsplit(hallmark[i],"\t"))[-2]
	for(j in 1:length(non_lnc)){
		index1 <- which(diff_non_CNA_lncRNA_gene1[,1]==non_lnc[j])
		gene <- vector()
		gene <- unique(diff_non_CNA_lncRNA_gene1[index1,2])
		sam_gene <- intersect(gene,hallmark1[-1])
		if(length(sam_gene)>0){
			x <- length(sam_gene)
			m <- length(gene)
			n <- 20351           #length(c(non_gene,hallmark_gene))       
			l <- length(unique(hallmark1[-1]))
			dhype_p <- vector()
			dhype_p <- phyper(x-1, m, n-m, l, lower.tail = FALSE, log.p = FALSE)
			#lnc_dhype <- vector()
			lnc_dhype <- matrix(0,x,5)
			lnc_dhype[,1] <- hallmark1[1]
			lnc_dhype[,2] <- non_lnc[j]
			lnc_dhype[,3] <- sam_gene
			lnc_dhype[,4] <- dhype_p
			#lnc_dhype[,5] <- 2
			#lnc_dhype <- c(hallmark1[1],non_lnc[j],dhype_p,"2")
			dhype_non_CNA <- rbind(dhype_non_CNA,lnc_dhype)
		}
	}
}
dhype_non_CNA <- dhype_non_CNA[c(-1,-2),]

dhype_lnc_CNA1 <- unique(dhype_non_CNA[,c(1,2,4,5)])

#dhype_lnc_CNA1[,4] <- p.adjust(as.numeric(dhype_lnc_CNA1[,3]),method="BH")
sig_dhype_lnc_CNA <- dhype_lnc_CNA1[as.numeric(dhype_lnc_CNA1[,3])<0.05,]
sig_dhype_lnc_pcg_CNA <- dhype_non_CNA[as.numeric(dhype_non_CNA[,4])<0.05,]

setwd("F:/LUAD/5-1miRNA/hallmark_analysis/result")
write.table(sig_dhype_lnc_CNA, "sig_dhype_hallmark_CNA_lncRNA(0.5).txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(sig_dhype_lnc_pcg_CNA, "sig_dhype_hallmark_CNA_lncRNA_pcg(0.5).txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)





