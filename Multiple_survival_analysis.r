### --------------------------------------------------------------------------------------------------------
# Multiple regression analysis
### --------------------------------------------------------------------------------------------------------
rm(list=ls())
setwd("F:/LUAD/Multiple_regression_analysis/data")

# 处理训练集
train__survival1 <- as.matrix(read.csv("survival_information_result.txt",sep="\t"))
train_miRNA_exp1 <- as.matrix(read.csv("tumor_mir_exp.txt",sep="\t"))
train_gene_exp1 <- as.matrix(read.csv("LUAD_diff_pcg.txt",sep="\t"))
train_TF_exp1 <- as.matrix(read.csv("LUAD_lncRNA_exp.txt",sep="\t"))
miRNA_TF_gene <- as.matrix(read.table("significant_dhype_CNA_lncRNA_mir_pcg(0.5).txt",sep="\t",header=TRUE))

colnames(train_miRNA_exp1) <- substr(colnames(train_miRNA_exp1),1,12)
colnames(train_gene_exp1) <- substr(colnames(train_gene_exp1),1,12)
rownames(train_TF_exp1) <- substr(rownames(train_TF_exp1),1,15)

sam_sample1 <- intersect(colnames(train_miRNA_exp1),colnames(train_gene_exp1))
sam_sample2 <- intersect(sam_sample1,colnames(train_TF_exp1))
sam_sample3 <- intersect(sam_sample2,rownames(train__survival1))

train_miRNA_exp <- train_miRNA_exp1[,sam_sample3]
train_gene_exp <- train_gene_exp1[,sam_sample3]
train_TF_exp <- train_TF_exp1[,sam_sample3]
train__survival <- train__survival1[sam_sample3,]

train__survival <- train__survival[order(rownames(train__survival)),]
train_TF_exp <- train_TF_exp[,order(colnames(train_TF_exp))]
train_miRNA_exp <- train_miRNA_exp[,order(colnames(train_miRNA_exp))]
train_gene_exp <- train_gene_exp[,order(colnames(train_gene_exp))]

# 读取miRNA_TF_gene motif
cnt_len <- nrow(miRNA_TF_gene)
miRNA_TF_gene_risk <- matrix(0,cnt_len,10)
colnames(miRNA_TF_gene_risk) <- c("miRNA","lncRNA","gene","p_miRNA","coef_miRNA","p_TF","coef_TF","p_gene","coef_gene","train_p_value")
miRNA_TF_gene_risk[,1:3] <- miRNA_TF_gene
cnt_len1 <- nrow(miRNA_TF_gene_risk)

setwd("F:/LUAD/Multiple_regression_analysis/result")
for(j in 1:cnt_len1){
hub1 <- miRNA_TF_gene[j,2]
hub2 <- miRNA_TF_gene[j,1]
hub3 <- miRNA_TF_gene[j,3]

# 获取hub miRNA的训练集表达值
# miRNA_index <- which(train_miRNA_exp[,1] == hub1)
train_hub_miRNA <- as.character(train_miRNA_exp[hub1,])

# 获取hub tf的训练集表达谱
# tf_index <- which(train_TF_exp[,1] == hub2)
train_hub_tf <- as.character(train_TF_exp[hub2,])

# 获取hub gene的训练集表达谱
# gene_index <- which(train_gene_exp[,1] == hub3)
train_hub_gene <- as.character(train_gene_exp[hub3,])


# 加载cox回归以及生存分析的R包
library(splines)
library(survival)
train_time <- as.numeric(as.character(train__survival[,2]))
train_status <- as.numeric(as.character(train__survival[,1]))
x <- as.numeric(as.character(train_hub_miRNA))
y <- as.numeric(as.character(train_hub_tf))
z <- as.numeric(as.character(train_hub_gene))
train1 <- list(train_time,train_status,x) 
train2 <- list(train_time,train_status,y)
train3 <- list(train_time,train_status,z)

# 进行cox回归分析			  
q1<-summary(coxph(Surv(train_time, train_status) ~ x, train1)) 
cox_result1 <- q1$coef
q2<-summary(coxph(Surv(train_time, train_status) ~ y, train2)) 
cox_result2 <- q2$coef
q3<-summary(coxph(Surv(train_time, train_status) ~ z, train3)) 
cox_result3 <- q3$coef

miRNA_TF_gene_risk[j,4] <- cox_result1[1,5]
miRNA_TF_gene_risk[j,5] <- cox_result1[1,1]
miRNA_TF_gene_risk[j,6] <- cox_result2[1,5]
miRNA_TF_gene_risk[j,7] <- cox_result2[1,1]
miRNA_TF_gene_risk[j,8] <- cox_result3[1,5]
miRNA_TF_gene_risk[j,9] <- cox_result3[1,1]
# miRNA_TF_gene_risk[j,11] <-  q[[10]][3]
train_score <- matrix(0,2,ncol(train_TF_exp))

# 计算训练集的risk_score值
train_score[1,] <- cox_result1[1,1]*x + cox_result2[1,1]*y + cox_result3[1,1]*z
colnames(train_score) <- colnames(train_miRNA_exp)

# 计算risk_score的中位数，并把样本分成两类
cuff_train<-median(train_score[1,])
    index1<-train_score[1,]<=cuff_train
    index2<-!(index1)
    train_score[2,index1]<-1
    train_score[2,index2]<-2
train_score1 <- list(train_time,train_status,train_score[2,]) 

# 进行生存分析
name <- paste(unlist(strsplit(hub1,"hsa-"))[2],hub2,hub3,sep="_")
dif <- survdiff(Surv(train_time,train_status)~train_score[2,],train_score1)
    train_p <- 1-pchisq(dif[[5]],1)
	miRNA_TF_gene_risk[j,10] <- train_p
    pdf(file=paste(name,".pdf",sep=""))
    plot(survfit(Surv(train_time,train_status)~train_score[2,],train_score1),xscale=365.25,yaxt="n",col=c('forestgreen','orangered'),lwd=3,lty = "solid",xlab="time:years",main=name,ylab="Proportion survival")
    train_p1 <- paste("p-value=",round(train_p,4))
    text(16*0.7,0.9,train_p1,pos=4)
	train_low <- paste("low risk n=",as.numeric(table(train_score[2,]))[1])
	train_high <- paste("high risk n=",as.numeric(table(train_score[2,]))[2])
    legend(16*0.7,0.85,c(train_low,train_high),lwd=3,lty = "solid",col=c('forestgreen','orangered'))
   
	#xat=c(0,1000,2000,3000,4000,5000)
	#xlabels=c(0,1000,2000,3000,4000,5000)
	yat=c(0,0.2,0.4,0.6,0.8,1)
	ylabels=c(0.0,0.2,0.4,0.6,0.8,1.0)
	#axis(1, at = xat, labels = xlabels, tick = TRUE, lty = "solid",
          #lwd = 2, lwd.ticks = 2)
	axis(2, at = yat, labels = ylabels, tick = TRUE, lty = "solid",
          lwd = 2, lwd.ticks = 2)
		  		  
	dev.off()


}

write.table(miRNA_TF_gene_risk, "miRNA_lncRNA_gene_risk.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
### -------------------------------------------------------------------------------------------------------------------
