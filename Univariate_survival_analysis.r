"""
Univariate cox analysis for lncRNA
"""

rm(list=ls())
setwd("F:/LUAD/Survival_analysis/data")

# processing training dataset
survival_info1 <- as.matrix(read.csv("survival_information_result.txt",sep="\t"))
LUAD_lncRNA_exp1 <- as.matrix(read.csv("LUAD_lncRNA_exp.txt",sep="\t"))
hallmark_lncRNA1 <- as.matrix(read.csv("significant_dhype_CNA_lncRNA_mir_pcg(0.5).txt",sep="\t",header=TRUE))

rownames(LUAD_lncRNA_exp1) <- substr(rownames(LUAD_lncRNA_exp1),1,15)
same_sample <- intersect(rownames(survival_info1),colnames(LUAD_lncRNA_exp1))
survival_info <- survival_info1[same_sample,]
LUAD_lncRNA_exp <- LUAD_lncRNA_exp1[,same_sample]
survival_info <- survival_info[order(rownames(survival_info)),]
LUAD_lncRNA_exp <- LUAD_lncRNA_exp[,order(colnames(LUAD_lncRNA_exp))]
sample_num <- length(same_sample)


setwd("F:/LUAD/Survival_analysis/result")
write.table(survival_info, "survival_info.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(LUAD_lncRNA_exp, "LUAD_lncRNA_exp.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

hallmark_lncRNA <- unique(hallmark_lncRNA1[,1])
lncRNA_risk <- matrix(0,length(hallmark_lncRNA),4)
colnames(lncRNA_risk) <- c("hallmark","p_lncRNA","coef_lncRNA","p_value")
lncRNA_risk[,1] <- hallmark_lncRNA
cnt_len1 <- nrow(lncRNA_risk)


for(j in 1:length(hallmark_lncRNA)){

#j <- 6
#j <- j+1
#hallmark_index <- which(hallmark_lncRNA1[,1]==hallmark_lncRNA[j])
#ceRNA_lncRNA <- hallmark_lncRNA1[hallmark_index,2]
#ceRNA_lncRNA <- unique(ceRNA_lncRNA1[,1])
#cnt_len1 <- length(ceRNA_lncRNA

#hub <- ceRNA_lncRNA
# Obtaining the expressing of hub lncRNAs
hub <- hallmark_lncRNA[j]
lncRNA_exp <- LUAD_lncRNA_exp[hub,]

# Loading the R package for cox regression and survival analysis
library(splines)
library(survival)
train_time <- as.numeric(as.character(survival_info[,2]))
train_status <- as.numeric(as.character(survival_info[,1]))
#x <- as.numeric(as.character(lncRNA_exp[1,]))
#y <- as.numeric(as.character(lncRNA_exp[2,]))
#z <- as.numeric(as.character(lncRNA_exp[3,]))
#w <- as.numeric(as.character(lncRNA_exp[4,]))

x <- as.numeric(as.character(lncRNA_exp))
# train <- list(train_time,train_status,x,y,z,w)
train <- list(train_time,train_status,x)

# Beginning cox regression analysis
# q<-summary(coxph(Surv(train_time, train_status) ~ x + y+z+w , train)) 
q<-summary(coxph(Surv(train_time, train_status) ~ x , train)) 
q

cox_result <- q$coef
train_p1 <- round(q[[9]][3],4)
lncRNA_risk[j,2] <- train_p1
# lncRNA_risk[j,2] <- cox_result[1,5]
lncRNA_risk[j,3] <- cox_result[1]

# Calculating the risk_score value for training set
train_score <- matrix(0,2,sample_num)
#train_score[1,] <- cox_result[1,1]*x +cox_result[2,1]*y+cox_result[3,1]*z+cox_result[4,1]*w
train_score[1,] <- cox_result[1,1]*x 

colnames(train_score) <- colnames(LUAD_lncRNA_exp)
# Calculating the median of risk_score and dividing the samples into two groups
cuff_train<-median(train_score[1,])
    index1<-train_score[1,]<=cuff_train
    index2<-!(index1)
    train_score[2,index1]<-1
    train_score[2,index2]<-2

# Deleting the patient samples that surviving over 10 years
# sample_index <- which(train_time>=3652.5)
train_time2 <- train_time
train_status2 <- train_status
train_score2 <- train_score	
train_score1 <- list(train_time2,train_status2,train_score2[2,]) #进行生存分析

# Beginning survival analysis
# name <- paste("all_samples_survival",hallmark_lncRNA[j],sep="_")
name <- hallmark_lncRNA[j]
dif <- survdiff(Surv(train_time,train_status)~train_score[2,],train_score1)
    train_p <- 1-pchisq(dif[[5]],1)
	train_p <- round(train_p,4)
	lncRNA_risk[j,4] <- train_p
    pdf(file=paste(name,"_all_curve",".pdf",sep=""))
    plot(survfit(Surv(train_time2,train_status2)~train_score2[2,],train_score1),xscale=365.25,col=c('forestgreen','orangered'),lwd=2,xlab="train_time:days",main=name,ylab="Proportion survival")
    train_p1 <- paste("p-value=",round(train_p,4))
    text(16*0.7,0.9,train_p1,pos=4)
	train_low <- paste("low risk n=",as.numeric(table(train_score2[2,]))[1])
	train_high <- paste("high risk n=",as.numeric(table(train_score2[2,]))[2])
    legend(16*0.7,0.85,c(train_low,train_high),lty=1:1,col=c('forestgreen','orangered'))
    dev.off()

}

index1 <- which(as.numeric(lncRNA_risk[,2])<0.05)
index2 <- which(as.numeric(lncRNA_risk[,4])<0.05)
index <- intersect(index1,index2)
sig_lncRNA_risk <- lncRNA_risk[index,]
write.table(lncRNA_risk, "lncRNA_risk.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(sig_lncRNA_risk, "sig_lncRNA_risk.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
