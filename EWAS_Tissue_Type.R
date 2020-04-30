#Setting up parallel processors
library(doParallel)
library(lme4)
library(qqman)
library(gplots)
cl<-makeCluster(16)
registerDoParallel(cl)
clusterEvalQ(cl, library(lme4))
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)

setwd("/mnt/data1/Array_Projects/Tissue_Predictor/")

load("EUGEI_ERisk_merged.rdat")

rownames(Master_Pheno) <- Master_Pheno$Basename
Master_Pheno <- Master_Pheno[order(rownames(Master_Pheno)),]
Master_betas <- Master_betas[,order(colnames(Master_betas))]
print(identical(rownames(Master_Pheno), colnames(Master_betas)))

# splitting basenmae into chip and postion as Excel often mucks up Chip
Master_betas <- as.matrix(Master_betas)
Master_Pheno$Basename2<-Master_Pheno$Basename
Master_Pheno<-separate(data = Master_Pheno, col = Basename2, into = c("Chip", "Position"), sep="_")


Master_Pheno$Plate <- as.factor(Master_Pheno$Plate)
Master_Pheno$Sex <- as.factor(Master_Pheno$Sex)
Master_Pheno$Tissue <-as.factor(Master_Pheno$Tissue)
Master_Pheno$Phenotype <-as.factor(Master_Pheno$Phenotype)
Master_Pheno$Age <-as.numeric(as.character(Master_Pheno$Age))

Master_betas <- na.omit(Master_betas)

testCpG<-function(row, Master_Pheno){
  
  model<-lm(Master_betas[i,] ~ Tissue + Age + Sex + Plate +Phenotype, data = Master_Pheno, REML = FALSE)
  Beta<-coefficients(model)["TissueSaliva"]
  SE<-summary(model)$coefficients["TissueSaliva", 2]
  P<-summary(model)$coefficients["TissueSaliva", 4]
  return(c(Beta,SE,P))
}

## Running function and outputing file intp res 

res<-foreach(i= 1:nrow(Master_betas), .combine=rbind) %dopar%{
  testCpG(Master_betas[i,], Master_Pheno)	
}

for (i in 1:nrow(Master_betas)) {
  res <- testCpG(Master_betas[i,], Master_Pheno)
}

rownames(res)<-rownames(Master_betas)
colnames(res)<-c("Tissue_Beta", "Tissue_SE", "Tissue_P")

# As the beta values here represent proportion of DNA methylation (i.e. they lie between 0 and 1), the regression coefficients represent the change in proportion. 
# Typically we report our findings on the % scale therefore we will multiple the regression coefficients and SE by 100. 
# This will need to be done for all variables for which you saved the results.


res[,"Tissue_Beta"]<-res[,"Tissue_Beta"]*100
res[,"Tissue_SE"]<-res[,"Tissue_SE"]*100

### Annotate the output

# At this stage we will be interested in adding information regarding where each probe is located, and what genes or regulatory features it overlaps with. 
# There is a lot of annotation available, the code below only takes a subset of columns which I think may be most relevant. 


epicManifest<-read.csv("/gpfs/mrc0/projects/Research_Project-MRC193462/EWAS/epicManifest_hg38.csv", stringsAsFactors = F, header=T, sep=" ")


epicManifest<-epicManifest[match(rownames(res), epicManifest$Name),]
res<-cbind(res, epicManifest[,c("CHR", "MAPINFO")]) ## this can be edited to include additional columns or exclude as required

write.csv(res, file="EWAS_Tissue_Type.csv")

res <- read.csv("EWAS_Tissue_Type.csv", row.names = 1)

## QQ plot 
pdf("QQplot_of_EWAS_Age.pdf")
qq(res$Tissue_P)
dev.off()

##Manhattan Plot
res<-res[which(res$CHR != "Y"),] ## to also exclude the X chromosome repeat and edit this line
res<-res[which(res$CHR != "X"),] ## to also exclude the X chromosome repeat and edit this line


#res1$CHR<-as.character(res$CHR)
#res$CHR[which(res$CHR == "X")]<-23 ## to also recode Y chromosome repeat and edit this line
#res$CHR[which(res$CHR == "Y")]<-23 ## to also recode Y chromosome repeat and edit this line
res$CHR<-as.numeric(as.character(res$CHR))
res<-res[which(res$CHR != ""),]
res<-res[which(res$MAPINFO != ""),]

#bonfP<-0.05/nrow(res)
bonfP<-0.05/nrow(res)

pdf("Manhattan_plot_of_Tissue.pdf")
manhattan(res, p = "Tissue_P", bp = "MAPINFO", chr = "CHR", genomewide = -log10(bonfP), suggestiveline = -log10(5e-5), logp=T, col=c("blue","yellow"), ylim=c(0,50))
dev.off()

significant_SNPs <- res[which(res$Tissue_P < bonfP),]

write.csv(significant_SNPs, file = "Probes_associated_with_Tissue.csv", row.names = TRUE)
