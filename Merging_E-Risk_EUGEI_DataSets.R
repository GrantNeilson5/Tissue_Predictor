library(dplyr)
moveme <- function (invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], 
                                 ",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", 
                              "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      }
      else if (A == "after") {
        after <- match(ba, temp)
      }
    }
    else if (A == "first") {
      after <- 0
    }
    else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}

setwd("/mnt/data1/Array_Projects/Tissue_Predictor/")

datMiniAnnotation=read.csv("/mnt/data1/Array_Projects/Portland_Samples/Hovarth_Tissue_Predictor/datMiniAnnotation.csv")
###########################
######Loading EUGEI data###
###########################

########################################
##Formatting EUGEI Blood Spread sheet###
########################################


load("/mnt/data1/EuGEI/QC/GeorginasQC/All_Plates_Blood_WithRepeats/EuGEIBloodSamples_Normalised.rdat")
pheno$Tissue <- rep("Blood", nrow(pheno))
EUGI_pheno_blood <- pheno
EUGI_betas_blood <- betas


EUGI_pheno_blood1 <- EUGI_pheno_blood[,c("Basename","Eilis.Sample_Name","MethPlate","Tissue","Sex","Age","Phenotype.EuGEI")]

EUGI_betas_blood1 <- as.data.frame(EUGI_betas_blood)

EUGI_betas_blood1$ProbeID <- rownames(EUGI_betas_blood1)

EUGI_betas_blood1 <- EUGI_betas_blood1[moveme(names(EUGI_betas_blood1), "ProbeID first")]

#Making the beta matrix into 27k probes

match1=match(datMiniAnnotation[,1], EUGI_betas_blood1[,1] )
EUGI_betas_blood1Reduced=EUGI_betas_blood1[match1,]
EUGI_betas_blood1Reduced[,1]=as.character(EUGI_betas_blood1Reduced[,1])
EUGI_betas_blood1Reduced[is.na(match1),1]= as.character(datMiniAnnotation[is.na(match1),1])
EUGI_betas_blood1= EUGI_betas_blood1Reduced

#Removing Probe ID column 
EUGI_betas_blood1 <- EUGI_betas_blood1[,-1]

#########################################
##Formatting EUGEI Saliva Spread sheet###
#########################################


load("/mnt/data1/EuGEI/QC/GeorginasQC/All_Plates_Buccal/EuGEIBuccalSamples_Normalised.rdat")
colnames(pheno)[9] <- "Tissue"

EUGI_pheno_Saliva <-pheno
EUGI_betas_Saliva <- betas

EUGI_pheno_Saliva1 <- EUGI_pheno_Saliva[,c("Basename","Eilis.Sample_Name","MethPlate","Tissue","Sex","Age","Phenotype.EuGEI")]

EUGI_betas_Saliva1 <- as.data.frame(EUGI_betas_Saliva)

EUGI_betas_Saliva1$ProbeID <- rownames(EUGI_betas_Saliva1)

EUGI_betas_Saliva1 <- EUGI_betas_Saliva1[moveme(names(EUGI_betas_Saliva1), "ProbeID first")]


#Making the beta matrix into 27k probes

match1=match(datMiniAnnotation[,1], EUGI_betas_Saliva1[,1] )
EUGI_betas_Saliva1Reduced=EUGI_betas_Saliva1[match1,]
EUGI_betas_Saliva1Reduced[,1]=as.character(EUGI_betas_Saliva1Reduced[,1])
EUGI_betas_Saliva1Reduced[is.na(match1),1]= as.character(datMiniAnnotation[is.na(match1),1])
EUGI_betas_Saliva1= EUGI_betas_Saliva1Reduced

#Removing Probe ID column 
EUGI_betas_Saliva1 <- EUGI_betas_Saliva1[,-1]



##Mergeing EUGEI blood and saliva sheets
EUGEI_pheno <- rbind(EUGI_pheno_blood1, EUGI_pheno_Saliva1)

##Keeping Epic array probes
Matching_cpgs <- match(EUGI_betas_Saliva1[,1], EUGI_betas_Saliva1[,1])
EUGI_betas_BloodReduced=EUGI_betas_blood1[Matching_cpgs,]
EUGI_betas_Saliva1 <- EUGI_betas_Saliva1[,-1]
EUGI_betas_BloodReduced <- EUGI_betas_BloodReduced[,-1]

EUGEI_betas <- cbind(EUGI_betas_BloodReduced, EUGI_betas_Saliva1)

colnames(EUGEI_pheno)[c(2,3,7)]<- c("SampleID", "Plate", "Phenotype")

########################
###loading Erisk data### 
########################

load("/mnt/data1/E-risk/LongitudinalTwinQC_Nov17/QC/AllPlates/LongitudinalERisk_Normalised_28022018.rdat")

ERisk_Blood_pheno <- pheno[which(pheno$Tissue =="Blood"),]

ERisk_buccal_pheno <- pheno[which(pheno$Tissue =="Buccal" & pheno$Age == "Age 18"),]
ERisk_pheno <- rbind(ERisk_Blood_pheno,ERisk_buccal_pheno)

##Making ERisk Betas sheet

ERisk_betas <- betas[,ERisk_pheno$Basename]

ERisk_betas1 <- as.data.frame(ERisk_betas)

ERisk_betas1$ProbeID <- rownames(ERisk_betas1)

ERisk_betas1 <- ERisk_betas1[moveme(names(ERisk_betas1), "ProbeID first")]


#Giving the beta matrix the correct Number of probes

match1=match(datMiniAnnotation[,1], ERisk_betas1[,1] )
ERisk_betas1Reduced=ERisk_betas1[match1,]
ERisk_betas1Reduced[,1]=as.character(ERisk_betas1Reduced[,1])
ERisk_betas1Reduced[is.na(match1),1]= as.character(datMiniAnnotation[is.na(match1),1])
ERisk_betas1= ERisk_betas1Reduced

#Removing Probe ID column 
ERisk_betas1 <- ERisk_betas1[,-1]



## Making ERisk spread sheet same as EUGEI

ERisk_pheno <- ERisk_pheno[,c("Basename","SampleID","Tissue","SEX","Age", "Plate.name")]
colnames(ERisk_pheno)[c(4,6)] <- c("Sex","Plate")
ERisk_pheno$Phenotype <- rep(NA, nrow(ERisk_pheno))
ERisk_pheno <- ERisk_pheno [moveme(names(ERisk_pheno), "Plate before Tissue")]


#######################
##Merging Everything###
#######################

Master_Pheno <- rbind(ERisk_pheno,EUGEI_pheno)

##Keeping Epic array probes

EUGEI_betas$ProbeID <- rownames(EUGEI_betas)

EUGEI_betas <- EUGEI_betas[moveme(names(EUGEI_betas), "ProbeID first")]

Matching_cpgs <- match(ERisk_betas1[,1], EUGEI_betas[,1])
EUGI_betas_Reduced=EUGEI_betas[Matching_cpgs,]
ERisk_betas1 <- ERisk_betas1[,-1]
EUGI_betas_Reduced <- EUGI_betas_Reduced[,-1]



Master_betas <- cbind(ERisk_betas1, EUGI_betas_Reduced)

row.names(Master_Pheno) <- Master_Pheno$Basename

Master_Pheno <- Master_Pheno[order(rownames(Master_Pheno)),]
Master_betas <- Master_betas[,order(colnames(Master_betas))]

print(identical(rownames(Master_Pheno), colnames(Master_betas)))

save(Master_betas, Master_Pheno, file = "EUGEI_ERisk_merged_EPIC_Probes.rdat")
