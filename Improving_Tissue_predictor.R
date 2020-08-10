### validating tissue predictor 
library(glmnet)
library(caret)
setwd("/mnt/data1/Array_Projects/Tissue_Predictor/")
load("/mnt/data1/reference_files/Tissue_Predictor/glmnet_class_predictor.rdat")

##Getting lsit of coeffs for each tissue 
tmp_coeffs <- coef(alpha1.fit, s = "lambda.min")
Blood_coef<- data.frame(name = tmp_coeffs$Blood@Dimnames[[1]][tmp_coeffs$Blood@i+1], coefficient = tmp_coeffs$Blood@x)
Blood_coef <- Blood_coef[-1,]
Buccal_coef <- data.frame(name = tmp_coeffs$Buccal@Dimnames[[1]][tmp_coeffs$Buccal@i+1], coefficient = tmp_coeffs$Buccal@x)
Buccal_coef <- Buccal_coef[-1,]
Saliva_coef<-data.frame(name = tmp_coeffs$Saliva@Dimnames[[1]][tmp_coeffs$Saliva@i+1], coefficient = tmp_coeffs$Saliva@x)
Saliva_coef <- Saliva_coef[-1,]
#Merging and saving  all coeffs
all_coefs <- rbind(Blood_coef, Buccal_coef, Saliva_coef)

write.csv(all_coefs, "new_coeffs.csv", row.names = F)

load("EUGEI_ERisk_merged_EPIC_Probes.rdat")

Master_Pheno$Tissue <-as.factor(Master_Pheno$Tissue)
Master_betas <- na.omit(Master_betas)
rownames(Master_Pheno) <- Master_Pheno$Basename

Master_Pheno <- Master_Pheno[order(rownames(Master_Pheno)),]
Master_betas <- Master_betas[,order(colnames(Master_betas))]
print(identical(rownames(Master_Pheno), colnames(Master_betas)))

Probes <- read.csv("new_coeffs.csv", stringsAsFactors = F)
ProbeIDs <- Probes$name


testbetas <- Master_betas[ProbeIDs,]

## Making training and validation dataset

validation_index <- createDataPartition(Master_Pheno$Tissue, p=0.80, list=FALSE)
validation_pheno <- Master_Pheno[-validation_index,]
# use the remaining 80% of data to training and testing the models
training_pheno <- Master_Pheno[validation_index,]

validation_betas <- testbetas[,-validation_index]
training_betas <- testbetas[,validation_index]

training_betas <- t(training_betas)


validation_betas <- t(validation_betas)



## USing GLMnet, formatting the data correctly 
Training_Tissue <- training_pheno$Tissue
#Training_Tissue <- as.numeric(Training_Tissue)

validation_Tissue <- validation_pheno$Tissue
#validation_Tissue <- as.numeric(validation_Tissue)



#using GLMnet
## Fitting model
alpha0.fit = cv.glmnet(training_betas, Training_Tissue, type.measure = "class", alpha =0,  family = "multinomial")

alpha1.fit = cv.glmnet(training_betas, Training_Tissue, type.measure = "class", alpha =1,  family = "multinomial")

##testing model accuracy
predicted <- predict(alpha0.fit, s=alpha0.fit$lambda.min,  newx = validation_betas, type = "class")
predicted1 <- predict(alpha1.fit, s=alpha1.fit$lambda.min,  newx = validation_betas, type = "class" )


cm <- confusionMatrix(as.factor(predicted[,1]), as.factor(validation_Tissue))
cm1 <- confusionMatrix(as.factor(predicted1[,1]), as.factor(validation_Tissue))


### ALpha0.fit is best predictor and will be saved
save(alpha0.fit, file ="glmnet_class_predictor_new.rdat")
