library(doParallel)
library(caret)
library(glmnet)
library(ROCR)
cl<-makeCluster(12)
registerDoParallel(cl)

setwd("/mnt/data1/Array_Projects/Tissue_Predictor/")
load("EUGEI_ERisk_merged.rdat")

Master_Pheno$Tissue <-as.factor(Master_Pheno$Tissue)
Master_betas <- na.omit(Master_betas)
rownames(Master_Pheno) <- Master_Pheno$Basename

Master_Pheno <- Master_Pheno[order(rownames(Master_Pheno)),]
Master_betas <- Master_betas[,order(colnames(Master_betas))]
print(identical(rownames(Master_Pheno), colnames(Master_betas)))

significant_Probes <- read.csv("Probes_associated_with_Tissue.csv")
ProbeIDs <- significant_Probes$SNP
ProbeIDs <- as.character(ProbeIDs)

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


### ALpha1.fit is best predictor and will be saved
save(alpha1.fit, file ="glmnet_class_predictor.rdat")


predicted_rounded <- round(predicted1)


predicted_rounded <- gsub(1,"Blood",predicted_rounded)
predicted_rounded <- gsub(2,"Buccal",predicted_rounded)                   
predicted_rounded <- gsub(3,"Saliva",predicted_rounded)

validation_pheno$Tissue <- as.character(validation_pheno$Tissue)

test <-predicted_rounded == validation_pheno$Tissue

identical(predicted1, validation_pheno$Tissue)

save(alpha1.fit, file = "GLMnet_Tissue_predictor.rdat")

###Using Caret
## Formatting data

validation_betas <- testbetas[,-validation_index]
training_betas <- testbetas[,validation_index]


training_betas <- t(training_betas)
training_betas <- as.data.frame(training_betas)
training_betas$Tissue <- training_pheno$Tissue

validation_betas <- t(validation_betas)
validation_betas <- as.data.frame(validation_betas)
validation_betas$Tissue <- validation_pheno$Tissue

# Run algorithms using 10-fold cross validation
control <- trainControl(method="cv", number=10, savePredictions = TRUE, classProbs = TRUE)
metric <- "Accuracy"
##SVM
set.seed(7)
fit.svm_fewer_probes <- train(Tissue~., data=training_betas, method="svmRadial", metric=metric, trControl=control,  allowParallel=TRUE)

# Random Forest
set.seed(7)
fit.rf <- train(Tissue~., data=training_betas, method="rf", metric=metric, trControl=control,  allowParallel=TRUE)

#knn
set.seed(7)
fit.knn <- fit.rf <- train(Tissue~., data=training_betas, method="knn", metric=metric, trControl=control,  allowParallel=TRUE)


# summarize accuracy of models
results <- resamples(list(svm=fit.svm, svm_fewer_probes = fit.svm_fewer_probes))
summary(results)

pdf("Accuracy_of_model.pdf")
dotplot(results)
dev.off()


print(fit.svm)

##Validation

predictions <- predict(fit.svm, validation_betas)
confusionMatrix(predictions, validation_betas$Tissue)


#Checking variable importance for GBM

#Variable Importance
varImp(object=fit.svm)


##Save model
save(fit.svm, file = "Tissue_predictor.rdat")
prediction_probes <- ProbeIDs
write.csv(ProbeIDs, file = "Tissue_prediction_probes.csv",row.names = F)
load("Tissue_predictor.rdat")



