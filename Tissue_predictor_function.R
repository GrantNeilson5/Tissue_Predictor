tissue_predictor <- function(QCmetrics, msetEPIC){
  
  load("/mnt/data1/Array_Projects/Tissue_Predictor/glmnet_class_predictor.rdat")
 
  if (exists("msetEPIC")){
  betas <- betas(msetEPIC)
  }else{
    betas <- data.matrix(betas)
  }
  
  ###Organising data frame
  ## Pulling out probes for testing 
  significant_Probes <- read.csv("/mnt/data1/Array_Projects/Tissue_Predictor/Probes_associated_with_Tissue.csv", stringsAsFactors = FALSE)
  ProbeIDs <- significant_Probes$SNP
  
 present <- ProbeIDs %in% row.names(betas)
  if(FALSE %in% present){
    stop("Not all prediction Probes Present")
    
  }else{
  message("All prediction probes present")
  ##Organising data frame
  prediction_betas <- betas[ProbeIDs,]
  prediction_betas <- t(prediction_betas)
  QCmetrics$Organ <- as.factor(QCmetrics$Organ)
  Tissue <- QCmetrics$Organ
  Tissue <- as.numeric(Tissue)
  
  ##Making Predictions
  predictions <- predict(alpha1.fit, s=alpha1.fit$lambda.min,  newx = prediction_betas, type= "response")
  predictions <- as.data.frame(predictions)
  
 tmp <-predict(alpha1.fit, s=alpha1.fit$lambda.min,  newx = prediction_betas, type= "class")
 predictions$Predicted_Tissue <- tmp[,1]
  ##Merging Predictions with QC metrics
  
  
  QCmetrics <- merge(QCmetrics, predictions, by.x = "Basename", by.y= "row.names")
  row.names(QCmetrics) <-QCmetrics$Basename
  
  QCmetrics$Mismatch_Tissue<-QCmetrics$Predicted_Tissue!= QCmetrics$Organ
  
  return(QCmetrics)
  } 
  
}


