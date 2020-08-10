tissue_predictor <- function(betas){
  
  #load("/mnt/data1/reference_files/Tissue_Predictor/glmnet_class_predictor")
  load("/mnt/data1/Array_Projects/Tissue_Predictor/glmnet_class_predictor_new.rdat")
  
  betas <- data.matrix(betas)
  
  
  ###Organising data frame
  ## Pulling out probes for testing 
  #significant_Probes <- read.csv("/mnt/data1/reference_files/Tissue_Predictor/Probes_associated_with_Tissue.csv", stringsAsFactors = FALSE)
  #ProbeIDs <- significant_Probes$SNP
  
  Probes <- read.csv("/mnt/data1/Array_Projects/Tissue_Predictor/new_coeffs.csv", stringsAsFactors = F)
  ProbeIDs <- Probes$name

  present <- ProbeIDs %in% row.names(betas)
  if(FALSE %in% present){
    stop("Not all prediction Probes Present")
    
  }else{
    message("All prediction probes present")
    ##Organising data frame
    prediction_betas <- betas[ProbeIDs,]
    prediction_betas <- t(prediction_betas)
    
    
    ##Making Predictions
    predictions <- predict(alpha0.fit, s=alpha0.fit$lambda.min,  newx = prediction_betas, type= "response")
    predictions <- predictions[1:nrow(predictions),1:ncol(predictions),]
    predictions <- as.data.frame(predictions)
   
    
    
    tmp <-predict(alpha0.fit, s=alpha0.fit$lambda.min,  newx = prediction_betas, type= "class")
    predictions$Predicted_Tissue <- tmp[,1]
    
    return(predictions)
  } 
  
}


