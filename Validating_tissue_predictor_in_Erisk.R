### validating tissue predictor 
library(glmnet)
library(caret)
setwd("/mnt/data1/Array_Projects/Tissue_Predictor/")

load("/mnt/data1/Eilis/Projects/Asthma/CellTypeComparisons/Correlations_New/BetasSortedByCellType_NoY.rdat")
source("/mnt/data1/reference_files/Tissue_Predictor/Tissue_predictor_function.R")

betas <- allbetas[["Buccal"]]

Predicted_Tissue <-tissue_predictor(allbetas[["Buccal"]])
Predicted_Tissue <- rbind(Predicted_Tissue,tissue_predictor(allbetas[["B-cells"]]))
Predicted_Tissue <- rbind(Predicted_Tissue,tissue_predictor(allbetas[["CD4 T-cells"]]))
Predicted_Tissue <- rbind(Predicted_Tissue,tissue_predictor(allbetas[["CD8 T-cells"]]))
Predicted_Tissue <- rbind(Predicted_Tissue,tissue_predictor(allbetas[["Monocytes"]]))
Predicted_Tissue <- rbind(Predicted_Tissue,tissue_predictor(allbetas[["whole blood"]]))
Predicted_Tissue <- rbind(Predicted_Tissue,tissue_predictor(allbetas[["Granulocytes"]]))
Predicted_Tissue <- rbind(Predicted_Tissue,tissue_predictor(allbetas[["Nasal"]]))



Predicted_Tissue$Reported_Tissue <-c(rep("buccal", 30),rep("B_cells",30),rep("CD4",30),rep("CD8",30),rep("Mono",30),
                     rep("Whole_Blood",30), rep("Gran",30), rep("Nasal",30))

png("Validation in Erisk.png")
ggplot(data = Predicted_Tissue, aes(x=Reported_Tissue, fill= Predicted_Tissue))+
  geom_bar( position= position_dodge())+
  labs(x = "Reported Tissue", fill= "Predicted Tissue")+
  theme(axis.text.x = element_text(angle = 90))
dev.off() 
