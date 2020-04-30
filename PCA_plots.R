##PCAs
library(ggplot2)
library(tidyverse)
library(cowplot)
library(rgl)
library(car)
library(magick)
library(magrittr)

setwd("/mnt/data1/Array_Projects/Tissue_Predictor/")
load("EUGEI_ERisk_merged.rdat")

## Generating PCA's
Master_betas <- na.omit(Master_betas)
Master_betas <-t(Master_betas)

Master_betas.scaled <- data.frame(apply(Master_betas,2,scale))

intensity_pca <- prcomp(Master_betas, center = TRUE, scale = TRUE)

####organising data for plotting
PCA <- intensity_pca$x
PCA <- as.data.frame(PCA)
PCA1 <- select(PCA, PC1, PC2,PC3)
PCA1 <- PCA1[order(rownames(PCA1)),]
print(identical(rownames(Master_Pheno),rownames(PCA1)))
Tissue <- Master_Pheno$Tissue
PCA1 <- cbind(PCA1, Tissue)

#### plotting graphs
pdf("PCA_of_EUGEI_Erisk_betas.pdf")
ggplot(PCA1, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = factor(Master_Pheno$Tissue)))+
  labs(
    title = "PCAs of Merged Betas from EUGEI and E-risk"
  )
dev.off()




with(PCA1, scatter3d(PC1, PC2, PC3, groups = Tissue, grid =FALSE, surface = FALSE, ellipsoid = TRUE, 
                     xlab ="PC1", ylab = "PC2", zlab = "PC3",
                     axis.col = c("black", "black", "black")))

palette(rainbow(3))
with(PCA1, plot3d(PC1,PC2,PC3, col = as.integer(Tissue),type = "s", radius = 7, 
                  xlab ="PC1", ylab = "PC2", zlab = "PC3"))
legend3d("topright", legend = levels(PCA1$Tissue), col = c("red", "green", "blue"), pch=19)

# Save like gif
setwd("/mnt/data1/Array_Projects/Tissue_Predictor/Plots/")
movie3d(
  movie="pca_plot",
  spin3d( axis = c(0, 0, 1), rpm = 3),
  duration = 10,
  dir = "./",
  type = "gif",
  clean = TRUE
)





