library(tidyverse)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationHub)
library(GOSemSim)
library(DOSE)
library(clusterProfiler)
library(msigdbr)
library(ggplot2)
library(gridExtra)
library(stringr)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")
BiocManager::install("clusterProfiler")
BiocManager::install("AnnotationHub", force = TRUE)

data <- read_csv(file = "PPIs_all.csv")
#######计算CC/MF correlation
hsGO1 <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont="CC", computeIC=FALSE) 
hsGO2 <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont="MF", computeIC=FALSE) 
hsGO3 <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont="BP", computeIC=FALSE) 

genesim_CC <- as.vector(NULL)
genesim_MF <- as.vector(NULL)
genesim_BP <- as.vector(NULL)

for (i in 1:nrow(data)) {
  sim1 <- GOSemSim::geneSim(as.character(data[i,5]),as.character(data[i,6]), semData=hsGO1, measure="Wang", drop=NULL, combine="BMA")
  sim2 <- GOSemSim::geneSim(as.character(data[i,5]),as.character(data[i,6]), semData=hsGO2, measure="Wang", drop=NULL, combine="BMA")
  sim3 <- GOSemSim::geneSim(as.character(data[i,5]),as.character(data[i,6]), semData=hsGO3, measure="Wang", drop=NULL, combine="BMA")
  genesim_CC <- append(genesim_CC, sim1[[1]])
  genesim_MF <- append(genesim_MF, sim2[[1]])
  genesim_BP <- append(genesim_BP, sim3[[1]])
}
genesim_CC <- as.matrix(genesim_CC)
genesim_MF <- as.matrix(genesim_MF)
genesim_BP <- as.matrix(genesim_BP)

colnames(genesim_CC) <- ("CC")
colnames(genesim_MF) <- ("MF")
colnames(genesim_BP) <- ("BP")

data <- cbind(data, genesim_CC,genesim_MF,genesim_BP)
data[is.na(data)] <- 0

write.csv(data, file = "PPI_all.csv")