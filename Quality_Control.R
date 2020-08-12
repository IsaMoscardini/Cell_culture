#### Head ####
rm(list = ls())
options(stringsAsFactors = F)
library(dplyr)
library(Glimma)
library(edgeR)
library(ggplot2)
library(vegan)
library(mixOmics)
library(PCAtools)
library(DESeq2)

##### Analysis ####
##### Load tables 
desc <- read.csv("Data/VK2 table.csv")
count <- read.csv2("Data/VK2 all counts.csv")
#View(count)
#View(desc)

rownames(desc) <- desc$sample_name
rownames(count) <- count$Gene
count$Gene <- NULL

colnames(count) <- gsub("_.", "_", colnames(count))
colnames(count) <- gsub("UT", "PUT", colnames(count))

identical(colnames(count), rownames(desc)) # TRUE

# Creat DGE list
y <- DGEList(counts = count, genes = row.names(count), group= desc$treatment)
# filter
keep <- rowSums(cpm(y)>2)>=4
y.1 <- y[keep,]
dim(y.1) # 11741    28

# normalize
y.1 <- calcNormFactors(y.1) # normaliza as reads ver vst

cp <- log2(cpm(y.1)+1) # normalizou de vdd
sample.var <- apply(cp, 2, var)
sample.var <- sort(sample.var)
barplot(sample.var, las=2)    # NEEDS SCALATION! PLEASE ATTENTION! - Do quality control if it is possible

# PCA tools ####


# PCA ####
#### All samples
pca<-prcomp(t(cp))
pcapanel<- as.data.frame(pca$x)
percentVar <- pca$sdev^2/sum(pca$sdev^2)

percentVar<- round(100* percentVar)

group.colors <- c(NI_NS = "#66C2A5", NI_T = "#8DA0CB", I_NS ="#FC8D62", I_T = "#E78AC3")

ggplot(pcapanel, aes(PC1, PC2, color=desc$treatment))+
  geom_point(size = 5)+ theme_minimal()+ theme(legend.position = "top")+ labs(group = "Treatment")+
  #geom_text(label=desc$sample_name, nudge_x = 0.9, nudge_y = 1)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

 #### No drugs
desc_nd <- desc[desc$treatment %in% c("CTR","CAD_PUT"),]
count_nd <- count[, intersect(rownames(desc_nd), colnames(count))]
identical(colnames(count_nd), rownames(desc_nd))

y <- DGEList(counts = count_nd, genes = row.names(count_nd), group= desc_nd$treatment)
keep <- rowSums(cpm(y)>2)>=4
y.1 <- y[keep,]
dim(y.1) # 11741    28
y.1 <- calcNormFactors(y.1) # normaliza as reads ver vst

cp <- log2(cpm(y.1)+1) # normalizou de vdd
pca<-prcomp(t(cp))
pcapanel<- as.data.frame(pca$x)
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar<- round(100* percentVar)

ggplot(pcapanel, aes(PC2, PC3, color=desc_nd$treatment))+
  geom_point(size = 5)+ theme_minimal()+ theme(legend.position = "top")+ labs(group = "Treatment")+
  xlab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylab(paste0("PC3: ",percentVar[3],"% variance"))

ipca.res <- ipca(t(cp), ncomp = 3)
plotIndiv(ipca.res, comp = c(1,2), group = desc_nd$treatment, title = "IPCA comp 1 and 2")
View(ipca.res$rotation)


# IPCA ####
ipca.res <- ipca(t(cp), ncomp = 4)
barplot(ipca.res$kurtosis)
plotIndiv(ipca.res, comp = c(1, 2), group = desc$treatment)

sipca.res <- sipca(t(cp), ncomp = 4)
plotIndiv(sipca.res, comp = c(1, 2), group = desc$treatment)
View(sipca.res$rotation)

