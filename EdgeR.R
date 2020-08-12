#### Head
rm(list = ls())
options(stringsAsFactors = F)

sapply(c("Glimma", "edgeR", "ggplot2", "dplyr"), require, character.only=T)

####
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

Group <- factor(desc$treatment)
cbind(desc, Group=Group)
design <- model.matrix(~0+Group)
colnames(design) <- levels(Group)
View(design)
fit <- estimateDisp(y.1, design)
qlfit <- glmQLFit(fit, design)

my.contrasts <- makeContrasts(
  CAD_PUT = CAD_PUT-CTR,
  DPV = DPV-CTR,
  DRV = DRV-CTR,
  CAD_PUT_DPV = CAD_PUT_DPV-CTR,
  CAD_PUT_DRV = CAD_PUT_DRV-CTR,
  levels=design)
View(my.contrasts)

dir.create("Results/DEGs/EdgeR", recursive = T)

# CAD_PUT
ql_CAD_PUT <- glmQLFTest(qlfit, contrast=my.contrasts[,"CAD_PUT"])
CAD_PUT <- topTags(ql_CAD_PUT, n=Inf)
write.csv2(CAD_PUT, "Results/DEGs/EdgeR/CAD_PUT.csv")
summary (dup<- decideTestsDGE(ql_CAD_PUT))
# 1*CAD_PUT -1*CTR
# Down                 32
# NotSig            11672
# Up                   37

# DPV
ql_DPV <- glmQLFTest(qlfit, contrast=my.contrasts[,"DPV"])
DPV <- topTags(ql_DPV, n=Inf)
write.csv2(DPV, "Results/DEGs/EdgeR/DPV.csv")
summary (dup<- decideTestsDGE(ql_DPV))
# -1*CTR 1*DPV
# Down           3308
# NotSig         5000
# Up             3433

# DRV
ql_DRV <- glmQLFTest(qlfit, contrast=my.contrasts[,"DRV"])
DRV <- topTags(ql_DRV, n=Inf)
write.csv2(DRV, "Results/DEGs/EdgeR/DRV.csv")
summary (dup<- decideTestsDGE(ql_DRV))
# -1*CTR 1*DRV
# Down           2811
# NotSig         6043
# Up             2887

# CAD_PUT_DPV
ql_CAD_PUT_DPV <- glmQLFTest(qlfit, contrast=my.contrasts[,"CAD_PUT_DPV"])
CAD_PUT_DPV <- topTags(ql_CAD_PUT_DPV, n=Inf)
write.csv2(CAD_PUT_DPV, "Results/DEGs/EdgeR/CAD_PUT_DPV.csv")
summary (dup<- decideTestsDGE(ql_CAD_PUT_DPV))
# 1*CAD_PUT_DPV -1*CTR
# Down                   3429
# NotSig                 4802
# Up                     3510

# CAD_PUT_DRV
ql_CAD_PUT_DRV <- glmQLFTest(qlfit, contrast=my.contrasts[,"CAD_PUT_DRV"])
CAD_PUT_DRV <- topTags(ql_CAD_PUT_DRV, n=Inf)
write.csv2(CAD_PUT_DRV, "Results/DEGs/EdgeR/CAD_PUT_DRV.csv")
summary (dup<- decideTestsDGE(ql_CAD_PUT_DRV))
# 1*CAD_PUT_DRV -1*CTR
# Down                   3056
# NotSig                 5726
# Up                     2959
























