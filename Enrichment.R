## Header #####
rm(list = ls())
options(stringsAsFactors = F)

sapply(c("Glimma", "edgeR", "ggplot2", "dplyr", "tmod"), require, character.only=T)


cad_put <- read.csv2("Results/DEGs/EdgeR/CAD_PUT.csv")
dpv <- read.csv2("Results/DEGs/EdgeR/DPV.csv")
drv <- read.csv2("Results/DEGs/EdgeR/DRV.csv")
cp_dpv <- read.csv2("Results/DEGs/EdgeR/CAD_PUT_DPV.csv")
cp_drv <- read.csv2("Results/DEGs/EdgeR/CAD_PUT_DRV.csv")
bg <- read.csv2("Data/VK2 all counts.csv")
bg <- bg$Gene

## Analysis Reactome ####

msig=tmodImportMSigDB(file = "Reactome_2016.gmt", format = "gmt")
C2 <- msig$MODULES$ID
View(C2)

dir.create("Results/Reactome")

# cad_put_LI <- tmodCERNOtest(cad_put$genes, mset= msig[C2])
# View(cad_put_LI)  # AQUI NAO TEM NADA
# write.csv2(cad_put_LI, "Results/Reactome/cad_put_RE.csv", row.names = F)  

##
dpv_up <- dpv[which(dpv$FDR < 0.05 & dpv$logFC > 0.322),]
dpv_up_re <- tmodHGtest(fg = dpv_up$genes, bg = bg, mset = msig[C2])
View(dpv_up_re)  
write.csv2(dpv_up_re, "Results/Reactome/dpv_up_re.csv", row.names = F)

dpv_down <- dpv[which(dpv$FDR < 0.05 & dpv$logFC < -0.322),]
dpv_down_re <- tmodHGtest(fg = dpv_down$genes, bg = bg, mset = msig[C2])
View(dpv_down_re)  
write.csv2(dpv_down_re, "Results/Reactome/dpv_down_re.csv", row.names = F)

##
drv_up <- drv[which(drv$FDR < 0.05 & drv$logFC > 0.322),]
drv_up_re <- tmodHGtest(fg = drv_up$genes, bg = bg, mset = msig[C2])
View(drv_up_re)  
write.csv2(drv_up_re, "Results/Reactome/drv_up_re.csv", row.names = F)

drv_down <- drv[which(drv$FDR < 0.05 & drv$logFC < -0.322),]
drv_down_re <- tmodHGtest(fg = drv_down$genes, bg = bg, mset = msig[C2])
View(drv_down_re)  
write.csv2(drv_down_re, "Results/Reactome/drv_down_re.csv", row.names = F)

##
cp_dpv_up <- cp_dpv[which(cp_dpv$FDR < 0.05 & cp_dpv$logFC > 0.322),]
cp_dpv_up <- tmodHGtest(fg = cp_dpv_up$genes, bg = bg, mset = msig[C2])
View(cp_dpv_up)  
write.csv2(cp_dpv_up, "Results/Reactome/cp_dpv_up.csv", row.names = F)

cp_dpv_down <- cp_dpv[which(cp_dpv$FDR < 0.05 & cp_dpv$logFC < -0.322),]
cp_dpv_down <- tmodHGtest(fg = cp_dpv_down$genes, bg = bg, mset = msig[C2])
View(cp_dpv_down)
write.csv2(cp_dpv_down, "Results/Reactome/cp_dpv_down.csv", row.names = F)

##
cp_drv_up <- cp_drv[which(cp_drv$FDR < 0.05 & cp_drv$logFC > 0.322),]
cp_drv_up <- tmodHGtest(fg = cp_drv_up$genes, bg = bg, mset = msig[C2])
View(cp_drv_up)  
write.csv2(cp_drv_up, "Results/Reactome/cp_drv_up.csv", row.names = F)

cp_drv_down <- cp_drv[which(cp_drv$FDR < 0.05 & cp_drv$logFC < -0.322),]
cp_drv_down <- tmodHGtest(fg = cp_drv_down$genes, bg = bg, mset = msig[C2])
View(cp_drv_down)
write.csv2(cp_drv_down, "Results/Reactome/cp_drv_down.csv", row.names = F)
  
## Graphs ####




## Analysis tmod ####

cad_put <- cad_put[order(cad_put$FDR), ]
View(cad_put)
dpv <- dpv[order(dpv$FDR), ]
drv <- drv[order(drv$FDR), ]
cp_dpv <- cp_dpv[order(cp_dpv$FDR), ]
cp_drv <- cp_drv[order(cp_drv$FDR), ]
View(cad_put)

dir.create("Results/tmod")

cad_put_LI <- tmodCERNOtest(cad_put$genes)
View(cad_put_LI)  
write.csv2(cad_put_LI, "Results/tmod/cad_put_LI.csv", row.names = F)  

dpv_LI <- tmodCERNOtest(dpv$genes)
View(dpv_LI)  
write.csv2(dpv_LI, "Results/tmod/dpv_LI.csv", row.names = F)

drv_LI <- tmodCERNOtest(drv$genes)
View(drv_LI)  
write.csv2(drv_LI, "Results/tmod/drv_LI.csv", row.names = F)

cp_dpv_LI <- tmodCERNOtest(cp_dpv$genes)
View(cp_dpv_LI)  
write.csv2(cp_dpv_LI, "Results/tmod/cp_dpv_LI.csv", row.names = F)

cp_drv_LI <- tmodCERNOtest(cp_drv$genes)
View(cp_drv_LI)  
write.csv2(cp_drv_LI, "Results/tmod/cp_drv_LI.csv", row.names = F)  

## Graphs ####
View(cad_put)

pie_cp <- tmodDecideTests(cad_put$genes, lfc=cad_put$logFC, pval=cad_put$FDR, 
                         mset = "LI", pval.thr = 0.05, lfc.thr = 0)
pie_dpv <- tmodDecideTests(dpv$genes, lfc=dpv$logFC, pval=dpv$FDR, 
                          mset = "LI", pval.thr = 0.05, lfc.thr = 0)
pie_drv <- tmodDecideTests(drv$genes, lfc=drv$logFC, pval=drv$FDR, 
                          mset = "LI", pval.thr = 0.05, lfc.thr = 0)
pie_cp_dpv <- tmodDecideTests(cp_dpv$genes, lfc=cp_dpv$logFC, pval=cp_dpv$FDR, 
                          mset = "LI", pval.thr = 0.05, lfc.thr = 0)
pie_cp_drv <- tmodDecideTests(cp_drv$genes, lfc=cp_drv$logFC, pval=cp_drv$FDR, 
                          mset = "LI", pval.thr = 0.05, lfc.thr = 0)

# transform in data frame
pie_cp <- as.data.frame(pie_cp)
pie_dpv <- as.data.frame(pie_dpv)
pie_drv <- as.data.frame(pie_drv)
pie_cp_dpv <- as.data.frame(pie_cp_dpv)
pie_cp_drv <- as.data.frame(pie_cp_drv)

# Remove the X if there is any
colnames(pie_cp) <- gsub("X.*\\.", "", colnames(pie_cp))
colnames(pie_dpv) <- gsub("X.*\\.", "", colnames(pie_dpv))
colnames(pie_drv) <- gsub("X.*\\.", "", colnames(pie_drv))
colnames(pie_cp_dpv) <- gsub("X.*\\.", "", colnames(pie_cp_dpv))
colnames(pie_cp_drv) <- gsub("X.*\\.", "", colnames(pie_cp_drv))

##
pie_cp <- list("Cad_Put"=pie_cp)
pie_dpv <- list("DPV"=pie_dpv)
pie_drv <- list("DRV"=pie_drv)
pie_cp_dpv <- list("CP_DPV"=pie_cp_dpv)
pie_cp_drv <- list("CP_DRV"=pie_cp_drv)

##
CP_LI <- list("Cad_Put"=cad_put_LI)
DPV_LI <- list("DPV"=dpv_LI)
DRV_LI <- list("DRV"=drv_LI)
CP_DPV_LI <- list("CP_DPV"=cp_dpv_LI)
CP_DRV_LI <- list("CP_DRV"=cp_drv_LI)

# Create the panels
# I_NS NI_NS vs I_T NI_NS
panel_list <- c(CP_LI, DPV_LI, CP_DPV_LI, DRV_LI, CP_DRV_LI)
is.list(panel_list)
pie_list <- c(pie_cp, pie_dpv, pie_cp_dpv, pie_drv, pie_cp_drv)

x11()
tmodPanelPlot(panel_list, pval.thr = 0.05, 
              pval.cutoff= 10^-30, filter.unknown = T, 
              text.cex =0.60, clust = "qval", pie =pie_list, pie.style = "pie")



