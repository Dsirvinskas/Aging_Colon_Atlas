require(ggplot2)
require(dplyr)
require(Matrix)
require(cowplot)
require(Seurat)

load(file="/data/tmp/Dovydas/Single_Cell_RNAseq/2020_02_05_YvO_Single_Cell/Seurat/With_old_4_best_mice/Merged.RData")

# Loading cluster information for different cell populations
Initial = read.csv(file = "5kmeans_initial.csv") # Initial 5 k-means 
StemTa = read.csv(file = "4kmeans_Stem_TA.csv") # 4 k means of Stem/TA cluster
TuftEE = read.csv(file = "2kmeans_Tuft_EE.csv") # 2 k means of Tuft/EE cluster

## Overlaying the initial 5 k-means clustering info on cells
tenmeans = Initial
tenmeans[,2] = as.factor(tenmeans[,2])
row.names(tenmeans) = tenmeans[,1]
tenmeans = tenmeans[row.names(tenmeans) %in% row.names(SI.integrated@meta.data),]

Idents(SI.integrated) = tenmeans[,2]
levels(SI.integrated@active.ident) = c("Colonocyte_1", "Goblet", "Colonocyte_2", "Tuft-EE", "Stem-TA")

# Generating csv files for number of each cell type in different mice (Initial clustering). 
num = cbind(SI.integrated@meta.data, "Cluster" = SI.integrated@active.ident)
info = table(num[,c(5,6,9)])
write.csv(info, "Mouse_Cell_Numbers.csv")




#################################### Re-clustering Tuft/EE ##############################################
Tuft = subset(SI.integrated, idents = "Tuft-EE")
tenmeans = TuftEE
tenmeans[,2] = as.factor(tenmeans[,2])
row.names(tenmeans) = tenmeans[,1]
tenmeans = tenmeans[row.names(tenmeans) %in% row.names(Tuft@meta.data),]

Idents(Tuft) = tenmeans[,2]
levels(Tuft@active.ident) = c("EnteroEndocrine_Cell", "Tuft_Cell")

Tuft = NormalizeData(Tuft, verbose = F)
Tuft = FindVariableFeatures(Tuft, verbose = F)
Tuft = ScaleData(Tuft, features = row.names(SI.Seurat@assays$RNA@data), verbose = F)
Tuft = RunPCA(Tuft, features = row.names(SI.Seurat@assays$RNA@data), verbose = F)
Tuft = RunTSNE(Tuft, dims = 1:30)

pdf(file = "Tuft-EE-tSNE_Reclust.pdf")
DimPlot(Tuft, reduction = "tsne") + theme(legend.position = "bottom")
DimPlot(Tuft, reduction = "tsne", group.by = "Compartment") + theme(legend.position = "bottom")
DimPlot(Tuft, reduction = "tsne", group.by = "Age") + theme(legend.position = "bottom")
FeaturePlot(Tuft, features = c("Dclk1", "Cd24a", "Chga", "Chgb"))
dev.off()


Tuft1 = subset(Tuft, idents = "Tuft_Cell")

Tuft1 = FindNeighbors(Tuft1, dims = 1:10)
Tuft1 =  FindClusters(Tuft1, resolution = 0.1)
pdf(file = paste("Graph_based_Tuft1.pdf", sep = "") )
plot = DimPlot(Tuft1, reduction = "tsne") 
print(plot)
dev.off() 

Tuft1.markers <- FindAllMarkers(Tuft1, only.pos = F, min.pct = 0.1, logfc.threshold = 0.1) 
Tuft1.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) 

top10 <- Tuft1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) 
write.csv(top10, file= "Tuft1_compartment_genes.csv")

Tuft1@meta.data$Compartment = factor(Tuft1@meta.data$Compartment, levels = c("Cecum", "ProxCol", "DistCol"))



col.plot = c("Red", "Black", "Blue") 
pdf(file = "Heatmap_Tuft1.pdf" )
plot = DoHeatmap(Tuft1, features = top10$gene, size=2, angle=0, group.by = "Compartment") + NoLegend() + scale_fill_gradient2( low = "Blue", mid = "Black", high = "Red", midpoint = 0, guide = "colourbar", aesthetics = "fill") +theme(axis.text = element_text(size=5)) 
print(plot)
dev.off() 














EnteroE = subset(Tuft, idents = "EnteroEndocrine_Cell")
EnteroE= NormalizeData(EnteroE, verbose = F)
EnteroE = FindVariableFeatures(EnteroE, verbose = F)
EnteroE = ScaleData(EnteroE, features = row.names(SI.Seurat@assays$RNA@data), verbose = F)
EnteroE = RunPCA(EnteroE, features = row.names(SI.Seurat@assays$RNA@data), verbose = F)
EnteroE = RunTSNE(EnteroE, dims = 1:30)

pdf(file = "EnteroE-tSNE_Reclust.pdf")
DimPlot(EnteroE, reduction = "tsne") + theme(legend.position = "bottom")
DimPlot(EnteroE, reduction = "tsne", group.by = "Compartment") + theme(legend.position = "bottom")
DimPlot(EnteroE, reduction = "tsne", group.by = "Age") + theme(legend.position = "bottom")
dev.off()


pdf(file = "EE-subtype-markers.pdf") ##From Haber paper
FeaturePlot(EnteroE, features = c("Clca4", "Slc12a2", "Sox4", "Dll1"), reduction = "tsne") #Early + Mid progenitor
FeaturePlot(EnteroE, features = c("Neurod1", "Neurod2", "Serpina1c", "Ghrl"), reduction = "tsne") #Late + A progenitor
FeaturePlot(EnteroE, features = c("Cck", "Gal", "Pyy", "Gcg"), reduction = "tsne") #SILA + SIL-P
FeaturePlot(EnteroE, features = c("Bdnf", "Scgn", "Gip", "Fabp5"), reduction = "tsne") #SIK-P + SIK
FeaturePlot(EnteroE, features = c("Sst", "Iapp", "Nts", "Adgrd1"), reduction = "tsne") #SAKD + SIN
FeaturePlot(EnteroE, features = c("Tac1", "Gch1", "Reg4", "Afp"), reduction = "tsne") #EC + EC-Reg4
dev.off()



pdf(file = "Tuft-subtype-markers.pdf")
FeaturePlot(Tuft, features = c("Ptprc", "Dclk1", "Cd24a", "Il25"), reduction = "tsne")
FeaturePlot(Tuft, features = c("Tslp", "Rac2", "Ptgs1", "Irf7"), reduction = "tsne")
FeaturePlot(Tuft, features = c("Nradd", "Gng13", "Nrep", "Rgs2"), reduction = "tsne")
dev.off()

# Generating csv files for number of each cell type in different mice (Tuft). 
num = cbind(Tuft@meta.data, "Cluster" = Tuft@active.ident)
info = table(num[,c(5,6,9)])
write.csv(info, "Mouse_Tuft-EE_Numbers.csv")

# Wilcox
Tuft = subset(SI.integrated, idents = "Tuft-EE")
tenmeans = TuftEE
tenmeans[,2] = as.factor(tenmeans[,2])
row.names(tenmeans) = tenmeans[,1]
tenmeans = tenmeans[row.names(tenmeans) %in% row.names(Tuft@meta.data),]

Idents(Tuft) = tenmeans[,2]
levels(Tuft@active.ident) = c("EnteroEndocrine_Cell", "Tuft_Cell")


Tuft@meta.data$Age[which(Tuft@meta.data$Age == "Young_20w")] = "Young"
DE.Tuft = subset(Tuft, cells = which(Tuft@active.ident == "Tuft_Cell"))

a=c()
a = FindMarkers(object = DE.Tuft,test.use = "wilcox", group.by = 'Age', ident.1 = "Young", ident.2 = "Old")
a = na.omit(a[a$p_val_adj <0.05,] )
a = cbind(a, "Gene"=row.names(a))
a$Gene = as.character(a$Gene)
write.csv(x = a ,file = "Tuft_YvODE_wilcox_genes.csv")

for (i in c("Cecum", "ProxCol", "DistCol")){
  DE.Tuft_Comp = subset(DE.Tuft, cells = which(DE.ColProg@meta.data$Compartment == i))
  a=c()
  a = FindMarkers(object = DE.Tuft_Comp ,test.use = "wilcox", group.by = 'Age', ident.1 = "Young", ident.2 = "Old")
  a = na.omit(a[a$p_val_adj <0.05,] )
  a = cbind(a, "Gene"=row.names(a))
  a$Gene = as.character(a$Gene)
  write.csv(x = a ,file = paste(i,"_Tuft_YvODE_wilcox_genes.csv", sep = ""))
}




DE.EE = subset(Tuft, cells = which(Tuft@active.ident == "EnteroEndocrine_Cell"))
a=c()
a = FindMarkers(object = DE.EE,test.use = "wilcox", group.by = 'Age', ident.1 = "Young", ident.2 = "Old")
a = na.omit(a[a$p_val_adj <0.05,] )
a = cbind(a, "Gene"=row.names(a))
a$Gene = as.character(a$Gene)
write.csv(x = a ,file = "EE_YvODE_wilcox_genes.csv")

for (i in c("Cecum", "ProxCol", "DistCol")){
  DE.EE_Comp = subset(DE.EE, cells = which(DE.ColProg@meta.data$Compartment == i))
  a=c()
  a = FindMarkers(object = DE.Tuft_Comp ,test.use = "wilcox", group.by = 'Age', ident.1 = "Young", ident.2 = "Old")
  a = na.omit(a[a$p_val_adj <0.05,] )
  a = cbind(a, "Gene"=row.names(a))
  a$Gene = as.character(a$Gene)
  write.csv(x = a ,file = paste(i,"_EE_YvODE_wilcox_genes.csv", sep = ""))
}



#################################### Overlay Tuft-EE on initial ##############################################
reclust = cbind(SI.integrated@meta.data, "reclust" = c(rep(3,nrow(SI.integrated@meta.data))))
levels(reclust$reclust) = c(1:3)
tenmeans = TuftEE
tenmeans[,2] = as.factor(tenmeans[,2])
row.names(tenmeans) = tenmeans[,1]
tenmeans = tenmeans[row.names(tenmeans) %in% row.names(SI.integrated@meta.data),]

for ( i in 1:nrow(tenmeans)) {
  a = tenmeans[i,2]
  b = which(row.names(reclust) == row.names(tenmeans)[i])
  reclust[b,9] = a
}
Idents(SI.integrated) = reclust$reclust
levels(SI.integrated) = c(1:3)
SI.integrated@meta.data =  cbind(SI.integrated@meta.data, "cluster" = SI.integrated@active.ident)

levels(SI.integrated@meta.data$cluster) = c("Enteroendocrine", "Tuft", "Other")

pdf(file = "Tuft-EE-tSNE.pdf")
DimPlot(SI.integrated, reduction = "tsne", cols = c("Red", "Blue", "#E9E9E9" ), group.by = "cluster") + theme(legend.position = "bottom")
dev.off()

#################################### Overlay Stem-TA on initial ##############################################
reclust = cbind(SI.integrated@meta.data, "reclust" = c(rep(5,nrow(SI.integrated@meta.data))))


tenmeans = StemTa
tenmeans[,2] = as.factor(tenmeans[,2])
row.names(tenmeans) = tenmeans[,1]
tenmeans = tenmeans[row.names(tenmeans) %in% row.names(SI.integrated@meta.data),]

for ( i in 1:nrow(tenmeans)) {
  a = tenmeans[i,2]
  b = which(row.names(reclust) == row.names(tenmeans)[i])
  reclust[b,10] = a
}

Idents(SI.integrated) = reclust$reclust
levels(SI.integrated) = c(1:5)
### Remove cluster column if already added by Tuft code
##### SI.integrated@meta.data = SI.integrated@meta.data[,-which(colnames(SI.integrated@meta.data) == "cluster")]
SI.integrated@meta.data =  cbind(SI.integrated@meta.data, "cluster" = SI.integrated@active.ident)

levels(SI.integrated@meta.data$cluster) = c("Colonocyte_Precursor", "Goblet_Precursor", "Stem_cell", "TA-cell", "Other")

pdf(file = "Stem-TA-tSNE.pdf")
DimPlot(SI.integrated, reduction = "tsne", cols = c("Red", "Blue", "Green", "Magenta", "#E9E9E9" ), group.by = "cluster") + theme(legend.position = "bottom",legend.title = element_text(size=5))
dev.off()

######################## csv with Stem/TA numbers ################################
### Remove cluster column if already added by Tuft code
##### SI.integrated@meta.data = SI.integrated@meta.data[,-which(colnames(SI.integrated@meta.data) == "cluster")]
Stem = subset(SI.integrated, idents = "Stem-TA")

tenmeans = StemTa
tenmeans[,2] = as.factor(tenmeans[,2])
row.names(tenmeans) = tenmeans[,1]
tenmeans = tenmeans[row.names(tenmeans) %in% row.names(Stem@meta.data),]

Idents(Stem) = tenmeans[,2]
levels(Stem@active.ident) = c("Colonocyte_Precursor", "Goblet_Precursor", "Stem_cell", "TA-cell")

num = cbind(Stem@meta.data, "Cluster" = Stem@active.ident)
num$Cluster = droplevels(num$Cluster)
info = table(num[,c(5,6,9)])
write.csv(info, "Mouse_Stem-TA_Numbers.csv")

##### VLN plots ######

levels(Stem@active.ident) = c(1,2,3,4)

cols = c("Red", "Blue", "Green", "Magenta")
pdf(file= "4_mean_Stem-TA_Violin.pdf")
VlnPlot(Stem, features = c("Lgr5", "Sox4", "Slc12a2","Ascl2","Mki67", "Pcna"), pt.size = 0.001, cols = cols)
VlnPlot(Stem, features = c("Lgr5", "Sox4", "Slc12a2","Ascl2","Mki67", "Pcna"), pt.size = 0,cols = cols)

VlnPlot(Stem, features = c("Dll1", "Atoh1", "Reg4","Muc2","Selenbp1", "Lgals3"), pt.size = 0.001, cols = cols)
VlnPlot(Stem, features = c("Dll1", "Atoh1", "Reg4","Muc2","Selenbp1", "Lgals3"), pt.size = 0,cols = cols)

VlnPlot(Stem, features = c("Slc37a2", "Cyp2c55", "Car1", "Aqp4", "Muc3", "Ggh"), pt.size = 0.001, cols = cols)
VlnPlot(Stem, features = c("Slc37a2", "Cyp2c55", "Car1", "Aqp4", "Muc3", "Ggh"), pt.size = 0,cols = cols)
dev.off()

pdf(file="4_mean_Stem-TA_Ridge.pdf")
RidgePlot(Stem, features = c("Lgr5", "Sox4", "Slc12a2","Ascl2","Mki67", "Pcna"), cols = cols)
RidgePlot(Stem, features = c("Dll1", "Atoh1", "Reg4","Muc2","Selenbp1", "Lgals3"), cols = cols)
RidgePlot(Stem, features = c("Slc37a2", "Cyp2c55", "Car1", "Aqp4", "Muc3", "Ggh"), cols = cols)
dev.off()


## Wilcox
Stem@meta.data$Age[which(Stem@meta.data$Age == "Young_20w")] = "Young"

DE.Stem = subset(Stem, cells = which(Stem@active.ident == "Stem_cell"))
a=c()
a = FindMarkers(object = DE.Stem ,test.use = "wilcox", group.by = 'Age', ident.1 = "Young", ident.2 = "Old")
a = na.omit(a[a$p_val_adj <0.05,] )
a = cbind(a, "Gene"=row.names(a))
a$Gene = as.character(a$Gene)
write.csv(x = a ,file = "Stem_YvODE_wilcox_genes.csv")

for (i in c("Cecum", "ProxCol", "DistCol")){
  DE.StemComp = subset(DE.Stem, cells = which(DE.Stem@meta.data$Compartment == i))
  a=c()
  a = FindMarkers(object = DE.StemComp ,test.use = "wilcox", group.by = 'Age', ident.1 = "Young", ident.2 = "Old")
  a = na.omit(a[a$p_val_adj <0.05,] )
  a = cbind(a, "Gene"=row.names(a))
  a$Gene = as.character(a$Gene)
  write.csv(x = a ,file = paste(i,"_Stem_YvODE_wilcox_genes.csv", sep = ""))
}


DE.TA_cell = subset(Stem, cells = which(Stem@active.ident == "TA-cell"))
a=c()
a = FindMarkers(object = DE.TA_cell ,test.use = "wilcox", group.by = 'Age', ident.1 = "Young", ident.2 = "Old")
a = na.omit(a[a$p_val_adj <0.05,] )
a = cbind(a, "Gene"=row.names(a))
a$Gene = as.character(a$Gene)
write.csv(x = a ,file = "TA_cell_YvODE_wilcox_genes.csv")

for (i in c("Cecum", "ProxCol", "DistCol")){
  DE.TA_Comp = subset(DE.TA_cell, cells = which(DE.TA_cell@meta.data$Compartment == i))
  a=c()
  a = FindMarkers(object = DE.TA_Comp ,test.use = "wilcox", group.by = 'Age', ident.1 = "Young", ident.2 = "Old")
  a = na.omit(a[a$p_val_adj <0.05,] )
  a = cbind(a, "Gene"=row.names(a))
  a$Gene = as.character(a$Gene)
  write.csv(x = a ,file = paste(i,"_TA_YvODE_wilcox_genes.csv", sep = ""))
}



DE.ColProg = subset(Stem, cells = which(Stem@active.ident == "Colonocyte_Precursor"))
a=c()
a = FindMarkers(object = DE.ColProg ,test.use = "wilcox", group.by = 'Age', ident.1 = "Young", ident.2 = "Old")
a = na.omit(a[a$p_val_adj <0.05,] )
a = cbind(a, "Gene"=row.names(a))
a$Gene = as.character(a$Gene)
write.csv(x = a ,file = "ColProg_YvODE_wilcox_genes.csv")

for (i in c("Cecum", "ProxCol", "DistCol")){
  DE.ColProg_Comp = subset(DE.ColProg, cells = which(DE.ColProg@meta.data$Compartment == i))
  a=c()
  a = FindMarkers(object = DE.ColProg_Comp ,test.use = "wilcox", group.by = 'Age', ident.1 = "Young", ident.2 = "Old")
  a = na.omit(a[a$p_val_adj <0.05,] )
  a = cbind(a, "Gene"=row.names(a))
  a$Gene = as.character(a$Gene)
  write.csv(x = a ,file = paste(i,"_ColProg_YvODE_wilcox_genes.csv", sep = ""))
}



DE.GobProg = subset(Stem, cells = which(Stem@active.ident == "Goblet_Precursor"))
a=c()
a = FindMarkers(object = DE.GobProg ,test.use = "wilcox", group.by = 'Age', ident.1 = "Young", ident.2 = "Old")
a = na.omit(a[a$p_val_adj <0.05,] )
a = cbind(a, "Gene"=row.names(a))
a$Gene = as.character(a$Gene)
write.csv(x = a ,file = "GobProg_YvODE_wilcox_genes.csv")

for (i in c("Cecum", "ProxCol", "DistCol")){
  DE.GobProg_Comp = subset(DE.GobProg, cells = which(DE.ColProg@meta.data$Compartment == i))
  a=c()
  a = FindMarkers(object = DE.GobProg_Comp ,test.use = "wilcox", group.by = 'Age', ident.1 = "Young", ident.2 = "Old")
  a = na.omit(a[a$p_val_adj <0.05,] )
  a = cbind(a, "Gene"=row.names(a))
  a$Gene = as.character(a$Gene)
  write.csv(x = a ,file = paste(i,"_GobProg_YvODE_wilcox_genes.csv", sep = ""))
}


############ Separate tSNE by sample / age / batch ##################
tenmeans = Initial
tenmeans[,2] = as.factor(tenmeans[,2])
row.names(tenmeans) = tenmeans[,1]
tenmeans = tenmeans[row.names(tenmeans) %in% row.names(SI.integrated@meta.data),]

Idents(SI.integrated) = tenmeans[,2]
levels(SI.integrated@active.ident) = c("Colonocyte_1", "Goblet", "Colonocyte_2", "Tuft-EE", "Stem-TA")
pdf(file = "Separate_Sample_tSNE.pdf")
DimPlot(SI.integrated, reduction = "tsne", split.by = "Sample", ncol = 3, pt.size = 0.01) + theme(legend.position = "bottom")
dev.off()


SI.integrated@meta.data$Age[which(SI.integrated@meta.data$Age == "Young_20w")] = "Young"
pdf(file = "Separate_Age_tSNE.pdf")
DimPlot(SI.integrated, reduction = "tsne", split.by = "Age", ncol = 3, pt.size = 0.01) + theme(legend.position = "bottom")
dev.off()

SI.integrated@meta.data$Batch = 1
SI.integrated@meta.data$Batch[SI.integrated@meta.data$Sample %in% c("Old_1","Old_2", "Old_3", "Young_1", "Young_2", "Young_3")] = 2
pdf(file = "Separate_Experiment_tSNE.pdf")
DimPlot(SI.integrated, reduction = "tsne", split.by = "Batch", ncol = 2, pt.size = 0.01) + theme(legend.position = "bottom", lab) + labs(title = "Batch")
dev.off()



######################## Colonocyte by compartment ################################

Colonocyte = subset(SI.integrated, idents = c("Colonocyte_1","Colonocyte_2")) 

Idents(Colonocyte)= Colonocyte@meta.data$Compartment
Colonocyte@meta.data$Compartment = factor(Colonocyte@meta.data$Compartment, levels = c("Cecum", "ProxCol", "DistCol"))


pdf(file = "Colonocyte_TSNE.pdf")
DimPlot(Colonocyte, reduction = "tsne", pt.size = 0.1, group.by = "Compartment")
dev.off()

Colonocyte.markers <- FindAllMarkers(Colonocyte, only.pos = F, min.pct = 0.1, logfc.threshold = 0.1) 
Colonocyte.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) 

top10 <- Colonocyte.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) 
write.csv(top10, file= "Colonocyte_compartment_genes.csv")

Colonocyte@meta.data$Compartment = factor(Colonocyte@meta.data$Compartment, levels = c("Cecum", "ProxCol", "DistCol"))



col.plot = c("Red", "Black", "Blue") 
pdf(file = "Heatmap_Colonocyte1.pdf" )
plot = DoHeatmap(Colonocyte, features = top10$gene, size=2, angle=0, group.by = "Compartment") + NoLegend() + scale_fill_gradient2( low = "Blue", mid = "Black", high = "Red", midpoint = 0, guide = "colourbar", aesthetics = "fill") +theme(axis.text = element_text(size=5)) 
print(plot)
dev.off() 


Colonocyte= NormalizeData(Colonocyte, verbose = F)
Colonocyte = FindVariableFeatures(Colonocyte, verbose = F)
Colonocyte = ScaleData(Colonocyte, features = row.names(SI.Seurat@assays$RNA@data), verbose = F)
Colonocyte = RunPCA(Colonocyte, features = row.names(SI.Seurat@assays$RNA@data), verbose = F)
Colonocyte = RunTSNE(Colonocyte, dims = 1:30)

pdf(file = "Colonocyte_TSNE1.pdf")
DimPlot(Colonocyte, reduction = "tsne", pt.size = 0.1, group.by = "Compartment")
dev.off()

Colonocyte = FindNeighbors(Colonocyte, dims = 1:10)
Colonocyte =  FindClusters(Colonocyte, resolution = 0.035)
pdf(file = paste("Graph_based_Colonocyte.pdf", sep = "") )
plot = DimPlot(Colonocyte, reduction = "tsne") 
print(plot)
dev.off() 

Colonocyte.markers <- FindAllMarkers(Colonocyte, only.pos = F, min.pct = 0.1, logfc.threshold = 0.1) 
Colonocyte.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) 

top10 <- Colonocyte.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) 
col.plot = c("Red", "Black", "Blue") 
pdf(file = "Heatmap_Colonocyte_Graph.pdf" )
plot = DoHeatmap(Colonocyte, features = top10$gene, size=2, angle=0) + NoLegend() + scale_fill_gradient2( low = "Blue", mid = "Black", high = "Red", midpoint = 0, guide = "colourbar", aesthetics = "fill") +theme(axis.text = element_text(size=5)) 
print(plot)
dev.off() 

Colonocyte_clean = subset(Colonocyte, idents=c(0,1))
pdf(file = "Colonocyte_clean_TSNE.pdf")
DimPlot(Colonocyte_clean, reduction = "tsne", pt.size = 0.1, group.by = "Compartment")
dev.off()

Colonocyte_clean = FindNeighbors(Colonocyte_clean, dims = 1:10)
Colonocyte_clean =  FindClusters(Colonocyte_clean, resolution = 0.03)
pdf(file = paste("Graph_based_Colonocyte_clean.pdf", sep = "") )
plot = DimPlot(Colonocyte_clean, reduction = "tsne") 
print(plot)
dev.off() 

identities = cbind(Colonocyte_clean@meta.data, "Test" = Colonocyte_clean@active.ident)
levels(identities$Test) = c(0,1,"A","B")

Colonocyte_zero = subset(Colonocyte_clean, idents=0)
Colonocyte_zero = FindNeighbors(Colonocyte_zero, dims = 1:10)
Colonocyte_zero =  FindClusters(Colonocyte_zero, resolution = 0.03)
levels(Colonocyte_zero@active.ident) = c("A", "B")
pdf(file = paste("Graph_based_Colonocyte_zero.pdf", sep = "") )
plot = DimPlot(Colonocyte_zero, reduction = "tsne") 
print(plot)
dev.off() 

a = cbind(Colonocyte_zero@meta.data, "test" = Colonocyte_zero@active.ident)


for (i in 1:nrow(a)) {
  b = which(row.names(a)[i] == row.names(identities))
  identities$Test[b] = a$test[i]
}
identities$Test = droplevels(identities$Test)

Idents(Colonocyte_clean) = identities$Test
pdf(file = paste("Graph_based_Colonocyte_Merge.pdf", sep = "") )
plot = DimPlot(Colonocyte_clean, reduction = "tsne") 
print(plot)
dev.off() 

Colonocyte_clean.markers <- FindAllMarkers(Colonocyte_clean, only.pos = F, min.pct = 0.1, logfc.threshold = 0.1) 
Colonocyte_clean.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) 

top10 <- Colonocyte_clean.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) 
col.plot = c("Red", "Black", "Blue") 

Colonocyte_clean@meta.data = cbind(Colonocyte_clean@meta.data, "Clust" = Colonocyte_clean@active.ident)
levels(Colonocyte_clean@meta.data$Clust) = c("Cecum-Specific","DistCol-Specific", "ProxCol-Specific")
Colonocyte_clean@meta.data$Clust = factor(Colonocyte_clean@meta.data$Clust, levels = c("Cecum-Specific", "ProxCol-Specific","DistCol-Specific"))


genes = top10$gene[c(1:10,21:30,11:20)]
pdf(file = "Heatmap_Colonocyte_Graph_merge.pdf" )
plot = DoHeatmap(Colonocyte_clean, features = genes, size=2, angle=0, group.by = "Clust") + NoLegend() + scale_fill_gradient2( low = "Blue", mid = "Black", high = "Red", midpoint = 0, guide = "colourbar", aesthetics = "fill") +theme(axis.text = element_text(size=5)) 
print(plot)
dev.off() 


pdf("3colonocytes.pdf")
DimPlot(Colonocyte_clean, reduction = "tsne", group.by = "Clust")
DimPlot(Colonocyte_clean, reduction = "tsne", group.by = "Compartment")
dev.off()

write.csv(table(Colonocyte_clean@meta.data[,c(5,6,12)]), file = "3colonocytes_numbers.csv")

a = Colonocyte_clean@meta.data[,c(1,12)]
write.csv(a, file="3colonocytes_each_cell.csv")

######################## Colonocyte Y v O DEG's ##########################
#For Wilcox
Colonocyte_clean@meta.data$Age[which(Colonocyte_clean@meta.data$Age == "Young_20w")] = "Young"
DE.cecum = subset(Colonocyte_clean, cells = which(Colonocyte_clean@meta.data$Compartment == "Cecum"))

a=c()
a = FindMarkers(object = DE.cecum,test.use = "wilcox", group.by = 'Age', ident.1 = "Young", ident.2 = "Old")
a = na.omit(a[a$p_val_adj <0.05,] )
a = cbind(a, "Gene"=row.names(a))
a$Gene = as.character(a$Gene)
write.csv(x = a ,file = "Cecum_enrich_colonocyte_YvODE_wilcox_genes.csv")


DE.ProxCol = subset(Colonocyte_clean, cells = which(Colonocyte_clean@meta.data$Compartment == "ProxCol"))
a=c()
a = FindMarkers(object = DE.ProxCol,test.use = "wilcox", group.by = 'Age', ident.1 = "Young", ident.2 = "Old")
a = na.omit(a[a$p_val_adj <0.05,] )
a = cbind(a, "Gene"=row.names(a))
a$Gene = as.character(a$Gene)
write.csv(x = a ,file = "ProxCol_enrich_colonocyte_YvODE_wilcox_genes.csv")


DE.DistCol = subset(Colonocyte_clean, cells = which(Colonocyte_clean@meta.data$Compartment == "DistCol"))
a=c()
a = FindMarkers(object = DE.DistCol,test.use = "wilcox", group.by = 'Age', ident.1 = "Young", ident.2 = "Old")
a = na.omit(a[a$p_val_adj <0.05,] )
a = cbind(a, "Gene"=row.names(a))
a$Gene = as.character(a$Gene)
write.csv(x = a ,file = "DistCol_enrich_colonocyte_YvODE_wilcox_genes.csv")

######################## Colonocyte by compartment enrichment - on initial clust################################

reclust = cbind(SI.integrated@meta.data, "reclust" = c(rep(4,nrow(SI.integrated@meta.data))))
levels(reclust$reclust) = c(1:4)

tenmeans = data.frame("Cluster" = Colonocyte_clean@meta.data$Clust)
row.names(tenmeans) = row.names(Colonocyte_clean@meta.data)

for ( i in 1:nrow(tenmeans)) {
  a = tenmeans[i,1]
  b = which(row.names(reclust) == row.names(tenmeans)[i])
  reclust[b,9] = a
}
Idents(SI.integrated) = reclust$reclust
levels(SI.integrated) = c(1:4)
SI.integrated@meta.data =  cbind(SI.integrated@meta.data, "cluster" = SI.integrated@active.ident)


levels(SI.integrated@meta.data$cluster) = c("Cecum-enriched", "ProxCol-enriched", "DistCol-enriched", "Other")

pdf(file = "3colonocyte-on-initial.pdf")
DimPlot(SI.integrated, reduction = "tsne", cols = c("Red","Green", "Blue", "#E9E9E9" ), group.by = "cluster")+ theme(legend.position = "bottom") 

dev.off()





########### YvO all cells



a=c()
a = FindMarkers(object = SI.integrated,test.use = "wilcox", group.by = 'Age', ident.1 = "Young", ident.2 = "Old")
a = na.omit(a[a$p_val_adj <0.05,] )
a = cbind(a, "Gene"=row.names(a))
a$Gene = as.character(a$Gene)
write.csv(x = a ,file = "AllCells_YvODE_wilcox_genes.csv")

############# Goblet ridge plots

pdf(file= "Goblet_Ridge.pdf")
RidgePlot(SI.integrated, features = c("Muc2", "Agr2"), group.by="Age", idents = "Goblet")
RidgePlot(SI.integrated, features = c("Muc2"), group.by="Age")
RidgePlot(SI.integrated, features = c("Agr2"), group.by="Age")
dev.off()

SI.comp = subset(SI.integrated, cells = which(SI.integrated@meta.data$Compartment == "DistCol"))
SI.comp@meta.data$Age[which(SI.comp@meta.data$Age == "Young_20w")] = "Young" 
pdf(file= "Goblet_Ridge_DistCol.pdf")
VlnPlot(SI.comp, features = c("Muc2", "Agr2"), group.by="Age", idents = "Goblet", pt.size = 1)
dev.off()


############## Colonocytes all together - common markers ########

SI.Colonocyte = SI.integrated
SI.Colonocyte@meta.data = cbind(SI.Colonocyte@meta.data, "Clust" = SI.Colonocyte@active.ident)
SI.Colonocyte@meta.data$Clust[which(SI.Colonocyte@meta.data$Clust == "Colonocyte_2")] = "Colonocyte_1"
SI.Colonocyte@meta.data$Clust = droplevels(SI.Colonocyte@meta.data$Clust)
levels(SI.Colonocyte@meta.data$Clust) = c("Colonocyte", "Goblet", "Tuft-EE", "Stem-TA")
Idents(SI.Colonocyte) = SI.Colonocyte@meta.data$Clust

SI.Colonocyte.markers <- FindAllMarkers(SI.Colonocyte, only.pos = F, min.pct = 0.1, logfc.threshold = 0.1) 
SI.Colonocyte.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) 

top10 <- SI.Colonocyte.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) 
write.csv(top10, file= "Merged_Colonocyte_genes.csv")

col.plot = c("Red", "Black", "Blue") 
pdf(file = "Heatmap_Merged_Colonocyte.pdf" )
plot = DoHeatmap(SI.Colonocyte, features = top10$gene, size=2, angle=0) + NoLegend() + scale_fill_gradient2( low = "Blue", mid = "Black", high = "Red", midpoint = 0, guide = "colourbar", aesthetics = "fill") +theme(axis.text = element_text(size=5)) 
print(plot)
dev.off() 


pdf(file= "Colonocyte_Vln.pdf")
VlnPlot(SI.Colonocyte, features = c("Alpi", "Slc2a1", "Car2", "Krt20", "Slc26a3"), group.by="Clust", pt.size = 0.1)
VlnPlot(SI.Colonocyte, features = c("Alpi", "Slc2a1", "Car2", "Krt20", "Slc26a3"), group.by="Clust", pt.size = 0)
dev.off()

pdf(file="Colonocyte_YvO_Vln.pdf")
VlnPlot(SI.Colonocyte, features = c("Alpi", "Slc2a1", "Car2", "Krt20", "Slc26a3"), group.by="Age", idents = "Colonocyte", pt.size = 0)
dev.off()