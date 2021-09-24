require(ggplot2)
require(dplyr)
require(Matrix)
require(cowplot)
require(Seurat)


## Loading new data
load("/data/tmp/Dovydas/Single_Cell_RNAseq/2020_02_05_YvO_Single_Cell/Seurat/StartingData.RData")
meta.SI = hash
meta.SI = meta.SI[c(which(meta.SI$Compartment == "Cecum"),which(meta.SI$Compartment == "ProxCol"),which(meta.SI$Compartment == "DistCol")),]
row.names(meta.SI) = meta.SI[,1]
meta.SI = meta.SI[,-1]
meta.SI = cbind(meta.SI, "Age" = meta.SI[,4])
meta.SI$Age = as.character(meta.SI$Age)

meta.SI$Age[which(substr(x = meta.SI$Age, start = 1, stop = 3)=="You")] = "Young"
meta.SI$Age[which(substr(x = meta.SI$Age, start = 1, stop = 3)=="Old")] = "Old"

data1 = data[,colnames(data) %in% row.names(meta.SI)]
data1 = data1[,match(colnames(data1),row.names(meta.SI))]


nams = colnames(data1)
for (i in 1:length(nams)){
  nams[i] = paste(nams[i],"2", sep = "-")
}
colnames(data1) = nams
nams = c()

nams = row.names(meta.SI)
for (i in 1:length(nams)){
  nams[i] = paste(nams[i],"2", sep = "-")
}
row.names(meta.SI) = nams
nams = c()

meta.SI = meta.SI[,-1]
meta.SI = droplevels.data.frame(meta.SI) #levels for compartment = Cecum DistCol and ProxCol

## Loading old data
info = read.csv(file= "/data/tmp/Dovydas/Single_Cell_RNAseq/2019_05_03_Single-Cell_RNAseq/Cecum and Colon/p1.cryptsAll.hashtag.csv", header=TRUE)
info = cbind(info,c(rep(0,nrow(info))))
colnames(info)[16] = "Sample"
for (i in 1:nrow(info)){
  info[i,16] = substr(info[i,1], 18,18)
}

info_1 = info[,-(2:13)]
info_1 = info_1[which(info_1$type == "single"),]
info_1 = droplevels(info_1)

info_1[which(info_1$TYPE=="shash7"),3] = "shash6"
info_1 = droplevels(info_1)

levels = c('Prox_SI','Int_SI','Dist_SI','Cecum', 'Prox_Col', 'Dist_Col')
levels(info_1$TYPE) = levels

a = c(rep(0,nrow(info_1)))
info_1 = cbind(info_1, a)
colnames(info_1)[5] = "Age"

a= c(which(info_1$Sample == 1), which(info_1$Sample == 2), which(info_1$Sample==6))
info_1[a,5] = "Old"

a= c(which(info_1$Sample == 3 ), which(info_1$Sample ==4), which(info_1$Sample ==5))
info_1[a,5] = "Young"

a=c(which(info_1$Sample == 7)) 
info_1[a,5] = "Young_20w" 

b = c(which(info_1$Sample == 1)) 
info_1[b,4] = "946_90w" 

b= c(which(info_1$Sample ==2)) 
info_1[b,4] = "948_92w" 

b= c(which(info_1$Sample ==3)) 
info_1[b,4] = "952_10w" 
b= c(which(info_1$Sample ==4)) 
info_1[b,4] = "953_8w" 
b= c(which(info_1$Sample ==5)) 
info_1[b,4] = "954_10w" 
b= c(which(info_1$Sample ==6)) 
info_1[b,4] = "951_90w" 
b= c(which(info_1$Sample ==7)) 
info_1[b,4] = "965_20w" 

info_1 = info_1[,-2] 

colnames(info_1)[2] = "Compartment" 
info_1$Sample = as.factor(info_1$Sample) 


meta_data = info_1[,2:4] 
row.names(meta_data) = info_1$sampleID 

meta_data_SI = meta_data[c(which(meta_data$Compartment == "Cecum"),which(meta_data$Compartment == "Prox_Col"),which(meta_data$Compartment == "Dist_Col")),] 


data = readMM(file= "/data/tmp/Dovydas/Single_Cell_RNAseq/2019_05_03_Single-Cell_RNAseq/Cecum and Colon/matrix.mtx.gz") 
cellID = read.csv(file="/data/tmp/Dovydas/Single_Cell_RNAseq/2019_05_03_Single-Cell_RNAseq/Cecum and Colon/barcodes.tsv.gz", header = FALSE, sep="") 
genesID = read.csv(file="/data/tmp/Dovydas/Single_Cell_RNAseq/2019_05_03_Single-Cell_RNAseq/Cecum and Colon/features.tsv.gz", header = FALSE, sep="") 
colnames(data) = cellID$V1 
rownames(data) = genesID$V2 

print("Starting large intestine analysis") 

SI = data[,colnames(data) %in% row.names(meta_data_SI)] 
SI = SI[,match(colnames(SI),row.names(meta_data_SI))] 



nams1 = colnames(SI)
for (i in 1:length(nams1)){
  nams1[i] = paste(nams1[i],"1", sep = "-")
}
colnames(SI) = nams1
nams1 = c()

nams1 = row.names(meta_data_SI)
for (i in 1:length(nams1)){
  nams1[i] = paste(nams1[i], "1", sep = "-")
}
rownames(meta_data_SI) = nams1
nams1 = c()


Batch = c(rep(NA,nrow(meta_data_SI)))
meta_data_SI = cbind(Batch, meta_data_SI)
meta_data_SI = droplevels.data.frame(meta_data_SI)
levels(meta_data_SI$Compartment) = c("Cecum", "ProxCol", "DistCol")

## Combining the runs

meta.SI$Compartment = as.character(meta.SI$Compartment)
meta.SI$Compartment = factor(meta.SI$Compartment, levels = c("Cecum", "ProxCol", "DistCol"))

meta_data_SI$Compartment = as.character(meta_data_SI$Compartment)
meta_data_SI$Compartment = factor(meta_data_SI$Compartment, levels = c("Cecum", "ProxCol", "DistCol"))

data_comb = cbind(data1, SI)
meta_comb = rbind(meta.SI, meta_data_SI) ##need to check if compartments stayed correct after merge - corrected 2020_07_06

data_comb = data_comb[,colnames(data_comb) %in% row.names(meta_comb)] 

## Running Seurat code    
print("Starting large intestine analysis") 
    
    SI.Seurat = CreateSeuratObject(counts = data_comb, project = "Whole_LI", min.cells = 1, min.features = 200, meta.data = meta_comb) 
    print("Seurat object created") 
    
    SI.Seurat[["percent.mt"]] = PercentageFeatureSet(SI.Seurat, pattern = "^mt-")
    
    Idents(SI.Seurat) = SI.Seurat@meta.data$Sample
    SI.Seurat = subset(SI.Seurat, idents = c("952_10w", "954_10w"), invert = T)
    SI.Seurat = subset(SI.Seurat, subset = nFeature_RNA > 1000 & percent.mt < 15)   
    
    SI.list = SplitObject(SI.Seurat, split.by = "Sample")
    SI.list = lapply(X = SI.list, FUN = function(x){
      x = NormalizeData(x, verbose = FALSE)
      x = FindVariableFeatures(x, verbose = FALSE)
    })
    
    features = SelectIntegrationFeatures(object.list = SI.list)
    SI.list = lapply(X = SI.list, FUN = function(x){
      x = ScaleData(x, features = features, verbose = FALSE)
      x = RunPCA(x, features = features, verbose = FALSE)
    })
    
    #took Young_1 from the new data as reference
    anchors = FindIntegrationAnchors(object.list = SI.list, reference = 2, reduction = "rpca", dims = 1:30)
    SI.integrated = IntegrateData(anchorset = anchors, dims = 1:30)
    
    SI.integrated = ScaleData(SI.integrated, verbose = FALSE, do.scale = FALSE, do.center = FALSE)
    SI.integrated = RunPCA(SI.integrated, verbose = FALSE)
    SI.integrated = RunTSNE(SI.integrated, dims = 1:30)
    
    
    pdf(file = "rPCA.pdf")
    DimPlot(SI.integrated, reduction = "pca", dims = c(1,2), pt.size = 0.01, group.by = "Sample")
    DimPlot(SI.integrated, reduction = "tsne", dims = c(1,2), pt.size = 0.01, group.by = "Sample")
    
    DimPlot(SI.integrated, reduction = "tsne", dims = c(1,2), pt.size = 0.01, group.by = "Age")
    DimPlot(SI.integrated, reduction = "tsne", dims = c(1,2), pt.size = 0.01, group.by = "Compartment")
    dev.off()
    
    pdf(file = "Violin.pdf")
    VlnPlot(SI.Seurat, features = c("percent.mt","nFeature_RNA", "nCount_RNA"))
    dev.off()
    
 
    


pdf(file = "Features_1.pdf")
a = FeaturePlot(SI.integrated, features = c("Lgr5", "Ascl2", "Mki67","Pcna"), pt.size = 0.01, reduction = "tsne")
print(a)
a = FeaturePlot(SI.integrated, features = c("Lyz1", "Defa17", "Defa24", "Lyz2"), pt.size = 0.01, reduction = "tsne")
print(a)
a=FeaturePlot(SI.integrated, features = c("Muc2", "Atoh1", "Spdef", "Agr2"), pt.size = 0.01, reduction = "tsne")
print(a)
a=FeaturePlot(SI.integrated, features = c("Gsdmc4", "Selenbp1", "Lgals3", "Mgat4c"), pt.size = 0.01, reduction = "tsne")
print(a)
a=FeaturePlot(SI.integrated, features = c("Dclk1", "Cd24a", "Chga", "Chgb"), pt.size = 0.01, reduction = "tsne")
print(a)
a=FeaturePlot(SI.integrated, features = c("Reg4", "Hoxb13", "Hoxb6"), pt.size = 0.01, reduction = "tsne")
print(a)
a=FeaturePlot(SI.integrated, features = c("percent.mt", "nFeature_RNA", "nCount_RNA"), pt.size = 0.01, reduction = "tsne")
print(a)
dev.off()



Idents(SI.integrated) = SI.integrated@meta.data$Compartment
SI.integrated.markers <- FindAllMarkers(SI.integrated, only.pos = F, min.pct = 0.1, logfc.threshold = 0.1) 
SI.integrated.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) 

top10 <- SI.integrated.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) 
write.csv(top10, file= "CompartmentTop10_genes_rPCA.csv" )

col.plot = c("Red", "Black", "Blue") 
pdf(file = "Heatmap_Compartment_rPCA.pdf")
plot = DoHeatmap(SI.integrated, features = top10$gene, size=2, angle=0) + NoLegend() + scale_fill_gradient2( low = "Blue", mid = "Black", high = "Red", midpoint = 0, guide = "colourbar", aesthetics = "fill") +theme(axis.text = element_text(size=5)) 
print(plot)
dev.off() 




## K-means
    dat = SI.integrated@assays$integrated@data
    texpdf<-t(dat)
    tot<-data.frame(center=NULL,tot=NULL)
    for (c in 5:8){
      k<-kmeans(texpdf,center=c,nstart=5)
      write.csv(k$cluster,paste(c,"TESTallcels.csv",sep=""))
      print(c)
      
      tenmeans = read.csv(file = paste(c,"TESTallcels.csv", sep = ""))
      tenmeans[,2] = as.factor(tenmeans[,2])
      row.names(tenmeans) = tenmeans[,1]
      tenmeans = tenmeans[row.names(tenmeans) %in% row.names(SI.integrated@meta.data),]
      
      Idents(SI.integrated) = tenmeans[,2]
      cols = c("Red", "Blue", "grey30", "Green", "Yellow", "Magenta", "Cyan", "indianred", "Pink", "Purple", "orangered1", "wheat4", "darkolivegreen2", "slategray3", "lightsalmon2")
      pdf(file = paste("rPCA_PCA_", c,"means.pdf", sep = "") )
      plot = DimPlot(SI.integrated, reduction = "pca", cols = cols) 
      print(plot)
      dev.off() 
      
      pdf(file = paste("rPCA_TSNE_", c,"means.pdf", sep = "")) 
      plot = DimPlot(SI.integrated, reduction = "tsne", cols = cols) 
      print(plot)
      dev.off() 
      
      SI.integrated.markers <- FindAllMarkers(SI.integrated, only.pos = F, min.pct = 0.1, logfc.threshold = 0.1) 
      SI.integrated.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) 
      
      top10 <- SI.integrated.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) 
      write.csv(top10, file= paste(c,"meanTop10_genes_rPCA.csv", sep ="") )
      
      col.plot = c("Red", "Black", "Blue") 
      pdf(file = paste("Heatmap_",c,"means_rPCA.pdf", sep = "") )
      plot = DoHeatmap(SI.integrated, features = top10$gene, size=2, angle=0) + NoLegend() + scale_fill_gradient2( low = "Blue", mid = "Black", high = "Red", midpoint = 0, guide = "colourbar", aesthetics = "fill") +theme(axis.text = element_text(size=5)) 
      print(plot)
      dev.off() 
    }
    

save.image(file="Merged.RData")