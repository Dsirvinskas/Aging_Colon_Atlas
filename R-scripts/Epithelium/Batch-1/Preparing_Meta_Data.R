require(ggplot2)
require(cowplot)
require(dplyr)
require(Seurat, lib.loc = "/data/tmp/Dovydas/Rpackage/")
require(Matrix)

info = read.csv(file= "p1.cryptsAll.hashtag.csv", header=TRUE)

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

meta_data_SI = meta_data[c(which(meta_data$Compartment == "Prox_SI"), which(meta_data$Compartment == "Int_SI"), which(meta_data$Compartment == "Dist_SI")),]
meta_data_LargeInt = meta_data[c(which(meta_data$Compartment == "Cecum"), which(meta_data$Compartment == "Prox_Col"), which(meta_data$Compartment == "Dist_Col")),]

data = readMM(file= "matrix.mtx.gz")
cellID = read.csv(file="barcodes.tsv.gz", header = FALSE, sep="")
genesID = read.csv(file="features.tsv.gz", header = FALSE, sep="")
colnames(data) = cellID$V1
rownames(data) = genesID$V2

data = as.data.frame(as.matrix(data))

#For SI
print("Starting small intestine analysis")

SI = data[,colnames(data) %in% row.names(meta_data_SI)]
SI = SI[,match(colnames(SI),row.names(meta_data_SI))]

print("Pre-processing complete")
save.image(file="Calculations.RData")

SI.Seurat = CreateSeuratObject(counts = SI, project = "Whole_SI", min.cells = 1, min.features = 200, meta.data = meta_data_SI)
SI.Seurat[["percent.mt"]] <- PercentageFeatureSet(SI.Seurat, pattern = "^mt.")

print("Seurat object created")
save.image(file="Calculations.RData")

##Removal of cells with a higher than a specific mitochondrial gene percentage
pdf("Violin_Plot.pdf")
VlnPlot(object = SI.Seurat , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# chose to try working with 15% mito genes
SI.Seurat <- subset(SI.Seurat, subset = nFeature_RNA > 200 & percent.mt < 15)
SI.Seurat = NormalizeData(SI.Seurat)

##Variable feature finding
SI.Seurat = FindVariableFeatures(SI.Seurat, selection.method = "vst", nfeatures = 2000)
top10 = head(VariableFeatures(SI.Seurat), 10)

write.csv(top10, file="Top10VariableGenes.csv")

plot1 = VariableFeaturePlot(SI.Seurat)
plot2 = LabelPoints(plot=plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
pdf("Variable_Plot.pdf")
CombinePlots(plots = list(plot1, plot2))
dev.off()

##Data scaling
all.genes = rownames(SI.Seurat)
SI.Seurat = ScaleData(SI.Seurat, features = all.genes)

print("Data scaled")
save.image(file="Calculations.RData")

SI.Seurat <- RunPCA(SI.Seurat, features = all.genes)

pdf("PCA_Plot")
DimPlot(SI.Seurat, reduction = "pca", cols = c("Red"))
dev.off()

print("PCA generated")
save.image(file="Calculations.RData")


texpdf<-t(SI.Seurat@assays$RNA@data)
tot<-data.frame(center=NULL,tot=NULL)
for (c in 7:10){
  k<-kmeans(texpdf,center=c,nstart=5)
  write.csv(k$cluster,paste(c,"allcels.csv",sep=""))
  print(c)
  tot<-rbind(tot,data.frame(center=c,tot=k$tot.withinss))
  write.csv(tot,"tot.csv")
}
pdf("K-means_TOT_plot")
plot(tot)
dev.off()

print("K-means clustering complete")
save.image(file="Calculations.RData")

tenmeans = read.csv("10allcels.csv")
tenmeans[,2] = as.factor(tenmeans[,2])

Idents(SI.Seurat) = tenmeans[,2]
levels = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10")
levels(SI.Seurat@active.ident) = levels


pdf(file = "PCA_after_kmeans_12means.pdf")
DimPlot(SI.Seurat, reduction = "pca", cols = c("Red", "Blue", "Grey", "Green", "Yellow", "Magenta", "Cyan", "indianred", "Pink", "Purple"))
dev.off()

SI.Seurat.markers <- FindAllMarkers(SI.Seurat, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
SI.Seurat.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

top10 <- SI.Seurat.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
col.plot = c("Red", "Black", "Blue")

pdf(file = "Heatmap_10means.pdf")
DoHeatmap(SI.Seurat, features = top10$gene, size=2, angle=0) + NoLegend() + scale_fill_gradient2( low = "Blue", mid = "Black", high = "Red", midpoint = 0, guide = "colourbar", aesthetics = "fill") +theme(axis.text = element_text(size=5))
dev.off()

print("Kmeans clustering overlaid on PCA")
save.image(file="Calculations.RData")

pdf(file="PCA_Age.pdf")
Idents(SI.Seurat)= SI.Seurat@meta.data$Age
DimPlot(SI.Seurat, reduction = "pca", cols=c("red","blue","grey", "green", "yellow", "magenta"))
dev.off()
pdf(file="PCA_Compartment.pdf")
Idents(SI.Seurat)= SI.Seurat@meta.data$Compartment
DimPlot(SI.Seurat, reduction = "pca", cols=c("red","blue","grey", "green", "yellow", "magenta"))
dev.off()

Idents(SI.Seurat) = tenmeans[,2]
pdf("Violin_Plot_With_Clusters.pdf")
VlnPlot(object = SI.Seurat , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pdf(file="Gene_exp_PCA.pdf")
genes.viz = c("Alpi", "Apoa1", "Apoa4", "Fabp1") ##enterocytes
FeaturePlot(SI.Seurat, genes.viz, pt.size = 0.01, cols = c("Grey","Red"))
genes.viz = c("Muc2", "Clca1", "Tff3", "Agr2") ##goblet
FeaturePlot(SI.Seurat, genes.viz, pt.size = 0.01, cols = c("Grey","Red"))
genes.viz = c("Lyz1", "Defa17", "Defa22", "Defa24", "Ang4") ##Paneth
FeaturePlot(SI.Seurat, genes.viz, pt.size = 0.01, cols = c("Grey","Red"))
genes.viz = c("Chga", "Chgb", "Tac1", "Tph1", "Neurog3") ##EE
FeaturePlot(SI.Seurat, genes.viz, pt.size = 0.01, cols = c("Grey","Red"))
genes.viz = c("Dclk1", "Trpm5", "Gfi1b", "Il25") ##Tuft
FeaturePlot(SI.Seurat, genes.viz, pt.size = 0.01, cols = c("Grey","Red"))
genes.viz = c("Lgr5", "Ascl2","Slc12a2", "Axin2", "Olfm4", "Gkn3") ## Stem
FeaturePlot(SI.Seurat, genes.viz, pt.size = 0.01, cols = c("Grey","Red"))
genes.viz = c("Pcna", "Dll1", "Mki67") ##TA?
FeaturePlot(SI.Seurat, genes.viz, pt.size = 0.01, cols = c("Grey","Red"))
genes.viz = c("H2-Ab1", "H2-Aa", "Cd74", "Psmb8", "Psmb9")
FeaturePlot(SI.Seurat, genes.viz, pt.size = 0.01, cols = c("Grey","Red"))
dev.off()

print("Initial analysis done")
save.image(file="Calculations.RData")