require(ggplot2)
require(dplyr)
require(Matrix)
require(cowplot)
require(Seurat)

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
 
SI.Seurat = CreateSeuratObject(counts = SI, project = "Whole_LI", min.cells = 50, min.features = 200, meta.data = meta_data_SI) 
print("Seurat object created") 

save.image(file = "data.RData")
 