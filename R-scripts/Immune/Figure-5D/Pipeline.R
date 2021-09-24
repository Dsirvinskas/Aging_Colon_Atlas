#################################################################
############### R and package version ###########################
#################################################################
> sessionInfo()
R version 3.5.1 (2018-07-02)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS  10.14.4

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets 
[6] methods   base     

other attached packages:
 [1] Matrix_1.2-17      Hmisc_4.2-0       
 [3] Formula_1.2-3      survival_2.44-1.1 
 [5] ggforce_0.3.1      reshape2_1.4.3    
 [7] pheatmap_1.0.12    Rmisc_1.5         
 [9] plyr_1.8.4         lattice_0.20-38   
[11] tidyr_0.8.3        RColorBrewer_1.1-2
[13] ggpubr_0.2.3       magrittr_1.5      
[15] hrbrthemes_0.6.0   viridis_0.5.1     
[17] viridisLite_0.3.0  dplyr_0.8.0.1     
[19] stringr_1.4.0      Seurat_3.0.2      
[21] ggplot2_3.1.1     

loaded via a namespace (and not attached):
 [1] Rtsne_0.15          colorspace_1.4-1   
 [3] ggsignif_0.6.0      ggridges_0.5.1     
 [5] htmlTable_1.13.1    base64enc_0.1-3    
 [7] rstudioapi_0.10     listenv_0.7.0      
 [9] farver_2.0.3        npsurv_0.4-0       
[11] ggrepel_0.8.1       codetools_0.2-16   
[13] splines_3.5.1       R.methodsS3_1.7.1  
[15] extrafont_0.17      lsei_1.2-0         
[17] knitr_1.22          polyclip_1.10-0    
[19] jsonlite_1.6        Rttf2pt1_1.3.7     
[21] ica_1.0-2           cluster_2.0.8      
[23] png_0.1-7           R.oo_1.22.0        
[25] sctransform_0.2.0   compiler_3.5.1     
[27] httr_1.4.0          backports_1.1.4    
[29] assertthat_0.2.1    lazyeval_0.2.2     
[31] tweenr_1.0.1        acepack_1.4.1      
[33] htmltools_0.4.0     tools_3.5.1        
[35] rsvd_1.0.1          igraph_1.2.4.1     
[37] gtable_0.3.0        glue_1.3.1         
[39] RANN_2.6.1          Rcpp_1.0.1         
[41] gdata_2.18.0        ape_5.3            
[43] nlme_3.1-139        extrafontdb_1.0    
[45] gbRd_0.4-11         lmtest_0.9-37      
[47] xfun_0.6            globals_0.12.4     
[49] irlba_2.3.3         gtools_3.8.1       
[51] future_1.13.0       MASS_7.3-51.4      
[53] zoo_1.8-6           scales_1.0.0       
[55] parallel_3.5.1      yaml_2.2.0         
[57] reticulate_1.12     pbapply_1.4-0      
[59] gridExtra_2.3       gdtools_0.2.0      
[61] rpart_4.1-15        latticeExtra_0.6-28
[63] stringi_1.4.3       checkmate_1.9.1    
[65] caTools_1.17.1.2    bibtex_0.4.2       
[67] Rdpack_0.11-0       SDMTools_1.1-221.1 
[69] rlang_0.4.5         pkgconfig_2.0.2    
[71] systemfonts_0.1.1   bitops_1.0-6       
[73] evaluate_0.13       ROCR_1.0-7         
[75] purrr_0.3.2         htmlwidgets_1.5.1  
[77] cowplot_0.9.4       tidyselect_0.2.5   
[79] R6_2.4.0            gplots_3.0.1.1     
[81] pillar_1.3.1        foreign_0.8-71     
[83] withr_2.1.2         fitdistrplus_1.0-14
[85] nnet_7.3-12         tibble_2.1.1       
[87] future.apply_1.3.0  tsne_0.1-3         
[89] crayon_1.3.4        KernSmooth_2.23-15 
[91] plotly_4.9.0        rmarkdown_1.13     
[93] grid_3.5.1          data.table_1.12.2  
[95] metap_1.1           digest_0.6.18      
[97] R.utils_2.9.0       munsell_0.5.0   


#################################################################
############### Loading data and packages #######################
#################################################################

setwd("./")
library(Seurat)
library(stringr)
library(dplyr)
library(viridis)
library(hrbrthemes)
library(ggpubr)
library(RColorBrewer)
library(tidyr)
library(Rmisc)
library(pheatmap)
library(reshape2)
library(ggforce)
library(Hmisc)
library(Matrix)


#################################################################
############### Pre-Treatment ###################################
#################################################################
#### Sample decomplex NA
#### Do seurat

pfile<-read.csv("./data_in/p_file.csv",head=T)
immune.data <- Read10X(data.dir = "./data_in/filtered_feature_bc_matrix")
immune <- CreateSeuratObject(counts = immune.data, project = "immune", min.cells = 0, min.features = 0)
immune[["percent.mt"]] <- PercentageFeatureSet(object = immune, pattern = "^mt-")
sampo<-rownames(as.data.frame(immune$orig.ident))

pfile$type<-as.character(pfile$type)
pfile$TYPE<-as.character(pfile$TYPE)
pfile$TYPES<-as.character(pfile$TYPES)
pfile$age<-as.character(pfile$age)
pfile$pos<-as.character(pfile$pos)
pfile_sample<-pfile[match(as.character(sampo),as.character(pfile$sample_ID)),]
immune[["type"]]<-pfile_sample$type
immune[["TYPE"]]<-pfile_sample$TYPE
immune[["TYPES"]]<-pfile_sample$TYPES
immune[["age"]]<-pfile_sample$age
immune[["pos"]]<-pfile_sample$pos
immune$type[immune$type=='multible']<-'multiple'
immune$type[is.na(immune$type)]<-"not"
immune$type<-factor(immune$type,levels=c("single","double","trible","multiple","not"))
RidgePlot(immune, features = c("nFeature_RNA"),group.by="type")
RidgePlot(immune, features = c("nCount_RNA"),group.by="type")
RidgePlot(immune, features = c("percent.mt"),group.by="type")
immune_single <- subset(x = immune, subset = type == "single")
VlnPlot(object = immuneF_single, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
VlnPlot(object = immune_single, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
immuneF_single <- subset(x = immune_single, subset = percent.mt <= 8)
immuneF_single_F2 <- subset(x = immuneF_single, subset = nFeature_RNA > 200 & nFeature_RNA < 1500)
immuneF_single_F3 <- subset(x = immuneF_single_F2, subset = nCount_RNA < 5000)
VlnPlot(object = immuneF_single_F3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pfile_sample$index<-str_split_fixed(pfile_sample$sample_ID,"-",2)[,2]

for (i in 1:nrow(pfile_sample)) {
if (pfile_sample$index[i] == "3") {
if(pfile_sample$TYPES[i]=="shash1.single") {
pfile_sample$pos[i]<-"SI_P"
}
else if(pfile_sample$TYPES[i]=="shash8.single") {
pfile_sample$pos[i]<-"SI_I"
}
else if(pfile_sample$TYPES[i]=="shash3.single") {
pfile_sample$pos[i]<-"SI_D"
}
else if(pfile_sample$TYPES[i]=="shash4.single") {
pfile_sample$pos[i]<-"Cecum"
}
else if(pfile_sample$TYPES[i]=="shash5.single") {
pfile_sample$pos[i]<-"C_P"
}
else if(pfile_sample$TYPES[i]=="shash7.single") {
pfile_sample$pos[i]<-"C_D"
}
}
else if (pfile_sample$index[i] == "4") {
if(pfile_sample$TYPES[i]=="shash1.single") {
pfile_sample$pos[i]<-"SI_P"
}
else if(pfile_sample$TYPES[i]=="shash2.single") {
pfile_sample$pos[i]<-"SI_I"
}
else if(pfile_sample$TYPES[i]=="shash3.single") {
pfile_sample$pos[i]<-"SI_D"
}
else if(pfile_sample$TYPES[i]=="shash4.single") {
pfile_sample$pos[i]<-"Cecum"
}
else if(pfile_sample$TYPES[i]=="shash5.single") {
pfile_sample$pos[i]<-"C_P"
}
else if(pfile_sample$TYPES[i]=="shash6.single") {
pfile_sample$pos[i]<-"C_D"
}
}
}

sampo_F3<-rownames(as.data.frame(immuneF_single_F3$orig.ident))
immune$pos_new<-pfile_sample$pos
pfile_sample_F3<-pfile_sample[match(sampo_F3,pfile_sample$sample_ID),]
pfile_sample_F3$sample_ID<-as.character(pfile_sample_F3$sample_ID)
pfile_sample_F3$pos<-as.character(pfile_sample_F3$pos)
immuneF_single_F3$pos<-pfile_sample_F3$pos
immuneF_single_F3 <- NormalizeData(object = immuneF_single_F3, normalization.method = "LogNormalize", scale.factor = 1e4)
immuneF_single_F3 <- FindVariableFeatures(object = immuneF_single_F3,selection.method = 'mvp',mean.cutoff = c(0.0125, 2), dispersion.cutoff = c(0.9, Inf))
all.genes <- rownames(x = immuneF_single_F3)
immuneF_single_F3 <- ScaleData(object = immuneF_single_F3)
immuneF_single_F3 <- RunPCA(object = immuneF_single_F3, features = VariableFeatures(object = immuneF_single_F3))
DimPlot(object = immuneF_single_F3, reduction = 'pca')
immuneF_single_F3 <- FindNeighbors(object = immuneF_single_F3, dims = 1:20)
immuneF_single_F3 <- FindClusters(object = immuneF_single_F3, resolution = 0.7)
immuneF_single_F3 <- RunTSNE(object = immuneF_single_F3, dims = 1:20)
DimPlot(object = immuneF_single_F3, reduction = 'tsne')
immuneF_single_F3 <- FindClusters(object = immuneF_single_F3, resolution = 0.3)
immuneF_single_F3 <- RunTSNE(object = immuneF_single_F3, dims = 1:20)
DimPlot(object = immuneF_single_F3, reduction = 'tsne')
immuneF_single_F3_samples<-rownames(as.data.frame(immuneF_single_F3$orig.ident))
immuneF_single_F3[["sampleID"]]<-str_split_fixed(immuneF_single_F3_samples,"-",1)[,1]
immuneF_single_F3[["Index"]]<-str_split_fixed(immuneF_single_F3_samples,"-",2)[,2]
my.data<-FetchData(immuneF_single_F3,c("orig.ident","sampleID","Index","nCount_RNA","nFeature_RNA","percent.mt","seurat_clusters","type","TYPE","TYPES","age","pos"))
pca<-Embeddings(immuneF_single_F3, reduction = "pca")[, 1:20]
tsne<-Embeddings(immuneF_single_F3, reduction = "tsne")[, 1:2]
my.data<-merge(my.data,tsne,by="row.names",all.x=T)
rownames(my.data)<-my.data$Row.names
my.data<-my.data[,-1]
my.data<-merge(my.data,pca,by="row.names",all.x=T)
rownames(my.data)<-my.data$Row.names
my.data<-my.data[,-1]
my.data$sampleIDs<-paste(my.data$Index,my.data$age,sep="_")
colnames(my.data)[35]<-"tmp"
colnames(my.data)[2]<-"Indexs"
colnames(my.data)[3]<-"sampleID"
colnames(my.data)[2]<-"Index"
colnames(my.data)[35]<-"sampleIDs"
my.data$batch<-"batch"
for (i in 1:nrow(my.data)) {if(my.data$sampleID[i] %in% c("1","2","3","4")){my.data$batch[i]<-"batch_A"} else if(my.data$sampleID[i] %in% c("5","6")) {my.data$batch[i]<-"batch_B"} }
print(ggplot(my.data, aes(tSNE_1, tSNE_2)) +geom_point(aes(color=batch), alpha = 0.7, size = 0.3) + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank()) + facet_wrap(~batch))
immuneF_single_F3.markers <- FindAllMarkers(object = immuneF_single_F3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
immuneF_single_F3.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC) -> top20
name<-c("S_B-Cell_0","S_CD4_1","S_CD4-Naive_2","S_CD4-TEM_3","S_UNKNOWN_4","S_CD4-Cytotoxic_5","M_Monocyte-Neutrophils-Machrophages_6","S_Plasma_7","S_CD8_8","S_B-Cell-Early_9","S_CD4-Gata3_10","S_CD4-Naive-Isg_11","S_NK_12","S_pDCs_13")
my.data$celltypes<-name[my.data$seurat_clusters]
my.data.new<-data.frame(table(my.data[,c(37,12,11,35)])) ## celltypes,pos,age,sampleIDs
my.data$celltypes<-factor(my.data$celltypes,levels = names(num[order(num,decreasing = TRUE)]))
my.data$age<-factor(my.data$age,levels = c("young","old"))
immuneF_single_F3$clusters_age<-paste(immuneF_single_F3$seurat_clusters,immuneF_single_F3$age,sep=".")
clusters<-unique(immuneF_single_F3$seurat_clusters)
for(i in clusters){assign(paste("age.DEG.wilcox",i,sep = "_") , FindMarkers(immuneF_single_F3, ident.1 = paste(i,"old",sep="."), ident.2 = paste(i,"young",sep="."), verbose = FALSE,group.by='clusters_age'))}
for(i in clusters){write.csv(get(paste("age.DEG.wilcox",i,sep = "_")),paste(i,"age.DEG.wilcox.csv",sep = "_"))}
immuneF_single_F3$celltypes<-name[immuneF_single_F3$seurat_clusters]
immuneF_single_F3$age<-factor(immuneF_single_F3$age,levels=c("young","old"))
immuneF_single_F3$Index_age<-paste(immuneF_single_F3$Index,immuneF_single_F3$age,sep="_")

saveRDS(immuneF_single_F3,"immuneF_single_F3")

#################################################################
############### Atlas.pipALL ####################################
#################################################################

myData <- readRDS(file = "./data_in/immuneF_single_F3")
#### subset large intestine
myData_LI<-subset(myData, subset = pos %in% c("C_D","Cecum","C_P"))
rm(myData)
myData_LI$seurat_clusters<-NULL
myData_LI$RNA_snn_res.0.7<-NULL
myData_LI$RNA_snn_res.0.3<-NULL
myData_LI <- NormalizeData(object = myData_LI, normalization.method = "LogNormalize", scale.factor = 1e4)
myData_LI <- FindVariableFeatures(object = myData_LI,selection.method = 'mvp',mean.cutoff = c(0.0125, 2), dispersion.cutoff = c(0.9, Inf))
all.genes <- rownames(x = myData_LI)
myData_LI <- ScaleData(object = myData_LI)
myData_LI <- RunPCA(object = myData_LI, features = VariableFeatures(object = myData_LI))
DimPlot(object = myData_LI, reduction = 'pca')
myData_LI <- FindNeighbors(object = myData_LI, dims = 1:20)
#### res=0.4
myData_LI <- FindClusters(object = myData_LI, resolution = 0.4)
myData_LI <- RunTSNE(object = myData_LI, dims = 1:20)
DimPlot(object = myData_LI, reduction = 'tsne')
DimPlot(object = myData_LI, reduction = 'tsne',split.by = "age")+coord_fixed(ratio = 1)
DimPlot(object = myData_LI, reduction = 'tsne',split.by = "age",label = T)+coord_fixed(ratio = 1)
#### Calculate markers
myData_LI.markers <- FindAllMarkers(object = myData_LI, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(myData_LI.markers,"./data_out/myData_LI.markers.csv")
#### Compare with previous markers top20
markers_all<-read.csv("/Users/jlu/Desktop/Pro_SingleCell/Data_Hash/Hashtag7/immuneF_single_F3.markers.csv",header = T)
myData_LI.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC) -> top20_target
markers_all %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC) -> top20_reference
top20_target$cluster<-paste("t",top20_target$cluster,sep="_")
top20_reference$cluster<-paste("r",top20_reference$cluster,sep="_")
clusters_celltypes_top20<-merge(top20_target,top20_reference,by.x='gene',by.y='gene',all=T)
nb.cols <- 14
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)
myData_LI.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) -> top50_target
markers_all %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) -> top50_reference
top50_target$cluster<-paste("t",top50_target$cluster,sep="_")
top50_reference$cluster<-paste("r",top50_reference$cluster,sep="_")
clusters_celltypes_top50<-merge(top50_target,top50_reference,by.x='gene',by.y='gene',all=T)

write.csv(top20_target,"./data_out//markers_top20.csv")
write.csv(top50_target,"./data_out/markers_top50.csv")

celltypes<-c("B-cell","CD4_CD8_CD4-Cytotoxic","CD4-TEM_CD4-Gata3_CD4-Cytotoxic","UNKNOWN_CD4-TEM_NK_CD4-Gata3","CD4-Naive_CD4-Cytotoxic","Plasma","CD8_NK_CD4","ct7","Myeloblast","ct9","Myeloblast","CD4-Naive-Isg","B-cell-early")
celltypes[10]
myData_LI$celltypes_new<-celltypes[myData_LI$seurat_clusters]
myData_LI$celltypes_old<-myData_LI$celltypes
myData_LI$celltypes<-NULL
my.data<-FetchData(myData_LI,c("orig.ident","sampleID","Index","nCount_RNA","nFeature_RNA","percent.mt","seurat_clusters","type","TYPE","TYPES","age","pos","celltypes_old","celltypes_new"))
pca<-Embeddings(myData_LI, reduction = "pca")[, 1:20]
tsne<-Embeddings(myData_LI, reduction = "tsne")[, 1:2]
my.data<-merge(my.data,tsne,by="row.names",all.x=T)
rownames(my.data)<-my.data$Row.names
my.data<-my.data[,-1]
my.data<-merge(my.data,pca,by="row.names",all.x=T)
rownames(my.data)<-my.data$Row.names
my.data<-my.data[,-1]
write.csv(my.data,"./data_out//my.data.csv")
my.data$celltypes_new_2<-as.character(my.data$celltypes_new)
my.data[my.data$celltypes_new_2=="CD4_CD8_CD4-Cytotoxic","celltypes_new_2"]<-"CD4"
my.data[my.data$celltypes_new_2=="CD4-Naive_CD4-Cytotoxic","celltypes_new_2"]<-"CD4-Naive"
my.data[my.data$celltypes_new_2=="CD4-TEM_CD4-Gata3_CD4-Cytotoxic","celltypes_new_2"]<-"CD4-TEM"
my.data[my.data$celltypes_new_2=="CD8_NK_CD4","celltypes_new_2"]<-"CD8"
my.data[my.data$celltypes_new_2=="ct7","celltypes_new_2"]<-"CD4_2"
my.data[my.data$celltypes_new_2=="ct9","celltypes_new_2"]<-"B-cells-2"
my.data[my.data$celltypes_new_2=="CD4_2","celltypes_new_2"]<-"CD4-2"
my.data[my.data$celltypes_new_2=="UNKNOWN_CD4-TEM_NK_CD4-Gata3","celltypes_new_2"]<-"Th17"
ggplot(my.data, aes(celltypes_new_2,fill=celltypes_old))+geom_bar()+scale_fill_manual(values = mycolors)+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=10),legend.text = element_text(size = 10),legend.title = element_text())

#### plot celltypes ratio

my.data2<-my.data[,c("celltypes_new_2", "age","Index")]
my.data2$var<-paste(my.data2$celltypes_new_2,my.data2$age,my.data2$Index,sep="#")
my.data2<-as.data.frame(table(my.data2$var))
my.data2$celltypes_new_2<-str_split_fixed(my.data2$Var1,"#",3)[,1]
my.data2$age<-str_split_fixed(my.data2$Var1,"#",3)[,2]
my.data2$Index<-str_split_fixed(my.data2$Var1,"#",3)[,3]
my.data2<-as.data.frame(my.data2)
tg<-my.data2 %>% dplyr::group_by(age,Index)  %>% dplyr::mutate(ratio = Freq/sum(Freq))
tg$Var2<-paste(tg$Index,tg$age,sep="#")
sum(tg[tg$Var2=="5#old","Freq"])
tg[tg$Var2=="5#old",]
tg$age<-factor(tg$age,levels = c("young","old"))
p<-ggbarplot(tg, x = "celltypes_new_2", y = "ratio", add = "mean_sd",
color = "age", palette = "jco",
position = position_dodge(0.8))
p+stat_compare_means(aes(group = age),method = "t.test", label = "p.format",size=2)+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+ggtitle("t.test_pairedF")+scale_color_manual(values = c("grey","grey0"))
#### check 1
my.data2 %>% dplyr::group_by(age,Index)  %>% dplyr::mutate(x = sum(Freq))
myData_LI <- RunUMAP(object = myData_LI, dims = 1:20)
DimPlot(object = myData_LI, reduction = 'umap',split.by = "Index_age",ncol =4)+coord_fixed(ratio = 1)
myData_LI <- FindClusters(object = myData_LI, resolution = 0.1)
DimPlot(object = myData_LI, reduction = 'tsne',group.by = "RNA_snn_res.0.1")+coord_fixed(ratio = 1)
#### check 2
my.data.2<-my.data[,c("celltypes_old", "age","Index")]
my.data.2$var<-paste(my.data.2$celltypes_old,my.data.2$age,my.data.2$Index,sep="#")
my.data.2<-as.data.frame(table(my.data.2$var))
my.data.2$celltypes_old_2<-str_split_fixed(my.data.2$Var1,"#",3)[,1]
my.data.2$age<-str_split_fixed(my.data.2$Var1,"#",3)[,2]
my.data.2$Index<-str_split_fixed(my.data.2$Var1,"#",3)[,3]
my.data.2<-as.data.frame(my.data.2)
tg.2<-my.data.2 %>% dplyr::group_by(age,Index)  %>% dplyr::mutate(ratio = Freq/sum(Freq))
tg.2$Var2<-paste(tg.2$Index,tg.2$age,sep="#")
tg.2$age<-factor(tg.2$age,levels = c("young","old"))
ggbarplot(tg.2, x = "celltypes_old_2", y = "ratio", add = "mean_sd",
color = "age", palette = "jco",
position = position_dodge(0.8))+stat_compare_means(aes(group = age),method = "t.test", label = "p.format",size=2)+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+ggtitle("t.test_pairedF")+scale_color_manual(values = c("grey","grey0"))

write.csv(tg,"./data_out/tg.csv")
write.csv(tg.2,"./data_out/tg2.csv")

ggbarplot(tg[(tg$Var2!="6#old")&(tg$Var2!="6#young"),], x = "celltypes_new_2", y = "ratio", add = "mean_sd",
color = "age", palette = "jco",
position = position_dodge(0.8))+stat_compare_means(aes(group = age),method = "t.test", label = "p.format",size=2)+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+ggtitle("t.test_pairedF")+scale_color_manual(values = c("grey","grey0"))

#### show pos and age and sample in facet plot
my.data3<-my.data[,c("celltypes_new_2", "age","Index","pos")]
my.data3$var<-paste(my.data3$celltypes_new_2,my.data3$age,my.data3$Index,my.data3$pos,sep="#")
my.data3<-as.data.frame(table(my.data3$var))
my.data3$celltypes_new_2<-str_split_fixed(my.data3$Var1,"#",4)[,1]
my.data3$age<-str_split_fixed(my.data3$Var1,"#",4)[,2]
my.data3$Index<-str_split_fixed(my.data3$Var1,"#",4)[,3]
my.data3$pos<-str_split_fixed(my.data3$Var1,"#",4)[,4]
my.data3<-as.data.frame(my.data3)
tg3<-my.data3 %>% dplyr::group_by(age,Index,pos)  %>% dplyr::mutate(ratio = Freq/sum(Freq))
tg3$Var2<-paste(tg3$Index,tg3$age,sep="#")
mycolors <- colorRampPalette(brewer.pal(12,"Paired"))(12)
tg3$Var2<-factor(tg3$Var2,levels = c("3#young","4#young","5#young","6#young","1#old","2#old","5#old","6#old"))
tg3$pos<-factor(tg3$pos,levels = c("Cecum","C_P","C_D"))
ggplot(tg3, aes(fill=celltypes_new_2, y=ratio, x=pos)) +
geom_bar(position="fill", stat="identity")+scale_fill_manual(values = mycolors)+facet_wrap(~Var2,nrow=2)+
theme(strip.text.x = element_text(size = 10),axis.text.x = element_text(angle=90, hjust=1))
myData_LI <- FindClusters(object = myData_LI, resolution = 0.7)
myData_LI <- FindClusters(object = myData_LI, resolution = 0.8)
DimPlot(object = myData_LI, reduction = 'tsne',group.by="RNA_snn_res.0.8",split.by = "age",label = T)+coord_fixed(ratio = 1)
DimPlot(object = myData_LI, reduction = 'tsne',group.by="RNA_snn_res.0.7",split.by = "age",label = T)+coord_fixed(ratio = 1)
#### redo with new resolution
myData_LI$seurat_clusters<-myData_LI$RNA_snn_res.0.7
x<-as.data.frame(myData_LI$RNA_snn_res.0.7)
#### merge cluster8 to cluster0
x[,1]<-as.character(x[,1])
x[x[,1]==8,1]<-0
x[,1]<-factor(x[,1],levels = c("0","1","2","3","4","5","6","7","9","10","11","12","13","14","15","16"))
myData_LI$seurat_clusters<-x[,1]
DimPlot(object = myData_LI, reduction = 'tsne',group.by="seurat_clusters",split.by = "age",label = T)+coord_fixed(ratio = 1)
myData_LI$RNA_snn_res.0.8<-NULL
myData_LI <- FindClusters(object = myData_LI, resolution = 0.7)
DimPlot(object = myData_LI, reduction = 'tsne',split.by = "age",label = T)+coord_fixed(ratio = 1)
myData_LI$seurat_clusters<-x[,1]
#### markers of original resolution 0.7 
myData_LI.markers.redo <- FindAllMarkers(object = myData_LI, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(myData_LI.markers.redo,"./data_out//myData_LI.markers.redo.csv")
#### markers of merged cluster8 and cluster0
Idents(myData_LI)<-x[,1]
myData_LI.markers.redo2 <- FindAllMarkers(object = myData_LI, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(myData_LI.markers.redo2,"./data_out//myData_LI.markers.redo2.csv")

#### top 20 markers of this dataset with top20 markers of all dataset
myData_LI.markers.redo %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC) -> top20_target.redo
top20_target.redo$cluster<-paste("t",top20_target.redo$cluster,sep="_")
clusters_celltypes_top20.redo<-merge(top20_target.redo,top20_reference,by.x='gene',by.y='gene',all=T)
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(16)
ggplot(clusters_celltypes_top20.redo, aes(cluster.x))+geom_bar(aes(fill=cluster.y))+scale_fill_manual(values = mycolors)
#### top50
myData_LI.markers.redo %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) -> top50_target.redo
top50_target.redo$cluster<-paste("t",top50_target.redo$cluster,sep="_")
clusters_celltypes_top50.redo<-merge(top50_target.redo,top50_reference,by.x='gene',by.y='gene',all=T)
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(16)
ggplot(clusters_celltypes_top50.redo, aes(cluster.x))+geom_bar(aes(fill=cluster.y))+scale_fill_manual(values = mycolors)

celltypes.redo<-c("B-cell","CD4-1","CD4-TEM","CD8","CD4-2","Plasma","CD4-Naive","CD4-3","B-cell","Myeloblast_9","B-cell","Th17","CD4-Cytotoxic","CD4-Naive-Isg","Myeloblast_14","B-cell-early","Myeloblast_16")
immune.markers<-read.csv("/Users/jlu/Desktop/Pro_SingleCell/Data_Hash/Hashtag7/immune_markers_more.csv",head=T)
tmp<-merge(top20_target.redo,immune.markers,by.x = "gene",by.y = "Gene")
tmp[tmp$cluster=="t_14",]
tmp[tmp$cluster=="t_16",]
tmp[tmp$cluster=="t_9",]
tmp<-merge(top50_target.redo,immune.markers,by.x = "gene",by.y = "Gene")
tmp[tmp$cluster=="t_9",]
tmp[tmp$cluster=="t_14",]
tmp[tmp$cluster=="t_16",]
celltypes.redo<-c("B-cell-1","CD4-1","CD4-TEM","CD8","CD4-2","Plasma","CD4-Naive","CD4-3","B-cell-1","Neutrophil","B-cell-2","Th17","CD4-Cytotoxic","CD4-Naive-Isg","Monocyte","B-cell-early","Macrophage")
myData_LI$celltypes_new_redo<-celltypes.redo[myData_LI$RNA_snn_res.0.7]
my.data.redo<-FetchData(myData_LI,c("orig.ident","sampleID","Index","nCount_RNA","nFeature_RNA","percent.mt","seurat_clusters","type","TYPE","TYPES","age","pos","celltypes_old","celltypes_new","celltypes_new_redo"))
pca<-Embeddings(myData_LI, reduction = "pca")[, 1:20]
tsne<-Embeddings(myData_LI, reduction = "tsne")[, 1:2]
my.data.redo<-merge(my.data.redo,tsne,by="row.names",all.x=T)
rownames(my.data.redo)<-my.data.redo$Row.names
my.data.redo<-my.data.redo[,-1]
my.data.redo<-merge(my.data.redo,pca,by="row.names",all.x=T)
rownames(my.data.redo)<-my.data.redo$Row.names
my.data.redo<-my.data.redo[,-1]
write.csv(my.data.redo,"./data_out//my.data.redo.csv")
ggplot(my.data.redo, aes(celltypes_new_redo,fill=celltypes_old))+geom_bar()+scale_fill_manual(values = mycolors)+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=10),legend.text = element_text(size = 10),legend.title = element_text())

#### plot celltypes ratio in each sample-compartment
facet<-function(my.data){
my.data3<-my.data[,c("celltypes_new_redo","age","Index","pos")]
my.data3$var<-paste(my.data3$celltypes_new_redo,my.data3$age,my.data3$Index,my.data3$pos,sep="#")
my.data3<-as.data.frame(table(my.data3$var))
my.data3$celltypes_new_redo<-str_split_fixed(my.data3$Var1,"#",4)[,1]
my.data3$age<-str_split_fixed(my.data3$Var1,"#",4)[,2]
my.data3$Index<-str_split_fixed(my.data3$Var1,"#",4)[,3]
my.data3$pos<-str_split_fixed(my.data3$Var1,"#",4)[,4]
my.data3<-as.data.frame(my.data3)
tg3<-my.data3 %>% dplyr::group_by(age,Index,pos)  %>% dplyr::mutate(ratio = Freq/sum(Freq))
tg3$Var2<-paste(tg3$Index,tg3$age,sep="#")
mycolors <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(my.data$celltypes_new_redo)))
tg3$Var2<-factor(tg3$Var2,levels = c("3#young","4#young","5#young","6#young","1#old","2#old","5#old","6#old"))
tg3$pos<-factor(tg3$pos,levels = c("Cecum","C_P","C_D"))
ggplot(tg3, aes(fill=celltypes_new_redo, y=ratio, x=pos)) +
geom_bar(position="fill", stat="identity")+scale_fill_manual(values = mycolors)+facet_wrap(~Var2,nrow=2)+
theme(strip.text.x = element_text(size = 10),axis.text.x = element_text(angle=90, hjust=1))
}
facet(my.data.redo)
#### plot celltypes ratio in each age-sample-compartment
facet<-function(my.data){
my.data3<-my.data[,c("celltypes_new_redo","age","Index","pos")]
my.data3$var<-paste(my.data3$celltypes_new_redo,my.data3$age,my.data3$Index,my.data3$pos,sep="#")
my.data3<-as.data.frame(table(my.data3$var))
my.data3$celltypes_new_redo<-str_split_fixed(my.data3$Var1,"#",4)[,1]
my.data3$age<-str_split_fixed(my.data3$Var1,"#",4)[,2]
my.data3$Index<-str_split_fixed(my.data3$Var1,"#",4)[,3]
my.data3$pos<-str_split_fixed(my.data3$Var1,"#",4)[,4]
my.data3<-as.data.frame(my.data3)
tg3<-my.data3 %>% dplyr::group_by(age,Index,pos)  %>% dplyr::mutate(ratio = Freq/sum(Freq))
tg3$Var2<-paste(tg3$Index,tg3$age,sep="#")
library(RColorBrewer)
mycolors <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(my.data$celltypes_new_redo)))
tg3$Var2<-factor(tg3$Var2,levels = c("3#young","4#young","5#young","6#young","1#old","2#old","5#old","6#old"))
tg3$pos<-factor(tg3$pos,levels = c("Cecum","C_P","C_D"))
ggplot(tg3, aes(fill=celltypes_new_redo, y=ratio, x=pos)) +
geom_bar(position="fill", stat="identity")+scale_fill_manual(values = mycolors)+facet_wrap(~Var2,nrow=2)+
theme(strip.text.x = element_text(size = 10),axis.text.x = element_text(angle=90, hjust=1))
return(tg3)
}
tg4<-facet(my.data.redo)
tg4$age<-factor(tg4$age,levels = c("young","old"))
#### ggpubr did not work so i use ggplot instead
mus2 <- summarySE(tg4, measurevar="ratio",
groupvars=c("celltypes_new_redo", "pos", "age"), na.rm = TRUE)
#### bar
ggplot(mus2, aes(x=pos, y=ratio, fill=age)) +
geom_bar(stat="identity", colour="black",position="dodge") +
facet_wrap(~celltypes_new_redo,nrow=4) + geom_errorbar(aes(ymin=ratio-sd, ymax=ratio+sd), size=0.5, width=.25,position=position_dodge(.9)) +ggtitle("mean_sd_t.test")+theme_classic()+scale_fill_manual(values = c("gray","gray40"))+stat_compare_means(data=tg4,aes(group = age), label = "p.signif",method="t.test")
#### line
ggplot(mus2, aes(x=pos, y=ratio, group=age,color=age)) +
geom_line() + geom_point()+
facet_wrap(~celltypes_new_redo,nrow=4) + geom_errorbar(aes(ymin=ratio-sd, ymax=ratio+sd), size=0.5, width=.1, position=position_dodge(0.05)) +ggtitle("mean_sd_t.test")+theme_classic()+scale_color_manual(values = c("gray","gray0"))+stat_compare_means(data=tg4,aes(group = age), label = "p.signif",method="t.test")
#### odds ratio
#### general all positions
tg2<-tg4 %>% dplyr::group_by(celltypes_new_redo,age)  %>% dplyr::mutate(celltypes_sum = sum(Freq))
tg2<-unique(tg2[,c(3,4,9)])
tg2<-tg2 %>% dplyr::group_by(age)  %>% dplyr::mutate(celltypes_sum_ratio = celltypes_sum/sum(celltypes_sum))
tg2_wide<-dcast(tg2,celltypes_new_redo~tg2$age,value.var = 'celltypes_sum_ratio')
tg2_wide$odds_ratio<-(tg2_wide$old/(1-tg2_wide$old))/(tg2_wide$young/(1-tg2_wide$young))
tg2<-tg2 %>% dplyr::group_by(age) %>% dplyr::mutate(age_sum = sum(celltypes_sum))
tg2_wide2<-dcast(tg2,celltypes_new_redo~tg2$age,value.var = c("celltypes_sum"))
tg2_wide2_2<-dcast(tg2,celltypes_new_redo~tg2$age,value.var = c("age_sum"))
tg2_wide2<-merge(tg2_wide2,tg2_wide2_2,by="celltypes_new_redo",suffixes=c(".specific",".total"))
tg2_wide2$k<-(tg2_wide2$young.specific+tg2_wide2$old.specific)
fisher_conf<-function(q,m,n,x){tmp<-fisher.test(matrix(c(q,(m-q),x,(n-x)),nrow=2))
return(tmp$conf.int[c(1,2)])}
#### p value
colnames(tg2_wide2)<-c("celltypes_new_redo","x_young.specific","q_old.specific","n_young.total","m_old.total","k")
for (i in 1:nrow(tg2_wide2)){tg2_wide2$conf_1[i]<-fisher_conf(tg2_wide2[i,"q_old.specific"],tg2_wide2[i,"m_old.total"],tg2_wide2[i,"n_young.total"],tg2_wide2[i,"x_young.specific"])[1];tg2_wide2$conf_2[i]<-fisher_conf(tg2_wide2[i,"q_old.specific"],tg2_wide2[i,"m_old.total"],tg2_wide2[i,"n_young.total"],tg2_wide2[i,"x_young.specific"])[2]}
tg2_wide2_m<-merge(tg2_wide,tg2_wide2,by = "celltypes_new_redo",all = T)
####95% CI
tg2_wide2_m$hypergeometric_p_lower<-phyper(tg2_wide2$q_old.specific,tg2_wide2$m_old.total,tg2_wide2$n_young.total,tg2_wide2$k)
tg2_wide2_m$hypergeometric_p_upper<-phyper(tg2_wide2$q_old.specific,tg2_wide2$m_old.total,tg2_wide2$n_young.total,tg2_wide2$k,lower.tail = F)
tg2_wide2_m$hypergeometric_p_tails<-pmin(tg2_wide2_m$hypergeometric_p_lower,tg2_wide2_m$hypergeometric_p_upper)
ggplot(tg2_wide2_m, aes(x=celltypes_new_redo, y=odds_ratio)) +
geom_bar(stat="identity", color="black") +
geom_errorbar(aes(ymin=conf_1, ymax=conf_2), width=.2,
position=position_dodge(.9)) + coord_flip()+scale_x_discrete(limits = rev(levels(factor(tg2_wide2_m$celltypes_new_redo))))+theme_classic(base_size = 20) + geom_hline(yintercept=1, linetype="dashed", color = "black")+geom_text(aes(label=signif(hypergeometric_p_tails,digits = 3)),size = 3, position=position_dodge(width=1),hjust = -1,vjust=2)
#### other positions
odds_ratio<-function(tg,pos){
#### odds ratio
tg<-tg[tg$pos==pos,]
tg2<-tg %>% dplyr::group_by(celltypes_new_redo,age)  %>% dplyr::mutate(celltypes_sum = sum(Freq))
tg2<-unique(tg2[,c(3,4,9)])
tg2<-tg2 %>% dplyr::group_by(age)  %>% dplyr::mutate(celltypes_sum_ratio = celltypes_sum/sum(celltypes_sum))
require(reshape2)
tg2_wide<-dcast(tg2,celltypes_new_redo~tg2$age,value.var = 'celltypes_sum_ratio')
tg2_wide$odds_ratio<-(tg2_wide$old/(1-tg2_wide$old))/(tg2_wide$young/(1-tg2_wide$young))
#### p value by fisher exact test
tg2<-tg2 %>% dplyr::group_by(age) %>% dplyr::mutate(age_sum = sum(celltypes_sum))
tg2_wide2<-dcast(tg2,celltypes_new_redo~tg2$age,value.var = c("celltypes_sum"))
tg2_wide2_2<-dcast(tg2,celltypes_new_redo~tg2$age,value.var = c("age_sum"))
tg2_wide2<-merge(tg2_wide2,tg2_wide2_2,by="celltypes_new_redo",suffixes=c(".specific",".total"))
tg2_wide2$k<-(tg2_wide2$young.specific+tg2_wide2$old.specific)
colnames(tg2_wide2)<-c("celltypes_new_redo","x_young.specific","q_old.specific","n_young.total","m_old.total","k")
for (i in 1:nrow(tg2_wide2)){tg2_wide2$conf_1[i]<-fisher_conf(tg2_wide2[i,"q_old.specific"],tg2_wide2[i,"m_old.total"],tg2_wide2[i,"n_young.total"],tg2_wide2[i,"x_young.specific"])[1];tg2_wide2$conf_2[i]<-fisher_conf(tg2_wide2[i,"q_old.specific"],tg2_wide2[i,"m_old.total"],tg2_wide2[i,"n_young.total"],tg2_wide2[i,"x_young.specific"])[2]}
tg2_wide2_m<-merge(tg2_wide,tg2_wide2,by = "celltypes_new_redo",all = T)
#### 95%CI value hypergeometric test
tg2_wide2_m$hypergeometric_p_lower<-phyper(tg2_wide2$q_old.specific,tg2_wide2$m_old.total,tg2_wide2$n_young.total,tg2_wide2$k)
tg2_wide2_m$hypergeometric_p_upper<-phyper(tg2_wide2$q_old.specific,tg2_wide2$m_old.total,tg2_wide2$n_young.total,tg2_wide2$k,lower.tail = F)
tg2_wide2_m$hypergeometric_p_tails<-pmin(tg2_wide2_m$hypergeometric_p_lower,tg2_wide2_m$hypergeometric_p_upper)
require(ggplot2)
ggplot(tg2_wide2_m, aes(x=celltypes_new_redo, y=odds_ratio)) +
geom_bar(stat="identity", color="black") +
geom_errorbar(aes(ymin=conf_1, ymax=conf_2), width=.2,
position=position_dodge(.9)) + coord_flip()+scale_x_discrete(limits = rev(levels(factor(tg2_wide2_m$celltypes_new_redo))))+theme_classic(base_size = 20) + geom_hline(yintercept=1, linetype="dashed", color = "black")+geom_text(aes(label=signif(hypergeometric_p_tails,digits = 3)),size = 3, position=position_dodge(width=1),hjust = -1,vjust=2)
}
odds_ratio(tg4,"Cecum")
odds_ratio(tg4,"C_P")
odds_ratio(tg4,"C_D")

################### redo C_D
odds_ratio<-function(tg,pos){
#### odds ratio
tg<-tg[tg$pos==pos,]
tg2<-tg %>% dplyr::group_by(celltypes_new_redo,age)  %>% dplyr::mutate(celltypes_sum = sum(Freq))
tg2<-unique(tg2[,c(3,4,9)])
tg2<-tg2 %>% dplyr::group_by(age)  %>% dplyr::mutate(celltypes_sum_ratio = celltypes_sum/sum(celltypes_sum))
require(reshape2)
tg2_wide<-dcast(tg2,celltypes_new_redo~tg2$age,value.var = 'celltypes_sum_ratio')
tg2_wide$odds_ratio<-(tg2_wide$old/(1-tg2_wide$old))/(tg2_wide$young/(1-tg2_wide$young))
#### p value by fisher exact test
tg2<-tg2 %>% dplyr::group_by(age) %>% dplyr::mutate(age_sum = sum(celltypes_sum))
tg2_wide2<-dcast(tg2,celltypes_new_redo~tg2$age,value.var = c("celltypes_sum"))
tg2_wide2_2<-dcast(tg2,celltypes_new_redo~tg2$age,value.var = c("age_sum"))
tg2_wide2<-merge(tg2_wide2,tg2_wide2_2,by="celltypes_new_redo",suffixes=c(".specific",".total"))
tg2_wide2$k<-(tg2_wide2$young.specific+tg2_wide2$old.specific)
colnames(tg2_wide2)<-c("celltypes_new_redo","x_young.specific","q_old.specific","n_young.total","m_old.total","k")
for (i in 1:nrow(tg2_wide2)){tg2_wide2$conf_1[i]<-fisher_conf(tg2_wide2[i,"q_old.specific"],tg2_wide2[i,"m_old.total"],tg2_wide2[i,"n_young.total"],tg2_wide2[i,"x_young.specific"])[1];tg2_wide2$conf_2[i]<-fisher_conf(tg2_wide2[i,"q_old.specific"],tg2_wide2[i,"m_old.total"],tg2_wide2[i,"n_young.total"],tg2_wide2[i,"x_young.specific"])[2]}
tg2_wide2_m<-merge(tg2_wide,tg2_wide2,by = "celltypes_new_redo",all = T)
#### 95%CI value hypergeometric test
tg2_wide2_m$hypergeometric_p_lower<-phyper(tg2_wide2$q_old.specific,tg2_wide2$m_old.total,tg2_wide2$n_young.total,tg2_wide2$k)
tg2_wide2_m$hypergeometric_p_upper<-phyper(tg2_wide2$q_old.specific,tg2_wide2$m_old.total,tg2_wide2$n_young.total,tg2_wide2$k,lower.tail = F)
tg2_wide2_m$hypergeometric_p_tails<-pmin(tg2_wide2_m$hypergeometric_p_lower,tg2_wide2_m$hypergeometric_p_upper)
require(ggplot2)
ggplot(tg2_wide2_m, aes(x=celltypes_new_redo, y=odds_ratio)) +
geom_bar(stat="identity", color="black") +
geom_errorbar(aes(ymin=conf_1, ymax=conf_2), width=.2,
position=position_dodge(.9)) + coord_flip()+scale_x_discrete(limits = rev(levels(factor(tg2_wide2_m$celltypes_new_redo))))+theme_classic(base_size = 20) + geom_hline(yintercept=1, linetype="dashed", color = "black")+geom_text(aes(label=signif(hypergeometric_p_tails,digits = 3)),size = 3, position=position_dodge(width=1),hjust = -1,vjust=2)
return(tg2_wide2_m)
}

cd<-odds_ratio(tg4,"C_D")

p1<-ggplot(cd,aes(x=celltypes_new_redo, y=odds_ratio)) +
geom_bar(stat="identity", color="black") +
geom_errorbar(aes(ymin=conf_1, ymax=conf_2), width=.2,
position=position_dodge(.9)) + scale_x_discrete(limits = rev(levels(factor(tg2_wide2_m$celltypes_new_redo))))+theme_classic(base_size = 20) + geom_hline(yintercept=1, linetype="dashed", color = "black")+geom_text(aes(label=signif(hypergeometric_p_tails,digits = 3)),size = 3, position=position_dodge(width=1),hjust = -1,vjust=2)+coord_flip()

p2<-ggplot(cd,aes(x=celltypes_new_redo, y=odds_ratio)) +
geom_bar(stat="identity", color="black") +
geom_errorbar(aes(ymin=conf_1, ymax=conf_2), width=.2,
position=position_dodge(.9)) + scale_x_discrete(limits = rev(levels(factor(tg2_wide2_m$celltypes_new_redo))))+theme_classic(base_size = 20) + geom_hline(yintercept=1, linetype="dashed", color = "black")+geom_text(aes(label=signif(hypergeometric_p_tails,digits = 3)),size = 3, position=position_dodge(width=1),hjust = -1,vjust=2)+coord_flip()+ylim(0,10)
cowplot::plot_grid(p1, p2)
#### adjust color order in tsne, make it consistent with other figures
x<-myData_LI$seurat_clusters
x<-as.data.frame(myData_LI$seurat_clusters)
o<-c("0","3","9","10","4","14","7","5","13","1","15","6","8","12","2","11")
x$color_order<-o[x[,1]]
x[,2]<-factor(x[,2],levels = 0:15)
myData_LI$seurat_clusters_ColorOrder<-x$color_order
DimPlot(object = myData_LI, reduction = 'tsne',split.by = "age",group.by="seurat_clusters_ColorOrder",label = F,cols =mycolors)+coord_fixed(ratio = 1)
#### correct colors
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(16)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(16)
ggplot(tg4, aes(fill=celltypes_new_redo, y=ratio, x=pos)) +
geom_bar(position="fill", stat="identity")+scale_fill_manual(values = mycolors)+facet_wrap(~Var2,nrow=2)+
theme(strip.text.x = element_text(size = 10),axis.text.x = element_text(angle=90, hjust=1))
DimPlot(object = myData_LI, reduction = 'tsne',split.by = "age",group.by="seurat_clusters_ColorOrder",label = F,cols =mycolors)+coord_fixed(ratio = 1)
DimPlot(object = myData_LI, reduction = 'tsne',split.by = "age",group.by="seurat_clusters_ColorOrder",label = T,cols =mycolors)+coord_fixed(ratio = 1)




#### REDO2
#### redo1 merge cluster8 to B-cell
#### redo2 merge cluster B-cell and CD4
tmp[tmp$cluster=="t_7",]
tmp<-merge(top50_target.redo,immune.markers,by.x = "gene",by.y = "Gene")
tmp[tmp$cluster=="t_7",]
tmp[tmp$cluster=="t_15",]
tmp[tmp$cluster=="t_13",]
mycolors2 <- colorRampPalette(brewer.pal(12, "Paired"))(14)
celltypes.redo2<-c("B-cell","B-cell","Bcell-Early","CD4","CD4","CD4-Early","CD4-Cytotoxic","CD4-Naive","Bcell-Naive-Isg","CD4-TEM","CD8","Macrophage","Monocyte","Neutrophil","Plasma","Th17")
myData_LI$celltypes_new_redo2<-celltypes.redo2[myData_LI$seurat_clusters_ColorOrder]
myData_LI$celltypes_new_redo2<-factor(myData_LI$celltypes_new_redo2,levels = c("CD4-Early","CD4-Naive","CD4","CD4-Cytotoxic","CD4-TEM","Th17","CD8","Bcell-Early","Bcell-Naive-Isg","B-cell","Plasma","Neutrophil","Monocyte","Macrophage"))
#### TSNE
DimPlot(object = myData_LI, reduction = 'tsne',split.by = "age",group.by="celltypes_new_redo2",label = T,cols =mycolors2)+coord_fixed(ratio = 1)
DimPlot(object = myData_LI, reduction = 'tsne',split.by = "age",group.by="celltypes_new_redo2",label = F,cols =mycolors2)+coord_fixed(ratio = 1)
#### compare with old celltypes
my.data.redo2<-FetchData(myData_LI,c("orig.ident","sampleID","Index","nCount_RNA","nFeature_RNA","percent.mt","seurat_clusters","type","TYPE","TYPES","age","pos","celltypes_old","celltypes_new","celltypes_new_redo","celltypes_new_redo2"))
pca<-Embeddings(myData_LI, reduction = "pca")[, 1:20]
tsne<-Embeddings(myData_LI, reduction = "tsne")[, 1:2]
my.data.redo2<-merge(my.data.redo2,tsne,by="row.names",all.x=T)
rownames(my.data.redo2)<-my.data.redo2$Row.names
my.data.redo2<-my.data.redo2[,-1]
my.data.redo2<-merge(my.data.redo2,pca,by="row.names",all.x=T)
rownames(my.data.redo2)<-my.data.redo2$Row.names
my.data.redo2<-my.data.redo2[,-1]
write.csv(my.data.redo2,"./data_out//my.data.redo2.csv")
ggplot(my.data.redo2, aes(celltypes_new_redo2,fill=celltypes_old))+geom_bar()+scale_fill_manual(values = mycolors2)+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=10),legend.text = element_text(size = 10),legend.title = element_text())
#### plot celltypes ratio in each sample-compartment, faceted by age-sample
facet_2<-function(my.data){
my.data3<-my.data[,c("celltypes_new_redo2","age","Index","pos")]
my.data3$var<-paste(my.data3$celltypes_new_redo2,my.data3$age,my.data3$Index,my.data3$pos,sep="#")
my.data3<-as.data.frame(table(my.data3$var))
my.data3$celltypes_new_redo2<-str_split_fixed(my.data3$Var1,"#",4)[,1]
my.data3$age<-str_split_fixed(my.data3$Var1,"#",4)[,2]
my.data3$Index<-str_split_fixed(my.data3$Var1,"#",4)[,3]
my.data3$pos<-str_split_fixed(my.data3$Var1,"#",4)[,4]
my.data3<-as.data.frame(my.data3)
tg3<-my.data3 %>% dplyr::group_by(age,Index,pos)  %>% dplyr::mutate(ratio = Freq/sum(Freq))
tg3$Var2<-paste(tg3$Index,tg3$age,sep="#")
mycolors <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(my.data$celltypes_new_redo2)))
tg3$Var2<-factor(tg3$Var2,levels = c("3#young","4#young","5#young","6#young","1#old","2#old","5#old","6#old"))
tg3$pos<-factor(tg3$pos,levels = c("Cecum","C_P","C_D"))
p<-ggplot(tg3, aes(fill=celltypes_new_redo2, y=ratio, x=pos)) +
geom_bar(position="fill", stat="identity")+scale_fill_manual(values = mycolors)+facet_wrap(~Var2,nrow=2)+
theme(strip.text.x = element_text(size = 10),axis.text.x = element_text(angle=90, hjust=1))
newList <- list("figure" = p, "variable" = tg3)
return(newList)
}
xx<-facet_2(my.data.redo2)
xx$figure
yy<-xx$variable
yy$celltypes_new_redo2<-factor(yy$celltypes_new_redo2,levels = c("CD4-Early","CD4-Naive","CD4","CD4-Cytotoxic","CD4-TEM","Th17","CD8","Bcell-Early","Bcell-Naive-Isg","B-cell","Plasma","Neutrophil","Monocyte","Macrophage"))

ggplot(yy, aes(fill=celltypes_new_redo2, y=ratio, x=pos)) +
geom_bar(position="fill", stat="identity")+scale_fill_manual(values = mycolors2)+facet_wrap(~Var2,nrow=2)+
theme(strip.text.x = element_text(size = 10),axis.text.x = element_text(angle=90, hjust=1))

#### plot cell type ratio change in aging by compartment, faceted by celltypes, ggpubr did not work, using ggplot instead
yy$age<-factor(yy$age,levels = c("young","old"))
mus_2 <- summarySE(yy, measurevar="ratio",
groupvars=c("celltypes_new_redo2", "pos", "age"), na.rm = TRUE)

ggplot(mus_2, aes(x=pos, y=ratio, fill=age)) +
geom_bar(stat="identity", colour="black",position="dodge") +
facet_wrap(~celltypes_new_redo2,nrow=2) + geom_errorbar(aes(ymin=ratio-sd, ymax=ratio+sd), size=0.5, width=.25,position=position_dodge(.9)) +ggtitle("mean_sd_t.test")+theme_classic()+scale_fill_manual(values = c("gray","gray40"))+stat_compare_means(data=yy,aes(group = age), label = "p.signif",method="t.test")+theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1))


#### odds ratio
odds_ratio_redo2<-function(tg,pos){
#### odds ratio
if(pos!="ALL"){tg<-tg[tg$pos==pos,]}
tg2<-tg %>% dplyr::group_by(celltypes_new_redo2,age)  %>% dplyr::mutate(celltypes_sum = sum(Freq))
tg2<-unique(tg2[,c(3,4,9)])
tg2<-tg2 %>% dplyr::group_by(age)  %>% dplyr::mutate(celltypes_sum_ratio = celltypes_sum/sum(celltypes_sum))
require(reshape2)
tg2_wide<-dcast(tg2,celltypes_new_redo2~tg2$age,value.var = 'celltypes_sum_ratio')
tg2_wide$odds_ratio<-(tg2_wide$old/(1-tg2_wide$old))/(tg2_wide$young/(1-tg2_wide$young))
#### p value by fisher exact test
tg2<-tg2 %>% dplyr::group_by(age) %>% dplyr::mutate(age_sum = sum(celltypes_sum))
tg2_wide2<-dcast(tg2,celltypes_new_redo2~tg2$age,value.var = c("celltypes_sum"))
tg2_wide2_2<-dcast(tg2,celltypes_new_redo2~tg2$age,value.var = c("age_sum"))
tg2_wide2<-merge(tg2_wide2,tg2_wide2_2,by="celltypes_new_redo2",suffixes=c(".specific",".total"))
tg2_wide2$k<-(tg2_wide2$young.specific+tg2_wide2$old.specific)
colnames(tg2_wide2)<-c("celltypes_new_redo2","x_young.specific","q_old.specific","n_young.total","m_old.total","k")
for (i in 1:nrow(tg2_wide2)){tg2_wide2$conf_1[i]<-fisher_conf(tg2_wide2[i,"q_old.specific"],tg2_wide2[i,"m_old.total"],tg2_wide2[i,"n_young.total"],tg2_wide2[i,"x_young.specific"])[1];tg2_wide2$conf_2[i]<-fisher_conf(tg2_wide2[i,"q_old.specific"],tg2_wide2[i,"m_old.total"],tg2_wide2[i,"n_young.total"],tg2_wide2[i,"x_young.specific"])[2]}
tg2_wide2_m<-merge(tg2_wide,tg2_wide2,by = "celltypes_new_redo2",all = T)
#### 95%CI value hypergeometric test
tg2_wide2_m$hypergeometric_p_lower<-phyper(tg2_wide2$q_old.specific,tg2_wide2$m_old.total,tg2_wide2$n_young.total,tg2_wide2$k)
tg2_wide2_m$hypergeometric_p_upper<-phyper(tg2_wide2$q_old.specific,tg2_wide2$m_old.total,tg2_wide2$n_young.total,tg2_wide2$k,lower.tail = F)
tg2_wide2_m$hypergeometric_p_tails<-pmin(tg2_wide2_m$hypergeometric_p_lower,tg2_wide2_m$hypergeometric_p_upper)
require(ggplot2)
p<-ggplot(tg2_wide2_m, aes(x=celltypes_new_redo2, y=odds_ratio)) +
geom_bar(stat="identity", color="black") +
geom_errorbar(aes(ymin=conf_1, ymax=conf_2), width=.2,
position=position_dodge(.9)) + coord_flip()+scale_x_discrete(limits = rev(levels(factor(tg2_wide2_m$celltypes_new_redo2))))+theme_classic(base_size = 20) + geom_hline(yintercept=1, linetype="dashed", color = "black")+geom_text(aes(label=signif(hypergeometric_p_tails,digits = 3)),size = 3, position=position_dodge(width=1),hjust = -1,vjust=2)
newList <- list("figure" = p, "variable" = tg2_wide2_m)
return(newList)
}

odds_ratio_redo2(yy,"Cecum")[,1]
test<-odds_ratio_redo2(yy,"Cecum")
test$figure
test<-odds_ratio_redo2(yy,"C_P")
test$figure
test<-odds_ratio_redo2(yy,"C_D")
cd2<-test$variable

p1<-ggplot(cd2,aes(x=celltypes_new_redo2, y=odds_ratio)) +
geom_bar(stat="identity", color="black") +
geom_errorbar(aes(ymin=conf_1, ymax=conf_2), width=.2,
position=position_dodge(.9)) + scale_x_discrete(limits = rev(levels(factor(cd2$celltypes_new_redo2))))+theme_classic(base_size = 20) + geom_hline(yintercept=1, linetype="dashed", color = "black")+geom_text(aes(label=signif(hypergeometric_p_tails,digits = 3)),size = 3, position=position_dodge(width=1),hjust = -1,vjust=2)+coord_flip()
p2<-ggplot(cd2,aes(x=celltypes_new_redo2, y=odds_ratio)) +
geom_bar(stat="identity", color="black") +
geom_errorbar(aes(ymin=conf_1, ymax=conf_2), width=.2,
position=position_dodge(.9)) + scale_x_discrete(limits = rev(levels(factor(cd2$celltypes_new_redo2))))+theme_classic(base_size = 20) + geom_hline(yintercept=1, linetype="dashed", color = "black")+geom_text(aes(label=signif(hypergeometric_p_tails,digits = 3)),size = 3, position=position_dodge(width=1),hjust = -1,vjust=2)+coord_flip()+ylim(0,10)
cowplot::plot_grid(p1, p2)
test<-odds_ratio_redo2(yy,"ALL")
test$figure

rm(batch.anchors)
rm(batch.anchors.2)
rm(immune.combined.2)
#### pheatmap with top20 markers
immune_markers<-read.csv("./data_in/immune_markers.csv",head=T,stringsAsFactors = F)
immune_markers$scm<-paste(immune_markers$source,immune_markers$cell_type,immune_markers$markers,sep="#")
#### recall markers with latest defined celltypes
Idents(myData_LI)<-myData_LI$celltypes_new_redo2
myData_LI_markers_redo2 <- FindAllMarkers(object = myData_LI, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(myData_LI_markers_redo2,"./data_out/myData_LI_markers_redo2.csv")
myData_LI_markers_redo2 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC) -> top20_target.redo2
top20_immune<-merge(top20_target.redo2,immune_markers,by.x='gene',by.y='markers',all =T)
top20_immune2<-top20_immune[(!is.na(top20_immune$cluster))&(!is.na(top20_immune$scm)),]
write.csv(top20_immune2,"./data_out/top20_immune2.csv")
#### seperate pdf files
for(i in unique(top20_immune2$gene)){pdf(paste(i,".pdf",sep=""),width=3,height = 3.5);plots<-VlnPlot(myData_LI, features = i,  group.by = "celltypes_new_redo2", pt.size = 0,cols=rep("grey",length(unique(myData_LI$celltypes_new_redo2))))+NoLegend()+ggtitle(top20_immune2[top20_immune2$gene==i,"scm"])+theme( plot.title = element_text(size = 8, face = "bold"),axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 8))+xlab("")+ylab("");print(plots);dev.off()}
#### combined pdf files
pdf("all_markers_vln")
for(i in unique(top20_immune2$gene)){plots<-VlnPlot(myData_LI, features = i,  group.by = "celltypes_new_redo2", pt.size = 0,cols=rep("grey",length(unique(myData_LI$celltypes_new_redo2))))+NoLegend()+ggtitle(top20_immune2[top20_immune2$gene==i,"scm"])+theme( plot.title = element_text(size = 8, face = "bold"),axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 8))+xlab("")+ylab("");print(plots)}
dev.off()

#### heatmap with top20 markers
pp<-DoHeatmap(myData_LI, features = top20_target.redo2$gene,size = 2)
pp+ NoLegend()+ theme(axis.text.y = element_text(size = 8),
legend.title = element_text(size = 4),
legend.text = element_text(size = 8)
) + scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))


#### calling aging DEG
myData_LI$celltype_age<-paste(myData_LI$celltypes_new_redo2,myData_LI$age,sep=".")
celltypes<-unique(myData_LI$celltypes_new_redo2)
for(i in celltypes){assign(paste("age.DEG.wilcox",i,sep = "_") , FindMarkers(myData_LI, ident.1 = paste(i,"old",sep="."), ident.2 = paste(i,"young",sep="."), verbose = FALSE,group.by='celltype_age'))}
for(i in celltypes){write.csv(get(paste("age.DEG.wilcox",i,sep = "_")),paste(i,"age.DEG.wilcox.csv",sep = "_"))}
rm(myData_LI_Combined)

save.image(file = "./data_in/Atlas.RData")
