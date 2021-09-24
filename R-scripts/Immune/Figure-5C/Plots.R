#################################################################
############### R and Packages version ##########################
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
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ggpubr_0.2.3    magrittr_1.5    Seurat_3.0.2    forcats_0.4.0   pheatmap_1.0.12 stringr_1.4.0   dplyr_0.8.0.1  
[8] ggplot2_3.1.1  

loaded via a namespace (and not attached):
 [1] nlme_3.1-139        tsne_0.1-3          bitops_1.0-6        RColorBrewer_1.1-2  httr_1.4.0         
 [6] sctransform_0.2.0   tools_3.5.1         R6_2.4.0            irlba_2.3.3         KernSmooth_2.23-15 
[11] lazyeval_0.2.2      colorspace_1.4-1    npsurv_0.4-0        withr_2.1.2         tidyselect_0.2.5   
[16] gridExtra_2.3       compiler_3.5.1      plotly_4.9.0        caTools_1.17.1.2    scales_1.0.0       
[21] lmtest_0.9-37       ggridges_0.5.1      pbapply_1.4-0       digest_0.6.18       R.utils_2.9.0      
[26] pkgconfig_2.0.2     htmltools_0.4.0     bibtex_0.4.2        htmlwidgets_1.5.1   rlang_0.4.5        
[31] rstudioapi_0.10     zoo_1.8-6           jsonlite_1.6        ica_1.0-2           gtools_3.8.1       
[36] R.oo_1.22.0         Matrix_1.2-17       Rcpp_1.0.1          munsell_0.5.0       ape_5.3            
[41] reticulate_1.12     R.methodsS3_1.7.1   stringi_1.4.3       yaml_2.2.0          gbRd_0.4-11        
[46] MASS_7.3-51.4       gplots_3.0.1.1      Rtsne_0.15          plyr_1.8.4          grid_3.5.1         
[51] parallel_3.5.1      gdata_2.18.0        listenv_0.7.0       ggrepel_0.8.1       crayon_1.3.4       
[56] lattice_0.20-38     cowplot_0.9.4       splines_3.5.1       SDMTools_1.1-221.1  pillar_1.3.1       
[61] igraph_1.2.4.1      ggsignif_0.6.0      future.apply_1.3.0  reshape2_1.4.3      codetools_0.2-16   
[66] glue_1.3.1          lsei_1.2-0          metap_1.1           data.table_1.12.2   png_0.1-7          
[71] Rdpack_0.11-0       gtable_0.3.0        RANN_2.6.1          purrr_0.3.2         tidyr_0.8.3        
[76] future_1.13.0       assertthat_0.2.1    rsvd_1.0.1          survival_2.44-1.1   viridisLite_0.3.0  
[81] tibble_2.1.1        cluster_2.0.8       globals_0.12.4      fitdistrplus_1.0-14 ROCR_1.0-7   

#################################################################
############### Load Data and Packages ##########################
#################################################################
setwd("./")
load("./data_in/Atlas.RData")
library(dplyr)
library(stringr)
library(pheatmap)
library(forcats)
library(Seurat)
library(ggpubr)
library(Seurat)
library(stringr)

#################################################################
############### REORDER AND RECOLOUR THE CELLTYPES ##############
#################################################################


logic_order<-c("Bcell-Early","B-cell","Bcell-Naive-Isg","Plasma","CD4-Early","CD4-Naive","CD4","CD4-Cytotoxic","CD4-TEM","CD8","Th17","Monocyte","Macrophage","Neutrophil")
population_order<-c("B-cell","CD4","CD4-TEM","CD8","Plasma","CD4-Naive","CD4-Early","Neutrophil","Th17","CD4-Cytotoxic","Bcell-Naive-Isg","Monocyte","Bcell-Early","Macrophage")
> mycolors3
  "#84BF96" "#6DBD57" "#7F9D55" "#3385BB" "#DE9E83" "#FBB268" "#FE8D19"
  "#F57C7C" "#E42622" "#891c1e" "#F3E587" "#977899" "#9D7BBA" "#B15928"



#### MARKERS HEATMAP AND VLNPLOT
#### top 20 markers dispalayed their expression by raw counts

Idents(myData_LI)<-myData_LI$celltypes_new_redo2
Idents(myData_LI)<-factor(Idents(myData_LI),levels = c("B-cell","CD4","CD4-TEM","CD8","Plasma","CD4-Naive","CD4-Early","Neutrophil","Th17","CD4-Cytotoxic","Bcell-Naive-Isg","Monocyte","Bcell-Early","Macrophage"))
top20_target.redo2$cluster<-factor(top20_target.redo2$cluster,levels = c("B-cell","CD4","CD4-TEM","CD8","Plasma","CD4-Naive","CD4-Early","Neutrophil","Th17","CD4-Cytotoxic","Bcell-Naive-Isg","Monocyte","Bcell-Early","Macrophage"))
top20_target.redo2<-top20_target.redo2[order(top20_target.redo2$cluster),]

#### check known markers
immune_markers_3<-read.csv("./data_in/immune_markers_3.csv",head=T,stringsAsFactors = F)
Idents(myData_LI)<-factor(Idents(myData_LI),levels = c("Bcell-Early","B-cell","Bcell-Naive-Isg","Plasma","CD4-Early","CD4-Naive","CD4","CD4-Cytotoxic","CD4-TEM","CD8","Th17","Monocyte","Macrophage","Neutrophil"))
#### change color scheme to mycolors3, how mycolors3 was created 
mycolors3<-mycolors2
mycolors3<-mycolors3[c(3,4,5,2,10,8,9,6,7,1,13,11,12,14)]
mycolors3[10]<-"#891c1e"
mycolors3<-mycolors3[c(1:11,13,12,14)]
VlnPlot(myData_LI, features = "Ifng",pt.size = 0.2)+scale_fill_manual(values = mycolors3)
#### all in one file, good idea for presenataion but not for publication. this will generate repetive for those gene that are not in the object, maybe print(plots) ought to be in "try()" 
pdf("reorder_all_known_markers_vln.pdf")
for(i in unique(immune_markers_3$markers)){try(plots<-VlnPlot(myData_LI, features = i, pt.size = 0.2)+scale_fill_manual(values = mycolors3))
print(plots)}
dev.off()
#### /Users/jlu/Desktop/Pro_Mine/Paper_Atlas/Immune/figure/Reorder/immune_markers_check_vln/reorder_known_markers_vln_gene.pdf
#### removed repetitive ones manually 
## Figure-S5B
for(i in unique(immune_markers_3$markers)){try(plots<-VlnPlot(myData_LI, features = i, pt.size = 0.2)+scale_fill_manual(values = mycolors3))
pdf(paste("reorder_known_markers_vln_",i,".pdf",sep=""))
print(plots)
dev.off()
}
##

## Figure-5A, UMAP PLOT, ONLY YOUNG
logic_order<-c("Bcell-Early","B-cell","Bcell-Naive-Isg","Plasma","CD4-Early","CD4-Naive","CD4","CD4-Cytotoxic","CD4-TEM","CD8","Th17","Monocyte","Macrophage","Neutrophil")
population_order<-c("B-cell","CD4","CD4-TEM","CD8","Plasma","CD4-Naive","CD4-Early","Neutrophil","Th17","CD4-Cytotoxic","Bcell-Naive-Isg","Monocyte","Bcell-Early","Macrophage")
myData_LI$pos<-factor(myData_LI$pos,levels = c("Cecum","C_P","C_D"))
myData_LI$celltypes_new_redo2<-factor(myData_LI$celltypes_new_redo2,levels = logic_order)
myData_LI$pos_age<-paste(myData_LI$pos,myData_LI$age,sep="#")
myData_LI$pos_age<-factor(myData_LI$pos_age,levels = c("Cecum#young","C_P#young","C_D#young","Cecum#old","C_P#old","C_D#old"))
###  reorder_Figure5A_1.pdf, with cell type labelled
DimPlot(object = subset(myData_LI, subset = age == "young"), reduction = 'umap',split.by = "pos_age",group.by="celltypes_new_redo2",label = F,cols =mycolors2,ncol=3)+coord_fixed(ratio = 1)+scale_color_manual(values=mycolors3)
### reorder_Figure5A_2.pdf, without cell type labelled
DimPlot(object = subset(myData_LI, subset = age == "young"), reduction = 'umap',split.by = "pos_age",group.by="celltypes_new_redo2",label = T,cols =mycolors2,ncol=3)+coord_fixed(ratio = 1)+scale_color_manual(values=mycolors3)
##

#### ODDS RATIO
yy$celltypes_new_redo2<-factor(yy$celltypes_new_redo2,levels = logic_order)
## Figure-6A: reorder_oddsRatio_Cecum.pdf
test<-odds_ratio_redo2(yy,"Cecum")
test$figure
## Figure-6A: reorder_oddsRatio_CP.pdf
test<-odds_ratio_redo2(yy,"C_P")
test$figure
## Figure-6A: reorder_oddsRatio_CD.pdf 
test<-odds_ratio_redo2(yy,"C_D")
test$figure
# reorder_oddsRatio_CD2.pdf
# notice: with or without Bcell-Naive-Isg, the p value will be different, p value ought be referred to the one with Bcell-Naive-Isg.
# test<-odds_ratio_redo2(yy[yy$celltypes_new_redo2!="Bcell-Naive-Isg",],"C_D") # did not keep the empty space
# test$figure

yyx<-yy
yyx[yyx$celltypes_new_redo2=="Bcell-Naive-Isg","Freq"]<-0
yyx[yyx$celltypes_new_redo2=="Bcell-Naive-Isg","ratio"]<-0
test<-odds_ratio_redo2(yyx,"C_D")
test$figure

#### UMAP IN YOUNG AND OLD 
myData_LI$pos<-factor(myData_LI$pos,levels = c("Cecum","C_P","C_D"))
myData_LI$celltypes_new_redo2<-factor(myData_LI$celltypes_new_redo2,levels = logic_order)
myData_LI$pos_age<-paste(myData_LI$pos,myData_LI$age,sep="#")
myData_LI$pos_age<-factor(myData_LI$pos_age,levels = c("Cecum#young","C_P#young","C_D#young","Cecum#old","C_P#old","C_D#old"))
## Figure-6B: reorder_figures6_umap_age.pdf
DimPlot(object = myData_LI, reduction = 'umap',split.by = "age",group.by="celltypes_new_redo2",label = T,cols =mycolors3)+coord_fixed(ratio = 1)
## Figure-6B: reorder_figures6_umap_age2.pdf
DimPlot(object = myData_LI, reduction = 'umap',split.by = "age",group.by="celltypes_new_redo2",label = F,cols =mycolors3)+coord_fixed(ratio = 1)

#### FIGURES 6
Num_DEGs<-read.csv("./data_in/Num_DEGs_Sig.csv",head=T)
Num_DEGs$Pos<-factor(Num_DEGs$Pos,levels = c("Cecum","C_P","C_D"))
Num_DEGs$CellTypes<-factor(Num_DEGs$CellTypes,levels = logic_order)
## Figure-6C: reorder_figures6_DEGs_number.pdf
ggplot(data=Num_DEGs, aes(x=CellTypes, y=Num_DEGs_Sig)) + geom_bar(stat="identity")+facet_grid(Pos~.)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),text = element_text(size=8))

#### FIGURES S6
## Figure-S6B: reorder_figuresS6_barfill.pdf
ggplot(yy, aes(fill=celltypes_new_redo2, y=ratio, x=pos)) +
geom_bar(position="fill", stat="identity")+scale_fill_manual(values = mycolors3)+facet_wrap(~Var2,nrow=2)+
theme(strip.text.x = element_text(size = 10),axis.text.x = element_text(angle=90, hjust=1))
myData_LI$age_pos<-paste(myData_LI$age,myData_LI$pos,sep="#")
myData_LI$age_pos<-factor(myData_LI$age_pos,levels = c("young#Cecum","young#C_P","young#C_D","old#Cecum","old#C_P","old#C_D"))
#### scalling process could be done by seleted genes
myData_LI <- ScaleData(object = myData_LI, features = top20_target.redo2$gene, block.size = 2000)
#### SingleImmune_BcellMarkers_Compartment.tiff
DoHeatmap(myData_LI, features = as.character(top20_target.redo2[top20_target.redo2$cluster=="B-cell",]$gene),size = 2,group.by = "age_pos",slot = "data")

#### FIGURE3 
#### for customed plots from pipeline plot

cellphonedb_all<-read.table("./data_in/all_network_count.txt",sep="\t",head=T,stringsAsFactors = F)
cellphonedb_all$SOURCE_celltype<-str_split_fixed(cellphonedb_all$SOURCE,"\\.",3)[,3]
cellphonedb_all$TARGET_celltype<-str_split_fixed(cellphonedb_all$TARGET,"\\.",3)[,3]
cellphonedb_all$SOURCE_age<-str_split_fixed(cellphonedb_all$SOURCE,"\\.",3)[,2]
cellphonedb_all$TARGET_age<-str_split_fixed(cellphonedb_all$TARGET,"\\.",3)[,2]
cellphonedb_all$SOURCE_pos<-str_split_fixed(cellphonedb_all$SOURCE,"\\.",3)[,1]
cellphonedb_all$TARGET_pos<-str_split_fixed(cellphonedb_all$TARGET,"\\.",3)[,1]
B_cell_lineage<-c("Bcell-Early","B-cell","Bcell-Naive-Isg","Plasma")
T_cell_lineage<-c("CD4-Early","CD4-Naive","CD4","CD4-Cytotoxic","CD4-TEM","CD8","Th17")
Myeloid_lineage<-c("Monocyte","Macrophage","Neutrophil")
Epithelial<-c("ProxCol-enriched-Colonocyte","Tuft_Cell","Stem_cell","EnteroEndocrine_Cell","TA-cell","Goblet_Precursor","Goblet","Colonocyte_Precursor","Cecum-enriched-Colonocyte","DistCol-enriched-Colonocyte")
cellphonedb_all$SOURCE_celltype_super<-"SOURCE_celltype_super"
cellphonedb_all$TARGET_celltype_super<-"TARGET_celltype_super"
cellphonedb_all[cellphonedb_all$SOURCE_celltype %in% B_cell_lineage,"SOURCE_celltype_super"]<-"B_cell_lineage"
cellphonedb_all[cellphonedb_all$SOURCE_celltype %in% T_cell_lineage,"SOURCE_celltype_super"]<-"T_cell_lineage"
cellphonedb_all[cellphonedb_all$SOURCE_celltype %in% Myeloid_lineage,"SOURCE_celltype_super"]<-"Myeloid_lineage"
cellphonedb_all[cellphonedb_all$SOURCE_celltype %in% Epithelial,"SOURCE_celltype_super"]<-"Epithelial"
cellphonedb_all[cellphonedb_all$TARGET_celltype %in% B_cell_lineage,"TARGET_celltype_super"]<-"B_cell_lineage"
cellphonedb_all[cellphonedb_all$TARGET_celltype %in% T_cell_lineage,"TARGET_celltype_super"]<-"T_cell_lineage"
cellphonedb_all[cellphonedb_all$TARGET_celltype %in% Myeloid_lineage,"TARGET_celltype_super"]<-"Myeloid_lineage"
cellphonedb_all[cellphonedb_all$TARGET_celltype %in% Epithelial,"TARGET_celltype_super"]<-"Epithelial"
cellphonedb_all$pos.age.SOURCE_celltype_super<-paste(cellphonedb_all$SOURCE_pos,cellphonedb_all$SOURCE_age,cellphonedb_all$SOURCE_celltype_super,sep=".")
cellphonedb_all$pos.age.TARGET_celltype_super<-paste(cellphonedb_all$TARGET_pos,cellphonedb_all$TARGET_age,cellphonedb_all$TARGET_celltype_super,sep=".")
cellphonedb_all$pos.age.INTERACT_celltype_super<-paste(cellphonedb_all$pos.age.SOURCE_celltype_super,cellphonedb_all$pos.age.TARGET_celltype_super,sep="#")
cellphonedb_all_super<-aggregate(cellphonedb_all$count, by=list(Category=cellphonedb_all$pos.age.INTERACT_celltype_super), FUN=sum)
cellphonedb_all_super$SOURCE<-str_split_fixed(cellphonedb_all_super$Category,"#",2)[,1]
cellphonedb_all_super$TARGET<-str_split_fixed(cellphonedb_all_super$Category,"#",2)[,2]
cellphonedb_all_super$SOURCE_pos<-str_split_fixed(cellphonedb_all_super$SOURCE,"\\.",3)[,1]
cellphonedb_all_super$SOURCE_age<-str_split_fixed(cellphonedb_all_super$SOURCE,"\\.",3)[,2]
cellphonedb_all_super$SOURCE_celltype<-str_split_fixed(cellphonedb_all_super$SOURCE,"\\.",3)[,3]
cellphonedb_all_super$TARGET_celltype<-str_split_fixed(cellphonedb_all_super$TARGET,"\\.",3)[,3]
cellphonedb_all_super$TARGET_age<-str_split_fixed(cellphonedb_all_super$TARGET,"\\.",3)[,2]
cellphonedb_all_super$TARGET_pos<-str_split_fixed(cellphonedb_all_super$TARGET,"\\.",3)[,1]
cellphonedb_all_super[cellphonedb_all_super$SOURCE_pos=="ProxCol","SOURCE_pos"]<-"C_P"
cellphonedb_all_super[cellphonedb_all_super$SOURCE_pos=="DistCol","SOURCE_pos"]<-"C_D"
cellphonedb_all_super[cellphonedb_all_super$TARGET_pos=="DistCol","TARGET_pos"]<-"C_D"
cellphonedb_all_super[cellphonedb_all_super$TARGET_pos=="ProxCol","TARGET_pos"]<-"C_P"
cellphonedb_all_super[cellphonedb_all_super$TARGET_age=="Young","TARGET_age"]<-"young"
cellphonedb_all_super[cellphonedb_all_super$TARGET_age=="Old","TARGET_age"]<-"old"
cellphonedb_all_super[cellphonedb_all_super$SOURCE_age=="Old","SOURCE_age"]<-"old"
cellphonedb_all_super[cellphonedb_all_super$SOURCE_age=="Young","SOURCE_age"]<-"young"
cellphonedb_all_super$SOURCEvsTARGET<-paste(cellphonedb_all_super$SOURCE_celltype,cellphonedb_all_super$TARGET_celltype,sep="#")
cellphonedb_all_super$SOURCE_pos<-factor(cellphonedb_all_super$SOURCE_pos,levels = c("Cecum","C_P","C_D"))
cellphonedb_all_super$SOURCE_age<-factor(cellphonedb_all_super$SOURCE_age,levels = c("young","old"))
ggplot(data=cellphonedb_all_super, aes(x=SOURCEvsTARGET, y=x, fill=SOURCE_age)) +
geom_bar(stat="identity", position=position_dodge())+
geom_text(aes(label=x), vjust=-0.3, color="black",
position = position_dodge(0.9), size=2)+
scale_fill_manual(values=c("gray","gray0"))+facet_wrap(~cellphonedb_all_super$SOURCE_pos,ncol=1)+theme(text = element_text(size=8))+ylab("sum of significant interactions")+xlab("SOURCE#TARGET")+theme(axis.text.x=element_text(angle=90,vjust=0.5))
cellphonedb_all$SOURCEvsTARGET<-paste(cellphonedb_all$SOURCE_celltype,cellphonedb_all$TARGET_celltype,sep="#")
cellphonedb_all[cellphonedb_all$SOURCE_pos=="ProxCol","SOURCE_pos"]<-"C_P"
cellphonedb_all[cellphonedb_all$SOURCE_pos=="DistCol","SOURCE_pos"]<-"C_D"
cellphonedb_all[cellphonedb_all$TARGET_pos=="DistCol","TARGET_pos"]<-"C_D"
cellphonedb_all[cellphonedb_all$TARGET_pos=="ProxCol","TARGET_pos"]<-"C_P"
cellphonedb_all[cellphonedb_all$TARGET_age=="Young","TARGET_age"]<-"young"
cellphonedb_all[cellphonedb_all$TARGET_age=="Old","TARGET_age"]<-"old"
cellphonedb_all[cellphonedb_all$SOURCE_age=="Old","SOURCE_age"]<-"old"
cellphonedb_all[cellphonedb_all$SOURCE_age=="Young","SOURCE_age"]<-"young"
cellphonedb_all$SOURCE_pos<-factor(cellphonedb_all$SOURCE_pos,levels = c("Cecum","C_P","C_D"))
cellphonedb_all$SOURCE_age<-factor(cellphonedb_all$SOURCE_age,levels = c("young","old"))
colnames(cellphonedb_all)<- c("TARGET","SOURCE","count","TARGET_celltype","SOURCE_celltype","TARGET_age","SOURCE_age","TARGET_pos","SOURCE_pos","TARGET_celltype_super","SOURCE_celltype_super","pos.age.TARGET_celltype_super","pos.age.SOURCE_celltype_super","pos.age.INTERACT_celltype_super","SOURCEvsTARGET")
cellphonedb_all$SOURCE_celltype<-gsub("-",".",cellphonedb_all$SOURCE_celltype)
cellphonedb_all$TARGET_celltype<-gsub("-",".",cellphonedb_all$TARGET_celltype)
cellphonedb_all$SOURCEvsTARGET<-paste(cellphonedb_all$SOURCE_celltype,cellphonedb_all$TARGET_celltype,sep="#")
xxx<-cellphonedb_all
logic_order
cellphonedb_all$TARGET_celltype<-factor(cellphonedb_all$TARGET_celltype,levels = c("Stem_cell","TA.cell","Colonocyte_Precursor","Cecum.enriched.Colonocyte","ProxCol.enriched.Colonocyte","DistCol.enriched.Colonocyte","EnteroEndocrine_Cell","Goblet_Precursor","Goblet","Tuft_Cell","Bcell.Early","B.cell","Bcell.Naive.Isg","Plasma","CD4.Early","CD4.Naive","CD4","CD4.Cytotoxic","CD4.TEM","CD8","Th17","Monocyte","Macrophage","Neutrophil"))
cellphonedb_all$SOURCE_celltype<-factor(cellphonedb_all$SOURCE_celltype,levels = c("Stem_cell","TA.cell","Colonocyte_Precursor","Cecum.enriched.Colonocyte","ProxCol.enriched.Colonocyte","DistCol.enriched.Colonocyte","EnteroEndocrine_Cell","Goblet_Precursor","Goblet","Tuft_Cell","Bcell.Early","B.cell","Bcell.Naive.Isg","Plasma","CD4.Early","CD4.Naive","CD4","CD4.Cytotoxic","CD4.TEM","CD8","Th17","Monocyte","Macrophage","Neutrophil"))
cellphonedb_all$SOURCE_age<-factor(cellphonedb_all$SOURCE_age,levels = c("young","old"))
cellphonedb_all$SOURCE_pos<-factor(cellphonedb_all$SOURCE_pos,levels = c("Cecum","C_P","C_D"))
## Figure-7A: reorder_FN_msg3_Bar_Myeloid_Epithelial
tmp<-cellphonedb_all[cellphonedb_all$SOURCE_celltype_super=="Myeloid_lineage" & cellphonedb_all$TARGET_celltype_super=="Epithelial",]
ggplot(data=tmp, aes(x=TARGET_celltype, y=count, fill=SOURCE_age)) +
geom_bar(stat="identity", position=position_dodge())+
geom_text(aes(label=count), vjust=-0.3, color="black",
position = position_dodge(0.9), size=2)+
scale_fill_manual(values=c("gray","gray0"))+facet_grid(tmp$SOURCE_celltype~tmp$SOURCE_pos)+theme(text = element_text(size=8))+ylab("sum of significant interactions")+xlab("SOURCE#TARGET")+theme(axis.text.x=element_text(angle=90,vjust=0.5))
## Figure-7A: reorder_FN_msg3_Bar_Tcell_Epithelial.pdf
tmp<-cellphonedb_all[cellphonedb_all$SOURCE_celltype_super=="T_cell_lineage" & cellphonedb_all$TARGET_celltype_super=="Epithelial",]
ggplot(data=tmp, aes(x=TARGET_celltype, y=count, fill=SOURCE_age)) +
geom_bar(stat="identity", position=position_dodge())+
geom_text(aes(label=count), vjust=-0.3, color="black",
position = position_dodge(0.9), size=2)+
scale_fill_manual(values=c("gray","gray0"))+facet_grid(tmp$SOURCE_celltype~tmp$SOURCE_pos)+theme(text = element_text(size=8))+ylab("sum of significant interactions")+xlab("SOURCE#TARGET")+theme(axis.text.x=element_text(angle=90,vjust=0.5))
## Figure-7A: reorder_FN_msg3_Bar_Bcell_Epithelial.pdf
tmp<-cellphonedb_all[cellphonedb_all$SOURCE_celltype_super=="B_cell_lineage" & cellphonedb_all$TARGET_celltype_super=="Epithelial",]
ggplot(data=tmp, aes(x=TARGET_celltype, y=count, fill=SOURCE_age)) +
geom_bar(stat="identity", position=position_dodge())+
geom_text(aes(label=count), vjust=-0.3, color="black",
position = position_dodge(0.9), size=2)+
scale_fill_manual(values=c("gray","gray0"))+facet_grid(tmp$SOURCE_celltype~tmp$SOURCE_pos)+theme(text = element_text(size=8))+ylab("sum of significant interactions")+xlab("SOURCE#TARGET")+theme(axis.text.x=element_text(angle=90,vjust=0.5))
## Figure-7A: reorder_FN_msg3_Bar_Epithelial_Epithelial.pdf
tmp<-cellphonedb_all[cellphonedb_all$SOURCE_celltype_super=="Epithelial" & cellphonedb_all$TARGET_celltype_super=="Epithelial",]
ggplot(data=tmp, aes(x=TARGET_celltype, y=count, fill=SOURCE_age)) +
geom_bar(stat="identity", position=position_dodge())+
geom_text(aes(label=count), vjust=-0.3, color="black",
position = position_dodge(0.9), size=2)+
scale_fill_manual(values=c("gray","gray0"))+facet_grid(tmp$SOURCE_celltype~tmp$SOURCE_pos)+theme(text = element_text(size=8))+ylab("sum of significant interactions")+xlab("SOURCE#TARGET")+theme(axis.text.x=element_text(angle=90,vjust=0.5))

#### REMOVE IMMUNE RECEPTORS IN THE FIRST METHOD OF CELLPHONEDB
# FN msg3:
# reorder_Method1_Old_CD.pdf     reorder_Method1_Young_CD.pdf   reorder_Method1_Old_CP.pdf     reorder_Method1_Young_CP.pdf   reorder_Method1_Old_Cecum.pdf   reorder_Method1_Young_Cecum.pdf
# /Users/jlu/Desktop/Pro_Mine/Paper_Atlas/Immune/figure/Reorder/

heatmap_count<-function(count_matrix){
row.names(count_matrix)<-colnames(count_matrix)
df<-data.frame(coln=colnames(count_matrix),leve=str_split_fixed(colnames(count_matrix),"\\.",3)[,3])
df$leve<-as.character(df$leve)
df$coln<-as.character(df$coln)
orders_ie<-c("Stem_cell","TA.cell","Colonocyte_Precursor","Cecum.enriched.Colonocyte","ProxCol.enriched.Colonocyte","DistCol.enriched.Colonocyte","EnteroEndocrine_Cell","Goblet_Precursor","Goblet","Tuft_Cell","Bcell.Early","B.cell","Bcell.Naive.Isg","Plasma","CD4.Early","CD4.Naive","CD4","CD4.Cytotoxic","CD4.TEM","CD8","Th17","Monocyte","Macrophage","Neutrophil")
dfo<-df[match(orders_ie,df$leve),]
dfo<-na.omit(dfo)
count_matrix<-count_matrix[dfo$coln,dfo$coln]
colnames(count_matrix)<-str_split_fixed(colnames(count_matrix),"\\.",3)[,3]
row.names(count_matrix)<-str_split_fixed(row.names(count_matrix),"\\.",3)[,3]
breaksList<-seq(0, 150, by = 1)
col1 = "dodgerblue4"
col2 = 'peachpuff'
col3 = 'deeppink4'
col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( length(breaksList) )
pheatmap(count_matrix[,!colnames(count_matrix) %in% immune_receptors], show_rownames = T, show_colnames = T, scale="none", cluster_cols = F,
border_color='white', cluster_rows = F, fontsize_row = 10, fontsize_col = 10,
main = '', treeheight_row = 0, family = 'Arial',color = col.heatmap,breaks = breaksList, treeheight_col = 0)
}
## Figure-S7B: reorder_Method1_Old_CD.pdf     reorder_Method1_Young_CD.pdf   reorder_Method1_Old_CP.pdf     reorder_Method1_Young_CP.pdf   reorder_Method1_Old_Cecum.pdf   reorder_Method1_Young_Cecum.pdf
immune_receptors<-c("Bcell.Early","B.cell","Bcell.Naive.Isg","Plasma","CD4.Early","CD4.Naive","CD4","CD4.Cytotoxic","CD4.TEM","CD8","Th17","Monocyte","Macrophage","Neutrophil")
count_matrix<-read.table("./data_in/young_Cecum_interaction_count.txt",header = T,row.names=1)
heatmap_count(count_matrix)
count_matrix<-read.table("./data_in/young_C_P_interaction_count.txt",header = T,row.names=1)
heatmap_count(count_matrix)
count_matrix<-read.table("./data_in/young_C_D_interaction_count.txt",header = T,row.names=1)
heatmap_count(count_matrix)
count_matrix<-read.table("./data_in/old_Cecum_interaction_count.txt",header = T,row.names=1)
heatmap_count(count_matrix)
count_matrix<-read.table("./data_in/old_C_P_interaction_count.txt",header = T,row.names=1)
heatmap_count(count_matrix)
count_matrix<-read.table("./data_in/old_C_D_interaction_count.txt",header = T,row.names=1)
heatmap_count(count_matrix)

#### 3 Heatmaps of Cecum, ProxCol, DistCol with expression delta (old/young) (first method) , remove Immunce cells receptors.
# reorder_Method1_YO_CD.pdf     reorder_Method1_YO_CP.pdf      reorder_Method1_YO_Cecum.pdf
# /Users/jlu/Desktop/Pro_Mine/Paper_Atlas/Immune/figure/Reorder/
heatmap_count_2<-function(count_matrix){
df<-data.frame(coln=colnames(count_matrix),leve=colnames(count_matrix))
df$leve<-as.character(df$leve)
df$coln<-as.character(df$coln)
orders_ie<-c("Stem_cell","TA.cell","Colonocyte_Precursor","Cecum.enriched.Colonocyte","ProxCol.enriched.Colonocyte","DistCol.enriched.Colonocyte","EnteroEndocrine_Cell","Goblet_Precursor","Goblet","Tuft_Cell","Bcell.Early","B.cell","Bcell.Naive.Isg","Plasma","CD4.Early","CD4.Naive","CD4","CD4.Cytotoxic","CD4.TEM","CD8","Th17","Monocyte","Macrophage","Neutrophil")
dfo<-df[match(orders_ie,df$leve),]
dfo<-na.omit(dfo)
count_matrix<-count_matrix[dfo$coln,dfo$coln]
breaksList<-seq(-7, 7, by = 0.1)
col1 = "dodgerblue4"
col2 = 'peachpuff'
col3 = 'deeppink4'
col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( length(breaksList) )
pheatmap(count_matrix[,!colnames(count_matrix) %in% immune_receptors], show_rownames = T, show_colnames = T, scale="none", cluster_cols = F,
border_color='white', cluster_rows = F, fontsize_row = 10, fontsize_col = 10,
main = '', treeheight_row = 0, family = 'Arial',color = col.heatmap,breaks = breaksList, treeheight_col = 0)
}

## Figure-7B: reorder_Method1_YO_CP.pdf
cm_cp_old<-read.table("./data_in/old_C_P_interaction_count.txt",header = T,row.names=1)
cm_cp_young<-read.table("./data_in/young_C_P_interaction_count.txt",header = T,row.names=1)
colnames(cm_cp_young)<-str_split_fixed(colnames(cm_cp_young),"\\.",3)[,3]
colnames(cm_cp_old)<-str_split_fixed(colnames(cm_cp_old),"\\.",3)[,3]
row.names(cm_cp_young)<-colnames(cm_cp_young)
row.names(cm_cp_old)<-colnames(cm_cp_old)
cm_cp_old<-cm_cp_old[row.names(cm_cp_young),colnames(cm_cp_young)]
cm_cp_log2fc1<-log2((cm_cp_old+1)/(cm_cp_young+1))
heatmap_count_2(cm_cp_log2fc1)
## Figure-7B: reorder_Method1_YO_CD.pdf
cm_cd_old<-read.table("./data_in/old_C_D_interaction_count.txt",header = T,row.names=1)
cm_cd_young<-read.table("./data_in/young_C_D_interaction_count.txt",header = T,row.names=1)
colnames(cm_cd_young)<-str_split_fixed(colnames(cm_cd_young),"\\.",3)[,3]
colnames(cm_cd_old)<-str_split_fixed(colnames(cm_cd_old),"\\.",3)[,3]
row.names(cm_cd_young)<-colnames(cm_cd_young)
row.names(cm_cd_old)<-colnames(cm_cd_old)
cm_cd_old<-cm_cd_old[row.names(cm_cd_young),colnames(cm_cd_young)]
cm_cd_log2fc1<-log2((cm_cd_old+1)/(cm_cd_young+1))
heatmap_count_2(cm_cd_log2fc1)
## Figure-7B: reorder_Method1_YO_Cecum.pdf
cm_cecum_old<-read.table("./data_in/old_Cecum_interaction_count.txt",header = T,row.names=1)
cm_cecum_young<-read.table("./data_in/young_Cecum_interaction_count.txt",header = T,row.names=1)
colnames(cm_cecum_young)<-str_split_fixed(colnames(cm_cecum_young),"\\.",3)[,3]
colnames(cm_cecum_old)<-str_split_fixed(colnames(cm_cecum_old),"\\.",3)[,3]
row.names(cm_cecum_young)<-colnames(cm_cecum_young)
row.names(cm_cecum_old)<-colnames(cm_cecum_old)
cm_cecum_old<-cm_cecum_old[row.names(cm_cecum_young),colnames(cm_cecum_young)]
cm_cecum_log2fc1<-log2((cm_cecum_old+1)/(cm_cecum_young+1))
heatmap_count_2(cm_cecum_log2fc1)

#### CYTOKINE DOT PLOT REORDER
test4s<-read.csv("./data_in/test4s_regardlessCompartments.csv",head=T)
## Figure-S7E: reorder_FN_test_DotPlot_regardlessCompartments3.pdf
ggplot(test4s[test4s$cytokines %in% c("Ifng","Il1b"),],aes(x= cytokines,y= age, fill = avgExp_Age))+
geom_point(aes(size=ratio_Age), shape = 21)+theme(text = element_text(size=10),axis.text.x = element_text(angle=90, hjust=1)) +scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+ theme(legend.position="bottom")
## reorder_FN_test_DotPlot_regardlessCompartments.pdf
ggplot(test4s,aes(x= cytokines,y= age, fill = avgExp_Age))+
geom_point(aes(size=ratio_Age), shape = 21)+theme(text = element_text(size=10),axis.text.x = element_text(angle=90, hjust=1)) +scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+ theme(legend.position="bottom")
## Figure-7E: reorder_FN_test_DotPlot_regardlessCompartments2.pdf
test2_regardlessCompartment<-read.csv("./data_out/test2_regardlessCompartments.csv",head=T)
test2_regardlessCompartment$celltypes_age<-factor(test2_regardlessCompartment$celltypes_age,levels = c("Bcell_Early.young","Bcell_Early.old","B_cell.young","B_cell.old","Bcell_Naive_Isg.young","Bcell_Naive_Isg.old","Plasma.young","Plasma.old","CD4_Early.young","CD4_Early.old","CD4_Naive.young","CD4_Naive.old","CD4.young","CD4.old","CD4_Cytotoxic.young","CD4_Cytotoxic.old","CD4_TEM.young","CD4_TEM.old","CD8.young","CD8.old","Th17.young","Th17.old","Monocyte.young","Monocyte.old","Macrophage.young","Macrophage.old","Neutrophil.young","Neutrophil.old"))

## reorder_FN_test_DotPlot_regardlessCompartments.pdf
ggplot(test2_regardlessCompartment,aes(x= cytokines,y= fct_rev(celltypes_age), fill = avgExp))+
geom_point(aes(size=ratio_all_age), shape = 21)+theme(text = element_text(size=10),axis.text.x = element_text(angle=90, hjust=1)) +scale_fill_gradient(low = "yellow", high = "red", na.value = NA)
## Figure-7E: reorder_FN_test_DotPlot_regardlessCompartments2.pdf
ggplot(test2_regardlessCompartment[test2_regardlessCompartment$cytokines %in% c("Ifng","Il1b"),],aes(x= celltypes_age,y= cytokines, fill = avgExp))+
geom_point(aes(size=ratio_all_age), shape = 21)+theme(text = element_text(size=10),axis.text.x = element_text(angle=90, hjust=1)) +scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+ theme(legend.position="bottom")

#### reorder_FN_test_DotPlot2.pdf
test2<-read.csv("./data_in/test2.csv",head=T)
test2$pos<-factor(test2$pos,levels = c("Cecum","C_P","C_D"))
test2$celltypes_age<-factor(test2$celltypes_age,levels = c("Bcell_Early.young","Bcell_Early.old","B_cell.young","B_cell.old","Bcell_Naive_Isg.young","Bcell_Naive_Isg.old","Plasma.young","Plasma.old","CD4_Early.young","CD4_Early.old","CD4_Naive.young","CD4_Naive.old","CD4.young","CD4.old","CD4_Cytotoxic.young","CD4_Cytotoxic.old","CD4_TEM.young","CD4_TEM.old","CD8.young","CD8.old","Th17.young","Th17.old","Monocyte.young","Monocyte.old","Macrophage.young","Macrophage.old","Neutrophil.young","Neutrophil.old"))
ggplot(test2,aes(x= cytokines,y= fct_rev(celltypes_age), fill = avgExp))+
geom_point(aes(size=ratio_all_pos_age), shape = 21)+facet_wrap(~test2$pos,ncol=1)+theme(text = element_text(size=10),axis.text.x = element_text(angle=90, hjust=1)) +scale_fill_gradient(low = "yellow", high = "red", na.value = NA)




## Figure-5B: Reorder_Figure5B.pdf
my_comparisons <- list(c("Cecum", "C_P"),c("C_P", "C_D"),c("Cecum", "C_D"))
ggbarplot(yy2[yy2$age=="young",], x = "pos", y = "ratio", add = "mean_sd",error.plot = "upper_errorbar",
palette = "jco", facet.by = "celltypes_new_redo3",
position = position_dodge(0.8))+
stat_compare_means(comparisons=my_comparisons,method = "t.test", label = "p.signif",label.y = c(1.1, 1.1, 1.2))
#### Figure S5A
logic_order<-c("Bcell-Early","B-cell","Bcell-Naive-Isg","Plasma","CD4-Early","CD4-Naive","CD4","CD4-Cytotoxic","CD4-TEM","CD8","Th17","Monocyte","Macrophage","Neutrophil")
mycolors3<-c("#84BF96","#6DBD57","#7F9D55","#3385BB","#DE9E83","#FBB268","#FE8D19","#F57C7C","#E42622","#891c1e","#F3E587","#977899","#9D7BBA","#B15928")


## Figure-S5A: Reorder_top20markers_heatmap2_rasterF2.pdf
Idents(myData_LI)<-myData_LI$celltypes_new_redo2
Idents(myData_LI)<-factor(Idents(myData_LI),levels = c("B-cell","CD4","CD4-TEM","CD8","Plasma","CD4-Naive","CD4-Early","Neutrophil","Th17","CD4-Cytotoxic","Bcell-Naive-Isg","Monocyte","Bcell-Early","Macrophage"))
top20_target.redo2$cluster<-factor(top20_target.redo2$cluster,levels = c("B-cell","CD4","CD4-TEM","CD8","Plasma","CD4-Naive","CD4-Early","Neutrophil","Th17","CD4-Cytotoxic","Bcell-Naive-Isg","Monocyte","Bcell-Early","Macrophage"))
top20_target.redo2<-top20_target.redo2[order(top20_target.redo2$cluster),]

pp<-DoHeatmap(myData_LI, features = top20_target.redo2$gene,size = 2,slot = 'data',raster=FALSE)
pp+ NoLegend()+ theme(axis.text.y = element_text(size = 1),
legend.title = element_text(size = 4),
legend.text = element_text(size = 8)
) + scale_fill_gradientn(colors = colorRampPalette(c("#1e1e1e","#ff0000"))(256))

#### Reorder_top20markers_heatmap2_rasterF.pdf
Idents(myData_LI)<-factor(Idents(myData_LI),levels = logic_order)
top20_target.redo2$cluster<-factor(top20_target.redo2$cluster,levels = logic_order)
top20_target.redo2<-top20_target.redo2[order(top20_target.redo2$cluster),]
pp<-DoHeatmap(myData_LI, features = top20_target.redo2$gene,size = 2,slot = 'data',raster=FALSE)
pp+ NoLegend()+ theme(axis.text.y = element_text(size = 1),
legend.title = element_text(size = 4),
legend.text = element_text(size = 8)
) + scale_fill_gradientn(colors = colorRampPalette(c("#1e1e1e","#ff0000"))(256))

mycolors3<-c("#84BF96","#6DBD57","#7F9D55","#3385BB","#DE9E83","#FBB268","#FE8D19","#F57C7C","#E42622","#891c1e","#F3E587","#977899","#9D7BBA","#B15928")
logic_order<-c("Bcell-Early","B-cell","Bcell-Naive-Isg","Plasma","CD4-Early","CD4-Naive","CD4","CD4-Cytotoxic","CD4-TEM","CD8","Th17","Monocyte","Macrophage","Neutrophil")
Idents(myData_LI)<-myData_LI$celltypes_new_redo2
Idents(myData_LI)<-factor(Idents(myData_LI),levels = logic_order)
top20_target.redo2$cluster<-factor(top20_target.redo2$cluster,levels = logic_order)
top20_target.redo2<-top20_target.redo2[order(top20_target.redo2$cluster),]
#### check known markers: /Users/jlu/Desktop/Pro_Mine/Paper_Atlas/Immune/figure/Reorder/immune_markers_check_vln/InPaper/
## Figure-S5B
VlnPlot(myData_LI, features = "Mki67", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Pou2af1", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Bank1", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Fcer2a", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Cd22", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Fcer2a", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Ifit3", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Isg15", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Stat1", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Adgre1", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Apoe", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Ccr2", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Ccrl2", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Cd14", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Cd19", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Cd68", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Cd160", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Crip1", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Ctla2a", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Gzma", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Gzmb", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Gzmk", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Herpud1", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Hmox1", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Ifitm3", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Igfbp4", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Cd8a", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Irf7", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Isg20", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Jaml", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Klrd1", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Lef1", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Lgmn", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Msrb1", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Plbd1", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Rorc", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "S100a4", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "S100a10", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Selenop", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Sell", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Stat1", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Tcf7", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Cd3g", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Cd3d", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Bank1", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Cd3", pt.size = 0)+scale_fill_manual(values = mycolors3)
VlnPlot(myData_LI, features = "Cd3e", pt.size = 0)+scale_fill_manual(values = mycolors3)
##
myData_LI$age_pos<-paste(myData_LI$age,myData_LI$pos,sep="#")
myData_LI$age_pos<-factor(myData_LI$age_pos,levels = c("young#Cecum","young#C_P","young#C_D","old#Cecum","old#C_P","old#C_D"))
myData_LI <- ScaleData(object = myData_LI, features = top20_target.redo2$gene, block.size = 2000)
#### /Users/jlu/Desktop/Pro_Mine/Paper_Atlas/Immune/figure/Reorder/
## SingleImmune_BcellMarkers_Compartment_2.pdf SingleImmune_TEMMarkers_Compartment_2.pdf SingleImmune_CD4Markers_Compartment_2.pdf SingleImmune_CD8Markers_Compartment_2.pdf
DoHeatmap(myData_LI, features = as.character(top20_target.redo2[top20_target.redo2$cluster=="B-cell",]$gene),size = 2,group.by = "age_pos",slot = "scale.data",raster=FALSE,group.bar.height = 0.02)+ scale_fill_gradientn(colors = colorRampPalette(c("#1e1e1e","#ff0000"))(6))
DoHeatmap(myData_LI, features = as.character(top20_target.redo2[top20_target.redo2$cluster=="CD4-TEM",]$gene),size = 2,group.by = "age_pos",slot = "scale.data",raster=FALSE,group.bar.height = 0.02)+ scale_fill_gradientn(colors = colorRampPalette(c("#1e1e1e","#ff0000"))(6))
DoHeatmap(myData_LI, features = as.character(top20_target.redo2[top20_target.redo2$cluster=="CD4",]$gene),size = 2,group.by = "age_pos",slot = "scale.data",raster=FALSE,group.bar.height = 0.02)+ scale_fill_gradientn(colors = colorRampPalette(c("#1e1e1e","#ff0000"))(6))
DoHeatmap(myData_LI, features = as.character(top20_target.redo2[top20_target.redo2$cluster=="CD8",]$gene),size = 2,group.by = "age_pos",slot = "scale.data",raster=FALSE,group.bar.height = 0.02)+ scale_fill_gradientn(colors = colorRampPalette(c("#1e1e1e","#ff0000"))(6))
DoHeatmap(myData_LI[myData_LI$age=="young",], features = as.character(top20_target.redo2[top20_target.redo2$cluster=="CD8",]$gene),size = 2,group.by = "age_pos",slot = "scale.data",raster=FALSE,group.bar.height = 0.02)+ scale_fill_gradientn(colors = colorRampPalette(c("#1e1e1e","#ff0000"))(6))
#### /Users/jlu/Desktop/Pro_Mine/Paper_Atlas/Immune/figure/Reorder/
## Figure-5D: SingleImmune_CD8Markers_Compartment_2y.pdf SingleImmune_CD4Markers_Compartment_2y.pdf SingleImmune_TEMMarkers_Compartment_2y.pdf SingleImmune_BcellMarkers_Compartment_2y.pdf
DoHeatmap(subset(myData_LI,subset=age=="young"), features = as.character(top20_target.redo2[top20_target.redo2$cluster=="CD8",]$gene),size = 2,group.by = "age_pos",slot = "scale.data",raster=FALSE,group.bar.height = 0.02)+ scale_fill_gradientn(colors = colorRampPalette(c("#1e1e1e","#ff0000"))(6))
DoHeatmap(subset(myData_LI,subset=age=="young"), features = as.character(top20_target.redo2[top20_target.redo2$cluster=="CD4",]$gene),size = 2,group.by = "age_pos",slot = "scale.data",raster=FALSE,group.bar.height = 0.02)+ scale_fill_gradientn(colors = colorRampPalette(c("#1e1e1e","#ff0000"))(6))
DoHeatmap(subset(myData_LI,subset=age=="young"), features = as.character(top20_target.redo2[top20_target.redo2$cluster=="CD4-TEM",]$gene),size = 2,group.by = "age_pos",slot = "scale.data",raster=FALSE,group.bar.height = 0.02)+ scale_fill_gradientn(colors = colorRampPalette(c("#1e1e1e","#ff0000"))(6))
DoHeatmap(subset(myData_LI,subset=age=="young"), features = as.character(top20_target.redo2[top20_target.redo2$cluster=="B-cell",]$gene),size = 2,group.by = "age_pos",slot = "scale.data",raster=FALSE,group.bar.height = 0.02)+ scale_fill_gradientn(colors = colorRampPalette(c("#1e1e1e","#ff0000"))(6))
##

#### barchart of summary of Figure 7A
## Figure-S7A: SingImmune_Pie_NormalizedByCelltype.pdf SingImmune_Pie_NormalizedByCompartment.pdf
write.csv(cellphonedb_all,"./data_out/cellphonedb_all.csv")
write.csv(cellphonedb_all_super,"./data_out/cellphonedb_all_super.csv")
#### normalized by celltypes
cellphonedb_all_pie<-cellphonedb_all
all.equal(cellphonedb_all_pie$TARGET_pos,cellphonedb_all_pie$SOURCE_pos)
all.equal(cellphonedb_all_pie$TARGET_age,cellphonedb_all_pie$SOURCE_age)
cellphonedb_all_pie<-cellphonedb_all_pie[cellphonedb_all_pie$TARGET_celltype_super=="Epithelial",]

cellphonedb_all_pie_ct<-aggregate(cellphonedb_all_pie$count, by=list(Category=cellphonedb_all_pie$pos.age.SOURCE_celltype_super), FUN=sum)
cellphonedb_all_pie_ct$SourceCT<-str_split_fixed(cellphonedb_all_pie_ct$Category,"\\.",3)[,3]
cellphonedb_all_pie_ct$pos.age<-paste(str_split_fixed(cellphonedb_all_pie_ct$Category,"\\.",3)[,1],str_split_fixed(cellphonedb_all_pie_ct$Category,"\\.",3)[,2],sep=".")
cellphonedb_all_pie_ct<-as.data.frame(dplyr::group_by(cellphonedb_all_pie_ct, SourceCT) %>% dplyr::mutate(percent = x/sum(x)))
cellphonedb_all_pie_ct$percent_x<-paste(cellphonedb_all_pie_ct$percent,cellphonedb_all_pie_ct$x,sep="_")
cellphonedb_all_pie_ct[cellphonedb_all_pie_ct$pos.age=="ProxCol.Old","pos.age"]<-"C_P.old"
cellphonedb_all_pie_ct[cellphonedb_all_pie_ct$pos.age=="ProxCol.Young","pos.age"]<-"C_P.young"
cellphonedb_all_pie_ct[cellphonedb_all_pie_ct$pos.age=="DistCol.Old","pos.age"]<-"C_D.old"
cellphonedb_all_pie_ct[cellphonedb_all_pie_ct$pos.age=="DistCol.Young","pos.age"]<-"C_D.young"
cellphonedb_all_pie_ct[cellphonedb_all_pie_ct$pos.age=="Cecum.Old","pos.age"]<-"Cecum.old"
cellphonedb_all_pie_ct[cellphonedb_all_pie_ct$pos.age=="Cecum.Young","pos.age"]<-"Cecum.young"
unique(cellphonedb_all_pie_ct$pos.age)
cellphonedb_all_pie_ct$pos.age<-factor(cellphonedb_all_pie_ct$pos.age,levels = c("Cecum.young","C_P.young","C_D.young","C_D.old","C_P.old","Cecum.old"))
mycols <- c("steelblue1", "seagreen1","tomato1","tomato3", "seagreen3", "steelblue3")
cellphonedb_all_pie_ct$SourceCT<-factor(cellphonedb_all_pie_ct$SourceCT,levels = c("Epithelial","T_cell_lineage","B_cell_lineage","Myeloid_lineage"))
ggplot(data=cellphonedb_all_pie_ct, aes(x="", y=round(percent,3), fill=pos.age)) +
geom_bar(stat="identity",color = "white")+
geom_text(aes(label=round(percent,3)), vjust=-0.3, color="black",
position = position_stack(vjust = 0.5), size=2)+facet_grid(~cellphonedb_all_pie_ct$SourceCT)+theme(text = element_text(size=8),axis.text.x=element_blank())+theme(axis.text.x=element_text(angle=90,vjust=0.5))+ coord_polar("y", start=0)+scale_fill_manual(values = mycols)+theme_void()

## Normalized by compartment
cellphonedb_all_pie_cm<-aggregate(cellphonedb_all_pie$count, by=list(Category=cellphonedb_all_pie$pos.age.SOURCE_celltype_super), FUN=sum)
cellphonedb_all_pie_cm$SourceCM<-str_split_fixed(cellphonedb_all_pie_cm$Category,"\\.",3)[,1]
cellphonedb_all_pie_cm$age<-str_split_fixed(cellphonedb_all_pie_cm$Category,"\\.",3)[,2]
cellphonedb_all_pie_cm[cellphonedb_all_pie_cm$age=="Young","age"]<-"young"
cellphonedb_all_pie_cm[cellphonedb_all_pie_cm$age=="Old","age"]<-"old"
cellphonedb_all_pie_cm$age.ct<-paste(cellphonedb_all_pie_cm$age,str_split_fixed(cellphonedb_all_pie_cm$Category,"\\.",3)[,3],sep=".")
cellphonedb_all_pie_cm$age.ct<-factor(cellphonedb_all_pie_cm$age.ct,levels = c("young.Epithelial","young.T_cell_lineage","young.B_cell_lineage","young.Myeloid_lineage","old.Myeloid_lineage","old.B_cell_lineage","old.T_cell_lineage","old.Epithelial"))

# cellphonedb_all_pie_cm<-as.data.frame(dplyr::group_by(cellphonedb_all_pie_cm, SourceCM) %>% dplyr::mutate(percent = x/sum(x)))
cellphonedb_all_pie_cm[cellphonedb_all_pie_cm$SourceCM=="DistCol","SourceCM"]<-"C_D"
cellphonedb_all_pie_cm[cellphonedb_all_pie_cm$SourceCM=="ProxCol","SourceCM"]<-"C_P"
cellphonedb_all_pie_cm<-as.data.frame(dplyr::group_by(cellphonedb_all_pie_cm, SourceCM) %>% dplyr::mutate(percent = x/sum(x)))
cellphonedb_all_pie_cm$SourceCM<-factor(cellphonedb_all_pie_cm$SourceCM,levels = c("Cecum","C_P","C_D"))
mycols2 <- c("skyblue1", "palegreen1","orange1","mediumpurple1", "mediumpurple3", "orange3","palegreen3","skyblue3")
ggplot(data=cellphonedb_all_pie_cm, aes(x="", y=round(percent,3), fill=age.ct)) +
geom_bar(stat="identity",color = "white")+
geom_text(aes(label=round(percent,3)), vjust=-0.3, color="black",
position = position_stack(vjust = 0.5), size=2)+facet_grid(~cellphonedb_all_pie_cm$SourceCM)+theme(text = element_text(size=8),axis.text.x=element_blank())+theme(axis.text.x=element_text(angle=90,vjust=0.5))+ coord_polar("y", start=0)+scale_fill_manual(values = mycols2)+theme_void()

## double check
cellphonedb_all_pie_cm[cellphonedb_all_pie_cm$SourceCM=="C_D","x"]
sum(cellphonedb_all_pie_cm[cellphonedb_all_pie_cm$SourceCM=="C_D","x"])
1790*0.12737430
sum(cellphonedb_all_pie_cm[cellphonedb_all_pie_cm$SourceCM=="Cecum","x"])
head(cellphonedb_all_pie_cm[cellphonedb_all_pie_cm$SourceCM=="Cecum","x"])
head(cellphonedb_all_pie_cm[cellphonedb_all_pie_cm$SourceCM=="Cecum","percent"])
2213*0.22413014

write.csv(cellphonedb_all_pie,"./data_out/cellphonedb_all_pie.csv")
write.csv(cellphonedb_all_pie_cm,"./data_out/cellphonedb_all_pie_cm.csv")
write.csv(cellphonedb_all_pie_ct,"./data_out/cellphonedb_all_pie_ct.csv")


### Figure 5C
my.data.redo2<-FetchData(myData_LI,c("orig.ident","sampleID","Index","nCount_RNA","nFeature_RNA","percent.mt","seurat_clusters","type","TYPE","TYPES","age","pos","celltypes_old","celltypes_new","celltypes_new_redo2"))
pca<-Embeddings(myData_LI, reduction = "pca")[, 1:20]
tsne<-Embeddings(myData_LI, reduction = "tsne")[, 1:2]
my.data.redo2<-merge(my.data.redo2,tsne,by="row.names",all.x=T)
rownames(my.data.redo2)<-my.data.redo2$Row.names
my.data.redo2<-my.data.redo2[,-1]
my.data.redo2<-merge(my.data.redo2,pca,by="row.names",all.x=T)
rownames(my.data.redo2)<-my.data.redo2$Row.names
my.data.redo2<-my.data.redo2[,-1]
write.csv(my.data.redo2,"./data_out//my.data.redo2.csv")
facet2<-function(my.data){
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
library(RColorBrewer)
mycolors <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(my.data$celltypes_new_redo2)))
tg3$Var2<-factor(tg3$Var2,levels = c("3#young","4#young","5#young","6#young","1#old","2#old","5#old","6#old"))
tg3$pos<-factor(tg3$pos,levels = c("Cecum","C_P","C_D"))
ggplot(tg3, aes(fill=celltypes_new_redo2, y=ratio, x=pos)) +
geom_bar(position="fill", stat="identity")+scale_fill_manual(values = mycolors)+facet_wrap(~Var2,nrow=2)+
theme(strip.text.x = element_text(size = 10),axis.text.x = element_text(angle=90, hjust=1))
return(tg3)
}
library(RColorBrewer)
mycolors <- colorRampPalette(brewer.pal(12,"Paired"))(length(unique(my.data$celltypes_new_redo2)))
tg42<-facet2(my.data.redo2)
tg42$age<-factor(tg42$age,levels=c("young","old"))
# only young
# Figure5C: SingleImmune_line_Reorder_pvalue_forPos_young_t-eachTwo.pdf
my_comparisons <- list(c("Cecum", "C_P"),c("C_P", "C_D"),c("Cecum", "C_D"))
ggline(tg42[tg42$age=="young",], x = "pos", y = "ratio", add = "mean_sd",color = "age", paltte = c("gray","gray0"),facet.by = "celltypes_new_redo2",nrow=2,ylim=c(0,0.8))+stat_compare_means(method = "anova", label.y = 0.8,size=2)+      # Add global p-value
    stat_compare_means(label = "p.signif", method = "t.test",
                     size=2,comparisons=my_comparisons,label.y = c(0.7, 0.7, 0.75),paired=F)   +theme(axis.text.x=element_text(angle=90,hjust=1))+scale_colour_manual(values = c("gray","gray0")) 
