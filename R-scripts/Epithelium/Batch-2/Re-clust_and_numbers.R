require(ggplot2)
require(dplyr)
require(Matrix)
require(cowplot)
require(Seurat)

load("/data/tmp/Dovydas/Single_Cell_RNAseq/2020_02_05_YvO_Single_Cell/Seurat/With_old_4_best_mice/Merged.RData")

## Gives a .csv file with the number of cells in young vs old, Cecum / ProxCol / DistCol, Initial compartments. 
c=5
tenmeans = read.csv(file = paste(c,"TESTallcels.csv", sep = ""))
tenmeans[,2] = as.factor(tenmeans[,2])
row.names(tenmeans) = tenmeans[,1]
tenmeans = tenmeans[row.names(tenmeans) %in% row.names(SI.integrated@meta.data),]

Idents(SI.integrated) = tenmeans[,2]

#Info = cbind(SI.integrated@meta.data, "Cluster" = SI.integrated@active.ident)
#Info$Age[which(Info$Age == "Young_20w")] = "Young"
#A = c(which(colnames(Info)=="Compartment"), which(colnames(Info)=="Age"), which(colnames(Info)=="Cluster"))
#levels(Info$Cluster) = c("Colonocyte_1", "Goblet", "Colonocyte_2", "Tuft-EE", "Stem-TA")

#Info1 = table(Info[,A])
#write.csv(Info1, file = "Compartment_Numbers.csv")



## Separating tuft and enteroendocrine cells
levels(SI.integrated@active.ident) = c("Colonocyte_1", "Goblet", "Colonocyte_2", "Tuft-EE", "Stem-TA")
SI.integrated@meta.data = cbind(SI.integrated@meta.data, "Orig_Cluster" = SI.integrated@active.ident)

  # For Tuft_EE cells
  #Tuft_EE = subset(SI.integrated, idents = "Tuft-EE")
  ## K-means
  #  dat = Tuft_EE@assays$integrated@data
  #  texpdf<-t(dat)
  #  tot<-data.frame(center=NULL,tot=NULL)
  #  c=2
  #    k<-kmeans(texpdf,center=c,nstart=5)
  #    write.csv(k$cluster,paste(c,"mean_Tuft_EE.csv",sep=""))
  #    print(c)
  #    
  #    tenmeans = read.csv(file = paste(c,"mean_Tuft_EE.csv", sep = ""))
  #    tenmeans[,2] = as.factor(tenmeans[,2])
  #    row.names(tenmeans) = tenmeans[,1]
  #    tenmeans = tenmeans[row.names(tenmeans) %in% row.names(Tuft_EE@meta.data),]
  
  #    Idents(Tuft_EE) = tenmeans[,2]
  #    cols = c("Red", "Blue", "grey30", "Green", "Yellow", "Magenta", "Cyan", "indianred", "Pink", "Purple", "orangered1", "wheat4", "darkolivegreen2", "slategray3", "lightsalmon2")
  #    pdf(file = paste("Tuft_EE_rPCA_TSNE_", c,"means.pdf", sep = "")) 
  #    plot = DimPlot(Tuft_EE, reduction = "tsne", cols = cols) 
  #    print(plot)
  #    dev.off() 
  
  #    Tuft_EE.markers <- FindAllMarkers(Tuft_EE, only.pos = F, min.pct = 0.1, logfc.threshold = 0.1) 
  #    Tuft_EE.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) 
  
  #    top10 <- Tuft_EE.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) 
  #    write.csv(top10, file= paste("Tuft_EE", c,"meanTop10_genes_rPCA.csv", sep ="") )
  
  #    col.plot = c("Red", "Black", "Blue") 
  #    pdf(file = paste("Tuft_EE_Heatmap_",c,"means_rPCA.pdf", sep = "") )
  #    plot = DoHeatmap(Tuft_EE, features = top10$gene, size=2, angle=0) + NoLegend() + scale_fill_gradient2( low = "Blue", mid = "Black", high = "Red", midpoint = 0, guide = "colourbar", aesthetics = "fill") +theme(axis.text = element_text(size=5)) 
  #    print(plot)
  #    dev.off() 
  
  #    Info = cbind(Tuft_EE@meta.data, "Cluster" = Tuft_EE@active.ident)
  #    Info$Age[which(Info$Age == "Young_20w")] = "Young"
  #    A = c(which(colnames(Info)=="Compartment"), which(colnames(Info)=="Age"), which(colnames(Info)=="Cluster"))
  
  #    Info1 = table(Info[,A])
  #    write.csv(Info1, file = "Tuft_EE_Compartment_Numbers.csv")


#    #For Stem-TA cells
#    Stem_TA = subset(SI.integrated, idents = "Stem-TA")
#    ## K-means
#    dat = Stem_TA@assays$integrated@data
#    texpdf<-t(dat)
#    tot<-data.frame(center=NULL,tot=NULL)
#    for (c in 6:9) {
#      k<-kmeans(texpdf,center=c,nstart=5)
#      write.csv(k$cluster,paste(c,"mean_Stem_TA.csv",sep=""))
#      print(c)
#      
#      tenmeans = read.csv(file = paste(c,"mean_Stem_TA.csv", sep = ""))
#      tenmeans[,2] = as.factor(tenmeans[,2])
#      row.names(tenmeans) = tenmeans[,1]
#      tenmeans = tenmeans[row.names(tenmeans) %in% row.names(Stem_TA@meta.data),]
#      
#      Idents(Stem_TA) = tenmeans[,2]
#      cols = c("Red", "Blue", "grey30", "Green", "Yellow", "Magenta", "Cyan", "indianred", "Pink", "Purple", "orangered1", "wheat4", "darkolivegreen2", "slategray3", "lightsalmon2")
#      pdf(file = paste("Stem-TA_rPCA_TSNE_", c,"means.pdf", sep = "")) 
#      plot = DimPlot(Stem_TA, reduction = "tsne", cols = cols) 
#      print(plot)
#      plot = DimPlot(Stem_TA, reduction = "tsne", cols = cols, group.by = "Compartment") 
#      print(plot)
#      plot = FeaturePlot(Stem_TA, reduction = "tsne", features = c("Lgr5", "Sox4", "Slc12a2", "Atoh1"), pt.size = 0.01, )
#      print(plot)
#      dev.off() 
#      
#      Stem_TA.markers <- FindAllMarkers(Stem_TA, only.pos = F, min.pct = 0.1, logfc.threshold = 0.1) 
#      Stem_TA.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) 
#      
#      top10 <- Stem_TA.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) 
#      write.csv(top10, file= paste("Stem_TA", c,"meanTop10_genes_rPCA.csv", sep ="") )
#      
#      col.plot = c("Red", "Black", "Blue") 
#      pdf(file = paste("Stem_TA_Heatmap_",c,"means_rPCA.pdf", sep = "") )
#      plot = DoHeatmap(Stem_TA, features = top10$gene, size=2, angle=0) + NoLegend() + scale_fill_gradient2( low = "Blue", mid = "Black", high = "Red", midpoint = 0, guide = "colourbar", aesthetics = "fill") +theme(axis.text = element_text(size=5)) 
#      print(plot)
#      dev.off() 
#      
#      Info = cbind(Stem_TA@meta.data, "Cluster" = Stem_TA@active.ident)
#      Info$Age[which(Info$Age == "Young_20w")] = "Young"
#      A = c(which(colnames(Info)=="Compartment"), which(colnames(Info)=="Age"), which(colnames(Info)=="Cluster"))
#      
#      Info1 = table(Info[,A])
#      write.csv(Info1, file = paste(c,"mean_Stem_TA_Compartment_Numbers.csv", sep=""))
#      
#      pdf(file= paste(c,"mean_Stem_TA_Violin.pdf"))
#      a = VlnPlot(Stem_TA, features = c("Lgr5", "Sox4", "Slc12a2", "Atoh1"), pt.size = 0.01)
#      print(a)
#      a = VlnPlot(Stem_TA, features = c("Lgr5", "Sox4", "Slc12a2", "Atoh1"), pt.size = 0)
#      print(a)
#      a = VlnPlot(Stem_TA, features = c("Mki67", "Pcna", "Dll1", "Ascl2"), pt.size = 0.01)
#      print(a)
#      a = VlnPlot(Stem_TA, features = c("Mki67", "Pcna", "Dll1", "Ascl2"), pt.size = 0)
#      print(a)
#      dev.off()
#      
#    }

## For Goblet cells
Goblet = subset(SI.integrated, idents = "Goblet")
## K-means
    dat = Goblet@assays$integrated@data
    texpdf<-t(dat)
    tot<-data.frame(center=NULL,tot=NULL)
   for (c in 2:6) {
      k<-kmeans(texpdf,center=c,nstart=5)
      write.csv(k$cluster,paste(c,"mean_Goblet.csv",sep=""))
      print(c)
      
      tenmeans = read.csv(file = paste(c,"mean_Goblet.csv", sep = ""))
      tenmeans[,2] = as.factor(tenmeans[,2])
      row.names(tenmeans) = tenmeans[,1]
      tenmeans = tenmeans[row.names(tenmeans) %in% row.names(Goblet@meta.data),]
      
      Idents(Goblet) = tenmeans[,2]
      cols = c("Red", "Blue", "grey30", "Green", "Yellow", "Magenta", "Cyan", "indianred", "Pink", "Purple", "orangered1", "wheat4", "darkolivegreen2", "slategray3", "lightsalmon2")
      pdf(file = paste("Goblet_rPCA_TSNE_", c,"means.pdf", sep = "")) 
      plot = DimPlot(Goblet, reduction = "tsne", cols = cols) 
      print(plot)
      plot = DimPlot(Goblet, reduction = "tsne", cols = cols, group.by = "Compartment") 
      print(plot)
      plot = FeaturePlot(Goblet, reduction = "tsne", features = c("Muc2", "Spink1", "Reg4", "Atoh1"), pt.size = 0.01, )
      print(plot)
      dev.off() 
      
      Goblet.markers <- FindAllMarkers(Goblet, only.pos = F, min.pct = 0.1, logfc.threshold = 0.1) 
      Goblet.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) 
      
      top10 <- Goblet.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) 
      write.csv(top10, file= paste("Goblet", c,"meanTop10_genes_rPCA.csv", sep ="") )
      
      col.plot = c("Red", "Black", "Blue") 
      pdf(file = paste("Goblet_Heatmap_",c,"means_rPCA.pdf", sep = "") )
      plot = DoHeatmap(Goblet, features = top10$gene, size=2, angle=0) + NoLegend() + scale_fill_gradient2( low = "Blue", mid = "Black", high = "Red", midpoint = 0, guide = "colourbar", aesthetics = "fill") +theme(axis.text = element_text(size=5)) 
      print(plot)
      dev.off() 
      
      Info = cbind(Goblet@meta.data, "Cluster" = Goblet@active.ident)
      Info$Age[which(Info$Age == "Young_20w")] = "Young"
      A = c(which(colnames(Info)=="Compartment"), which(colnames(Info)=="Age"), which(colnames(Info)=="Cluster"))
      
      Info1 = table(Info[,A])
      write.csv(Info1, file = paste(c,"mean_Goblet_Compartment_Numbers.csv", sep=""))
      
    pdf(file= paste(c,"mean_Goblet_Violin.pdf"))
    a = VlnPlot(Goblet, features = c("Lgr5", "Sox4", "Slc12a2","Ascl2","Mki67", "Pcna"), pt.size = 0.01)
    print(a)
    a = VlnPlot(Goblet, features = c("Lgr5", "Sox4", "Slc12a2","Ascl2","Mki67", "Pcna"), pt.size = 0)
    print(a)
    
    a = VlnPlot(Goblet, features = c("Dll1", "Atoh1","Chga", "Chgb"), pt.size = 0.01)
    print(a)
    a = VlnPlot(Goblet, features = c("Dll1", "Atoh1","Chga", "Chgb"), pt.size = 0)
    print(a)
    
    a = VlnPlot(Goblet, features = c("Reg4", "Muc2", "Agr2", "Spink1","Cd24a", "Dclk1"), pt.size = 0.01)
    print(a)
    a = VlnPlot(Goblet, features = c("Reg4", "Muc2", "Agr2", "Spink1","Cd24a", "Dclk1"), pt.size = 0)
    print(a)
    
    a = VlnPlot(Goblet, features = c("Ptprc", "H2-aa", "H2-ab1"), pt.size = 0.01)
    print(a)
    a = VlnPlot(Goblet, features = c("Ptprc", "H2-aa", "H2-ab1"), pt.size = 0)
    print(a)
    
    dev.off()
     
   }