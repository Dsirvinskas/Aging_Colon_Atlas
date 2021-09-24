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
levels(tenmeans[,2]) = c("Colonocyte_1", "Goblet", "Colonocyte_2", "Tuft-EE", "Stem-TA")
Idents(SI.integrated) = tenmeans[,2]




SI.integrated@meta.data = cbind(SI.integrated@meta.data, "Cluster" = SI.integrated@active.ident)
#MAST
      require(MAST)

SI.integrated@meta.data$Age[which(SI.integrated@meta.data$Age == "Young_20w")] = "Young"
SI.integrated@active.assay = "RNA"

ScaleData(SI.integrated, verbose = FALSE, do.scale = FALSE, do.center = FALSE)
 
      
      for (var in c("Cecum", "ProxCol", "DistCol")){
       DE1.Seurat = subset(SI.integrated, cells = which(SI.integrated@meta.data$Compartment == var))
        for (clust in c("Colonocyte_1", "Goblet", "Colonocyte_2", "Tuft-EE", "Stem-TA")){
          if (length(which(DE1.Seurat@meta.data$Cluster == clust)) > 0){
            DE.Seurat = subset(DE1.Seurat, idents = clust)
             if (length(which(DE.Seurat@meta.data$Age =="Old"))>=3 & length(which(DE.Seurat@meta.data$Age =="Young"))>=3){
                  Idents(DE.Seurat) = DE.Seurat@meta.data$Age
                  a = FindMarkers(object = DE.Seurat,test.use = "MAST", group.by = 'Age', ident.1 = "Young", ident.2 = "Old")
                  a = na.omit(a[a$p_val_adj <0.05,] )
                  
                      if(nrow(a)>0){
                  a = cbind(a, "Gene"=row.names(a))
                  a$Gene = as.character(a$Gene)
                  top = rbind(a %>% top_n(n = 10, wt = avg_logFC),a %>% top_n(n=-10, wt = avg_logFC))
                  
                  
                  write.csv(x = a ,file = paste(var,clust,"YvODE_MAST_genes.csv",sep=""))
                      }else{
                        print(paste("No significantly DE genes in", var, clust,sep=" "))
                      }
             } else {
               print(paste("Could not perform DE analysis on", var, clust, "due to missing cells",sep=" "))
             }
          } else {
            print(paste("Could not perform DE analysis on", var, clust, "due to missing cells",sep=" "))
        }
}
}
      
# Calculating compartment specific YvO 
for (var in c("Cecum", "ProxCol", "DistCol")){
  DE1.Seurat = subset(SI.integrated, cells = which(SI.integrated@meta.data$Compartment == var))
  
  a = FindMarkers(object = DE1.Seurat,test.use = "MAST", group.by = 'Age', ident.1 = "Young", ident.2 = "Old")

a = na.omit(a[a$p_val_adj <0.05,] )

if(nrow(a)>0){
  a = cbind(a, "Gene"=row.names(a))
  a$Gene = as.character(a$Gene)
  top = rbind(a %>% top_n(n = 10, wt = avg_logFC),a %>% top_n(n=-10, wt = avg_logFC))
  
  
  write.csv(x = a ,file = paste(var,"YvODE_MAST_genes.csv",sep=""))
}else{
  print(paste("No significantly DE genes in", var, sep=" "))
}
}


# Calculating Cell type specific YvO 

  DE1.Seurat = SI.integrated
  for (clust in c("Colonocyte_1", "Goblet", "Colonocyte_2", "Tuft-EE", "Stem-TA")){
    if (length(which(DE1.Seurat@meta.data$Cluster == clust)) > 0){
      DE.Seurat = subset(DE1.Seurat, idents = clust)
      if (length(which(DE.Seurat@meta.data$Age =="Old"))>=3 & length(which(DE.Seurat@meta.data$Age =="Young"))>=3){
        Idents(DE.Seurat) = DE.Seurat@meta.data$Age
        a = FindMarkers(object = DE.Seurat,test.use = "MAST", group.by = 'Age', ident.1 = "Young", ident.2 = "Old")
        a = na.omit(a[a$p_val_adj <0.05,] )
        
        if(nrow(a)>0){
          a = cbind(a, "Gene"=row.names(a))
          a$Gene = as.character(a$Gene)
          top = rbind(a %>% top_n(n = 10, wt = avg_logFC),a %>% top_n(n=-10, wt = avg_logFC))
          
          
          write.csv(x = a ,file = paste(clust,"YvODE_MAST_genes.csv",sep=""))
        }else{
          print(paste("No significantly DE genes in", clust,sep=" "))
        }
      } else {
        print(paste("Could not perform DE analysis on", var, clust, "due to missing cells",sep=" "))
      }
    } else {
      print(paste("Could not perform DE analysis on", var, clust, "due to missing cells",sep=" "))
    }
  }

#For Wilcox

  for (var in c("Cecum", "ProxCol", "DistCol")){
    DE1.Seurat = subset(SI.integrated, cells = which(SI.integrated@meta.data$Compartment == var))
    for (clust in c("Colonocyte_1", "Goblet", "Colonocyte_2", "Tuft-EE", "Stem-TA")){
      if (length(which(DE1.Seurat@meta.data$Cluster == clust)) > 0){
        DE.Seurat = subset(DE1.Seurat, idents = clust)
        if (length(which(DE.Seurat@meta.data$Age =="Old"))>=3 & length(which(DE.Seurat@meta.data$Age =="Young"))>=3){
          Idents(DE.Seurat) = DE.Seurat@meta.data$Age
          a = FindMarkers(object = DE.Seurat,test.use = "wilcox", group.by = 'Age', ident.1 = "Young", ident.2 = "Old")
          a = na.omit(a[a$p_val_adj <0.05,] )
          
          if(nrow(a)>0){
            a = cbind(a, "Gene"=row.names(a))
            a$Gene = as.character(a$Gene)
            top = rbind(a %>% top_n(n = 10, wt = avg_logFC),a %>% top_n(n=-10, wt = avg_logFC))
            
            
            write.csv(x = a ,file = paste(var,clust,"YvODE_wilcox_genes.csv",sep=""))
          }else{
            print(paste("No significantly DE genes in", var, clust,sep=" "))
          }
        } else {
          print(paste("Could not perform DE analysis on", var, clust, "due to missing cells",sep=" "))
        }
      } else {
        print(paste("Could not perform DE analysis on", var, clust, "due to missing cells",sep=" "))
      }
    }
  }
  
  # Calculating compartment specific YvO 
  for (var in c("Cecum", "ProxCol", "DistCol")){
    DE1.Seurat = subset(SI.integrated, cells = which(SI.integrated@meta.data$Compartment == var))
    
    a = FindMarkers(object = DE1.Seurat,test.use = "wilcox", group.by = 'Age', ident.1 = "Young", ident.2 = "Old")
    
    a = na.omit(a[a$p_val_adj <0.05,] )
    
    if(nrow(a)>0){
      a = cbind(a, "Gene"=row.names(a))
      a$Gene = as.character(a$Gene)
      top = rbind(a %>% top_n(n = 10, wt = avg_logFC),a %>% top_n(n=-10, wt = avg_logFC))
      
      
      write.csv(x = a ,file = paste(var,"YvODE_wilcox_genes.csv",sep=""))
    }else{
      print(paste("No significantly DE genes in", var, sep=" "))
    }
  }
  
  
  # Calculating Cell type specific YvO 
  
  DE1.Seurat = SI.integrated
  for (clust in c("Colonocyte_1", "Goblet", "Colonocyte_2", "Tuft-EE", "Stem-TA")){
    if (length(which(DE1.Seurat@meta.data$Cluster == clust)) > 0){
      DE.Seurat = subset(DE1.Seurat, idents = clust)
      if (length(which(DE.Seurat@meta.data$Age =="Old"))>=3 & length(which(DE.Seurat@meta.data$Age =="Young"))>=3){
        Idents(DE.Seurat) = DE.Seurat@meta.data$Age
        a = FindMarkers(object = DE.Seurat,test.use = "wilcox", group.by = 'Age', ident.1 = "Young", ident.2 = "Old")
        a = na.omit(a[a$p_val_adj <0.05,] )
        
        if(nrow(a)>0){
          a = cbind(a, "Gene"=row.names(a))
          a$Gene = as.character(a$Gene)
          top = rbind(a %>% top_n(n = 10, wt = avg_logFC),a %>% top_n(n=-10, wt = avg_logFC))
          
          
          write.csv(x = a ,file = paste(clust,"YvODE_wilcox_genes.csv",sep=""))
        }else{
          print(paste("No significantly DE genes in", clust,sep=" "))
        }
      } else {
        print(paste("Could not perform DE analysis on", var, clust, "due to missing cells",sep=" "))
      }
    } else {
      print(paste("Could not perform DE analysis on", var, clust, "due to missing cells",sep=" "))
    }
  }
  
  
  
  

save.image(file ="LI.RData")