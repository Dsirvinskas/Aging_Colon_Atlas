
########################
#### Cellphonedb #######
########################

library(dplyr)
library(Seurat)
library(janitor)
library(stringr)

myData_LI.VariableFeatures<-readRDS("data_in/myData_LI.VariableFeatures.rds")
feature.names = read.delim("data_in/features.tsv.gz",header = FALSE,stringsAsFactors = FALSE)
myData_LI.count_norm<-read.table("data_in/Immune_counts.txt",sep="\t",header = T,stringsAsFactors = F)
myData_LI.count_norm$ENSMUSG<-feature.names[match(myData_LI.count_norm$Gene,feature.names$V2),"V1"]
myData_LI.count_norm <- myData_LI.count_norm %>% dplyr::select(ENSMUSG, dplyr::everything())
#### Basic function to convert mouse to human gene names in ENSMUSG
convertMouseGeneList <- function(x){
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = x , mart = mouse, attributesL = c("ensembl_gene_id"), martL = human, uniqueRows=T)
humanx <- unique(genesV2[, 2])
# Print the first 6 genes found to the screen
print(head(humanx))
return(genesV2)
}
aa<-convertMouseGeneList(myData_LI.count_norm$ENSMUSG)
bb<-sapply(myData_LI.count_norm$ENSMUSG,function(x){ifelse(x %in% aa$Gene.stable.ID,aa[aa$Gene.stable.ID==x,"Gene.stable.ID.1"],x)})
bb<-as.data.frame(bb)
myData_LI.count_norm$bb<-bb$bb
myData_LI.count_norm <- myData_LI.count_norm %>% dplyr::select(bb, dplyr::everything())
myData_LI.count_norm.ensembl.h<-myData_LI.count_norm
myData_LI.count_norm.ensembl.h<-myData_LI.count_norm.ensembl.h[,-c(2,3)]
colnames(myData_LI.count_norm.ensembl.h)[1]<-"Gene"
myData_LI.count_norm.ensembl.h.var<-myData_LI.count_norm[myData_LI.count_norm$Gene %in% myData_LI.VariableFeatures,]
myData_LI.count_norm.ensembl.h.var<-myData_LI.count_norm.ensembl.h.var[,-c(2,3)]
write.table(myData_LI.count_norm.ensembl.h.var, './data_in/myData_LI.count_norm.ensembl.h.var.txt', sep='\t', quote=F,row.names = F)
write.table(myData_LI.count_norm.ensembl.h, './data_in/myData_LI.count_norm.ensembl.h.txt', sep='\t', quote=F,row.names = F)
myData_LI.meta<-read.table("data_in/Immune_meta.txt",head=T,stringsAsFactors = F)
all.equal(colnames(myData_LI.count_norm.ensembl.h),gsub("-",".",myData_LI.meta$Cell))
unique(colnames(myData_LI.count_norm.ensembl.h)[-1]==gsub("-",".",myData_LI.meta$Cell))
saveRDS(myData_LI.count_norm,"./data_in/myData_LI.count_norm.rds")

##############################
#### Epithelial in Odin ######
##############################
setwd("/u2/tmp/Jing/PIPs/CellphoneDB/epithelial_Dovydas")
load("Merged.RData")
#### Loading cluster information for different cell populations
Initial = read.csv(file = "5kmeans_initial.csv") # Initial 5 k-means 
StemTa = read.csv(file = "4kmeans_Stem_TA.csv") # 4 k means of Stem/TA cluster
TuftEE = read.csv(file = "2kmeans_Tuft_EE.csv") # 2 k means of Tuft/EE cluster
#### Overlaying the initial 5 k-means clustering info on cells
tenmeans = Initial
tenmeans[,2] = as.factor(tenmeans[,2])
row.names(tenmeans) = tenmeans[,1]
tenmeans = tenmeans[row.names(tenmeans) %in% row.names(SI.integrated@meta.data),]
Idents(SI.integrated) = tenmeans[,2]
levels(SI.integrated@active.ident) = c("Colonocyte_1", "Goblet", "Colonocyte_2", "Tuft-EE", "Stem-TA")
SI.integrated.counts_norm<-SI.integrated[["RNA"]]@counts
feature.names = read.delim("/data/tmp/Jing/PIPs/CellphoneDB/epithelial_Dovydas/features.tsv.gz",header = FALSE,stringsAsFactors = FALSE)
SI.integrated.counts_norm <- apply(SI.integrated.counts_norm, 2, function(x) (x/sum(x))*10000)
SI.integrated.counts_norm<-as.data.frame(SI.integrated.counts_norm)
SI.integrated.counts_norm$Symbol<-row.names(SI.integrated.counts_norm)
library(dplyr)
SI.integrated.counts_norm <- SI.integrated.counts_norm %>% dplyr::select(Symbol, dplyr::everything())
SI.integrated.counts_norm$ENSMUSG<-feature.names[match(SI.integrated.counts_norm$Symbol,feature.names$V2),"V1"]
SI.integrated.counts_norm <- SI.integrated.counts_norm %>% dplyr::select(ENSMUSG, dplyr::everything())
#### Basic function to convert mouse to human gene names in ENSMUSG
convertMouseGeneList <- function(x){
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = x , mart = mouse, attributesL = c("ensembl_gene_id"), martL = human, uniqueRows=T)
humanx <- unique(genesV2[, 2])
#### Print the first 6 genes found to the screen
print(head(humanx))
return(genesV2)
}
aa<-convertMouseGeneList(SI.integrated.counts_norm$ENSMUSG)
bb<-sapply(SI.integrated.counts_norm$ENSMUSG,function(x){ifelse(x %in% aa$Gene.stable.ID,aa[aa$Gene.stable.ID==x,"Gene.stable.ID.1"],x)})
bb<-as.data.frame(bb)
SI.integrated.counts_norm$Gene<-bb$bb
SI.integrated.counts_norm <- SI.integrated.counts_norm %>% dplyr::select(Gene, dplyr::everything())
SI.integrated.counts_norm$Gene<-as.character(SI.integrated.counts_norm$Gene)
SI.integrated.counts_norm<-SI.integrated.counts_norm[!is.na(SI.integrated.counts_norm$Gene),]
epithelial_meta<-read.csv("meta-data-forJing.csv",head=T,stringsAsFactors=F)
all.equal(epithelial_meta$X,colnames(SI.integrated.counts_norm)[-c(1,2,3)])
epithelial_meta$cell_type<-paste(epithelial_meta$Compartment,epithelial_meta$Age,epithelial_meta$Recluster,sep=".")
colnames(epithelial_meta)[1]<-"Cell"
SI.integrated.VariableFeatures<-VariableFeatures(SI.integrated)
write.table(SI.integrated.counts_norm[SI.integrated.counts_norm$Symbol %in% SI.integrated.VariableFeatures,-c(2,3)],"SI.integrated.counts_norm_variable.txt", sep='\t', quote=F,row.names = F)
write.table(SI.integrated.counts_norm[,-c(2,3)],"SI.integrated.counts_norm.txt", sep='\t', quote=F,row.names = F)
saveRDS(SI.integrated.VariableFeatures,"SI.integrated.VariableFeatures.rds")

###############################################
### combine epithelial and immune in Odin #####
###############################################

myData_LI.count_norm <-readRDS("/data/tmp/Jing/PIPs/CellphoneDB/myData_LI.count_norm.rds")
colnames(myData_LI.count_norm)[3]<-"Symbol"
colnames(myData_LI.count_norm)[1]<-"Gene"
IE.counts_norm<-merge(myData_LI.count_norm,SI.integrated.counts_norm,by="Symbol",all=T)

myData_LI.VariableFeatures <- readRDS("../myData_LI.VariableFeatures.rds")
IE.VariableFeatures<-unique(c(myData_LI.VariableFeatures,SI.integrated.VariableFeatures))
#### IE.counts_norm.2: variable genes in combined data
IE.counts_norm.2<-IE.counts_norm[IE.counts_norm$Symbol %in% IE.VariableFeatures,]
IE.counts_norm.2 <- IE.counts_norm.2[,-which(names(IE.counts_norm.2) %in% c("Gene.y","ENSMUSG.y","ENSMUSG.x","Symbol"))]
colnames(IE.counts_norm.2)[1]<-"Gene"
IE.counts_norm.2$Gene<-as.character(IE.counts_norm.2$Gene)
IE.counts_norm.2<-IE.counts_norm.2[!is.na(IE.counts_norm.2[,1]),]
IE.counts_norm.2[is.na(IE.counts_norm.2)]<-0
#### meta for conbined data
immune_meta<-read.table("../Immune_meta.txt",head=T,stringsAsFactors=F)
IE_meta<-rbind(immune_meta,epithelial_meta[,c(1,12)])
IE_meta$Cell<-gsub("\\.","-",IE_meta$Cell)
IE_meta <- IE_meta[match(colnames(IE.counts_norm.2)[-1],IE_meta$Cell),]

write.table(IE_meta,"IE_meta.txt", sep='\t', quote=F,row.names = F)
write.table(IE.counts_norm.2,"IE.counts_norm.2.txt", sep='\t', quote=F,row.names = F)
savehistory("cellphonedb_pip4.Rhistory")

################################################################
### Plot Plot Plot #############################################
################################################################

##### interaction
ie.interaction<-read.table("./data_out/all_0.1_out_var_HM/significant_means.txt",head=T,sep="\t",stringsAsFactors = F)
ie.interaction<-ie.interaction[ie.interaction$is_integrin=="False",]
ie.interaction<-t(ie.interaction)
ie.interaction<-as.data.frame(ie.interaction)

ie.interaction<-ie.interaction[-1,]
ie.interaction<-ie.interaction[-c(2,3,4,5,6,7,8,9,10,11),]
names(ie.interaction)<-as.matrix(ie.interaction[1,])
ie.interaction<-ie.interaction[-1,]
ie.interaction[] <- lapply(ie.interaction, function(x) type.convert(as.character(x)))
colheader<-read.delim2("./data_out/all_0.1_out_var_HM/colheader",head=T,sep="\t",stringsAsFactors = F)
all.equal(colheader$name2,row.names(ie.interaction))
unique(colheader$name2==row.names(ie.interaction))
unique(gsub("-",".",colheader$name2)==row.names(ie.interaction))
all.equal(gsub("-",".",colheader$name2),row.names(ie.interaction))

colheader$compartment1<-str_split_fixed(colheader$celltype1,"\\.",3)[,1]
colheader$age1<-str_split_fixed(colheader$celltype1,"\\.",3)[,2]
colheader$compartment2<-str_split_fixed(colheader$celltype2,"\\.",3)[,1]
colheader$age2<-str_split_fixed(colheader$celltype2,"\\.",3)[,2]

colheader[colheader$compartment1=="DistCol","compartment1"]<-"C_D"
colheader[colheader$compartment1=="ProxCol","compartment1"]<-"C_P"
colheader[colheader$compartment2=="ProxCol","compartment2"]<-"C_P"
colheader[colheader$compartment2=="DistCol","compartment2"]<-"C_D"
colheader[colheader$age1=="Young","age1"]<-"young"
colheader[colheader$age1=="Old","age1"]<-"old"
colheader[colheader$age2=="Young","age2"]<-"young"
colheader[colheader$age2=="Old","age2"]<-"old"
ie.interaction.2<-ie.interaction[(colheader$compartment1==colheader$compartment2) & (colheader$age1==colheader$age2) & (colheader$age1 != "Young_20w"),]

ie.interaction.2$colheader<-row.names(ie.interaction.2)
ie.interaction.2$colheader<-gsub("ProxCol","C_P",ie.interaction.2$colheader)
ie.interaction.2$colheader<-gsub("DistCol","C_D",ie.interaction.2$colheader)
ie.interaction.2$colheader<-gsub("Young","young",ie.interaction.2$colheader)
ie.interaction.2$colheader<-gsub("Old","old",ie.interaction.2$colheader)
ie.interaction.2$compartment<-str_split_fixed(ie.interaction.2$colheader,"\\.",3)[,1]
ie.interaction.2$age<-str_split_fixed(ie.interaction.2$colheader,"\\.",3)[,2]

colheader$subcelltype1<-str_split_fixed(colheader$celltype1,"\\.",3)[,3]
colheader$subcelltype2<-str_split_fixed(colheader$celltype2,"\\.",3)[,3]
colheader.2<-colheader[(colheader$compartment1==colheader$compartment2) & (colheader$age1==colheader$age2) & (colheader$age1 != "Young_20w"),]

all.equal(gsub("-",".",colheader.2$name2),row.names(ie.interaction.2))
ie.interaction.2$subcelltype1<-colheader.2$subcelltype1
ie.interaction.2$subcelltype2<-colheader.2$subcelltype2
ie.interaction.2$subcelltypes<-paste(ie.interaction.2$subcelltype1,ie.interaction.2$subcelltype2,sep="#")


ie.interaction.2$CountNoNA<-rowSums(!is.na(ie.interaction.2[,c(1:112)]))
ie.interaction.3<-ie.interaction.2[,c(1:112)]
ie.interaction.3[is.na(ie.interaction.3)]<-0
ie.interaction.2$SumMeans<-rowSums(ie.interaction.3)
ie.interaction.2$AvgMeans<-ie.interaction.2$SumMeans/ie.interaction.2$CountNoNA

bar_plot<-function(ie.interaction.2,age="young",compartment="Cecum",top=20,measurement="CountNoNA",size=20){
ie.interaction.2.tmp<-ie.interaction.2[ie.interaction.2$age==age & ie.interaction.2$compartment==compartment,]
print(unique(ie.interaction.2.tmp[,"age"]))
print(unique(ie.interaction.2.tmp[,"compartment"]))
print(measurement)
print(head(ie.interaction.2.tmp[,measurement]))
p<-ggplot(head(ie.interaction.2.tmp[order(-ie.interaction.2.tmp[,measurement]),],n=top), aes(x =reorder(subcelltypes,get(measurement)), y = get(measurement))) + geom_bar(stat = "identity")+coord_flip()+theme(text = element_text(size=size))+
xlab(paste("top20",age,compartment,sep="_"))+ylab(measurement)
print(p)
}

bar_plot(ie.interaction.2)
bar_plot(ie.interaction.2,age="young_20w",compartment = "C_P",measurement = "AvgMeans")


pdf("all_interaction.pdf")
for(i in c("young","old")){for(j in c("Cecum","C_P","C_D")){for(m in c("CountNoNA","SumMeans","AvgMeans")){bar_plot(ie.interaction.2,age=i,compartment = j,measurement = m,top=30,size=10)}}}
dev.off()


immune_c<-c(
"B-cell",
"Bcell-Early",
"Bcell-Naive-Isg",
"CD4",
"CD4-Cytotoxic",
"CD4-Early",
"CD4-Naive",
"CD4-TEM",
"CD8",
"Macrophage",
"Monocyte",
"Neutrophil",
"Plasma",
"Th17")
epithelial_c<-c(
"Cecum-enriched-Colonocyte",
"Colonocyte_Precursor",
"DistCol-enriched-Colonocyte",
"EnteroEndocrine_Cell",
"Goblet_Precursor",
"ProxCol-enriched-Colonocyte",
"Stem_cell",
"TA-cell",
"Tuft_Cell")

ie.interaction.ie<-ie.interaction.2[(ie.interaction.2$subcelltype1 %in% immune_c & ie.interaction.2$subcelltype2 %in% epithelial_c)|(ie.interaction.2$subcelltype1 %in% epithelial_c & ie.interaction.2$subcelltype2 %in% immune_c),]
savehistory("~/Desktop/Pro_Mine/Paper_Atlas/Immune/script/cellphonedb_pip5.Rhistory")


#### immune vs epithelial
pdf("all_interaction_ie_withoutY20w.pdf")
for(i in c("young","old")){for(j in c("Cecum","C_P","C_D")){for(m in c("CountNoNA","SumMeans","AvgMeans")){bar_plot(ie.interaction.2.ie,age=i,compartment = j,measurement = m,top=30,size=10)}}}
dev.off()
#### immune vs immune and epithelial vs epithelial
ie.interaction.2.ii<-ie.interaction.2[(ie.interaction.2$subcelltype1 %in% immune_c & ie.interaction.2$subcelltype2 %in% immune_c),]
ie.interaction.2.ee<-ie.interaction.2[(ie.interaction.2$subcelltype1 %in% epithelial_c & ie.interaction.2$subcelltype2 %in% epithelial_c),]
pdf("all_interaction_ii_withoutY20w.pdf")
for(i in c("young","old")){for(j in c("Cecum","C_P","C_D")){for(m in c("CountNoNA","SumMeans","AvgMeans")){bar_plot(ie.interaction.2.ii,age=i,compartment = j,measurement = m,top=30,size=10)}}}
dev.off()
pdf("all_interaction_ee_withoutY20w.pdf")
for(i in c("young","old")){for(j in c("Cecum","C_P","C_D")){for(m in c("CountNoNA","SumMeans","AvgMeans")){bar_plot(ie.interaction.2.ee,age=i,compartment = j,measurement = m,top=30,size=10)}}}
dev.off()

bar_plot_2<-function(ie.interaction.2,age="young",compartment="Cecum",top=20,measurement="CountNoNA",size=20){
ie.interaction.2.tmp<-ie.interaction.2[ie.interaction.2$age==age & ie.interaction.2$compartment==compartment,]
print(paste(unique(ie.interaction.2.tmp[,"age"]),unique(ie.interaction.2.tmp[,"compartment"]),measurement,sep="_"))
p<-head(ie.interaction.2.tmp[order(-ie.interaction.2.tmp[,measurement]),],n=top)
print("subcelltype1:")
print(table(p[,"subcelltype1"]))
print("subcelltype2:")
print(table(p[,"subcelltype2"]))
}
for(i in c("young","old")){for(j in c("Cecum","C_P","C_D")){for(m in c("CountNoNA","SumMeans","AvgMeans")){bar_plot_2(ie.interaction.2.ie,age=i,compartment = j,measurement = m,top=30,size=10)}}}
for(i in c("young","old")){for(j in c("Cecum","C_P","C_D")){for(m in c("CountNoNA")){bar_plot_2(ie.interaction.2.ie,age=i,compartment = j,measurement = m,top=30,size=10)}}}


#### check gene expression
IE.counts_norm.2<-read.table("/Users/jlu/Desktop/Pro_Mine/Paper_Atlas/Immune/data_out/all_0.1_out_var_HM/IE.counts_norm.2.txt",head=T,sep="\t",stringsAsFactors = F)
IE_meta<-read.table("/Users/jlu/Desktop/Pro_Mine/Paper_Atlas/Immune/data_out/all_0.1_out_var_HM/IE_meta.txt",head=T,sep="\t",stringsAsFactors = F)
all.equal(gsub("-",".",IE_meta$Cell),colnames(IE.counts_norm.2)[-1])
colnames(IE.counts_norm.2)[-1]<-IE_meta$cell_type
row.names(IE.counts_norm.2)<-IE.counts_norm.2$Gene
test_exp<-IE.counts_norm.2[IE.counts_norm.2$Gene=="ENSG00000019582",]
row.names(test_exp)<-test_exp[,1]
test_exp<-test_exp[,-1]
test_exp_2<-as.data.frame(t(test_exp))
test_exp_2$id<-colnames(test_exp)
test_exp_2$id<-colnames(IE.counts_norm.2)[-1]
test_exp_2$compartment<-str_split_fixed(test_exp_2$id,"\\.",3)[,1]
test_exp_2$age<-str_split_fixed(test_exp_2$id,"\\.",3)[,2]
test_exp_2$celltype<-str_split_fixed(test_exp_2$id,"\\.",3)[,3]
test_exp_2[test_exp_2$age=="Young","age"]<-"young"
test_exp_2[test_exp_2$age=="Old","age"]<-"old"
test_exp_2[test_exp_2$compartment=="ProxCol","compartment"]<-"C_P"
test_exp_2[test_exp_2$compartment=="DistCol","compartment"]<-"C_D"
test_exp_2[test_exp_2$celltype %in% immune_c,"superC"]<-"Immune"
test_exp_2[test_exp_2$celltype %in% epithelial_c,"superC"]<-"Epithelial"
test_exp_2$age<-factor(test_exp_2$age,levels = c("young","Young_20w","old"))
test_exp_2$compartment<-factor(test_exp_2$compartment,levels = c("Cecum","C_P","C_D"))
p<-ggplot(test_exp_2[!is.na(test_exp_2$superC),], aes(x=celltype, y=ENSG00000019582, fill=age)) + geom_violin() +facet_grid(compartment~superC,scales='free_x')+theme(text = element_text(size=10),axis.text.x = element_text(angle=90, hjust=1))
p+ geom_dotplot(binaxis='y', stackdir='center',dotsize=0.1,position=position_dodge(1))

####
exp_check<-function(IE.counts_norm.2,gene="ENSG00000019582"){
test_exp<-IE.counts_norm.2[IE.counts_norm.2$Gene==gene,]
row.names(test_exp)<-test_exp[,1]
test_exp<-test_exp[,-1]
test_exp_2<-as.data.frame(t(test_exp))
test_exp_2$id<-colnames(IE.counts_norm.2)[-1]
test_exp_2$compartment<-str_split_fixed(test_exp_2$id,"\\.",3)[,1]
test_exp_2$age<-str_split_fixed(test_exp_2$id,"\\.",3)[,2]
test_exp_2$celltype<-str_split_fixed(test_exp_2$id,"\\.",3)[,3]
test_exp_2[test_exp_2$age=="Young","age"]<-"young"
test_exp_2[test_exp_2$age=="Old","age"]<-"old"
test_exp_2[test_exp_2$compartment=="ProxCol","compartment"]<-"C_P"
test_exp_2[test_exp_2$compartment=="DistCol","compartment"]<-"C_D"
test_exp_2[test_exp_2$celltype %in% immune_c,"superC"]<-"Immune"
test_exp_2[test_exp_2$celltype %in% epithelial_c,"superC"]<-"Epithelial"
test_exp_2$age<-factor(test_exp_2$age,levels = c("young","Young_20w","old"))
test_exp_2$compartment<-factor(test_exp_2$compartment,levels = c("Cecum","C_P","C_D"))
p<-ggplot(test_exp_2[!is.na(test_exp_2$superC),], aes(x=celltype, y=get(gene), fill=age)) + geom_violin() +facet_grid(compartment~superC,scales='free_x')+theme(text = element_text(size=10),axis.text.x = element_text(angle=90, hjust=1))
p+ geom_dotplot(binaxis='y', stackdir='center',dotsize=0.1,position=position_dodge(1))+ylab(paste("normalized_exp",gene,sep=" "))
}
#### Check: gene expression check
exp_check(IE.counts_norm.2,gene="ENSG00000142192")
savehistory("~/Desktop/Pro_Mine/Paper_Atlas/Immune/script/cellphonedb_pip6.Rhistory")
for(i in c("young","old")){for(j in c("Cecum","C_P","C_D")){for(m in c("AvgMeans")){bar_plot_2(ie.interaction.2.ie,age=i,compartment = j,measurement = m,top=3000000,size=10)}}}
for(i in c("young","old")){for(j in c("Cecum","C_P","C_D")){for(m in c("AvgMeans")){bar_plot_2(ie.interaction.2.ie,age=i,compartment = j,measurement = m,top=30,size=10)}}}

#### Check: ligand-receptor partner expression level for B-cell ligand
unique(test_exp_2$celltype)
test_B<-ie.interaction.2.ie[ie.interaction.2.ie$subcelltype1=="B-cell",]
subs_B<-test_B[colSums(!is.na(test_B)) > 0]
test_B<-ie.interaction.2.ie[ie.interaction.2.ie$subcelltype1=="Bcell-Naive-Isg",]
subs_B<-test_B[colSums(!is.na(test_B)) > 0]
test_B<-ie.interaction.2.ie[ie.interaction.2.ie$subcelltype1=="Bcell-Naive-Isg"&ie.interaction.2.ie$age=="young",]
subs_B<-test_B[colSums(!is.na(test_B)) > 0]
test_B<-ie.interaction.2.ie[ie.interaction.2.ie$subcelltype1=="Monocyte"&ie.interaction.2.ie$age=="young",]
subs_B<-test_B[colSums(!is.na(test_B)) > 0]
test_B<-ie.interaction.2.ie[ie.interaction.2.ie$subcelltype1=="Bcell-Early"&ie.interaction.2.ie$age=="young",]
subs_B<-test_B[colSums(!is.na(test_B)) > 0]
test_B<-ie.interaction.2.ie[ie.interaction.2.ie$subcelltype1=="Bcell-Early"&ie.interaction.2.ie$age=="old",]
subs_B<-test_B[colSums(!is.na(test_B)) > 0]
test_B<-ie.interaction.2.ie[ie.interaction.2.ie$subcelltype1=="CD4-Early"&ie.interaction.2.ie$age=="young",]
subs_B<-test_B[colSums(!is.na(test_B)) > 0]
test_B<-ie.interaction.2.ie[ie.interaction.2.ie$subcelltype1=="Macrophage"&ie.interaction.2.ie$age=="young",]
subs_B<-test_B[colSums(!is.na(test_B)) > 0]
test_B<-ie.interaction.2.ie[ie.interaction.2.ie$subcelltype1=="B-cell"&ie.interaction.2.ie$age=="young",]
subs_B<-test_B[colSums(!is.na(test_B)) > 0]
test_B<-ie.interaction.2.ii[ie.interaction.2.ii$subcelltype1=="B-cell"&ie.interaction.2.ii$age=="young",]
subs_B<-test_B[colSums(!is.na(test_B)) > 0]
test_B<-ie.interaction.2[ie.interaction.2$subcelltype1=="B-cell"&ie.interaction.2$age=="young"&ie.interaction.2$compartment=="Cecum",]
subs_B<-test_B[colSums(!is.na(test_B)) > 0]
subs_Bx<-subs_B[subs_B$subcelltype1!="none" | subs_B$subcelltype2!="none",]
subs_B<-subs_B[!grepl("none",subs_B$subcelltypes),]
names(subs_B)

#### Check: ligand-receptor partner expression level for Endocrine ligand
test_E<-ie.interaction.2.ee[ie.interaction.2.ee$subcelltype1=="EnteroEndocrine_Cell"&ie.interaction.2.ee$age=="old"&ie.interaction.2.ee$compartment=="C_D",]
subs_E<-test_E[colSums(!is.na(test_E)) > 0]
View(subs_E)
test_E<-ie.interaction.2[ie.interaction.2$subcelltype1=="EnteroEndocrine_Cell"&ie.interaction.2$age=="old"&ie.interaction.2$compartment=="C_D",]
subs_E<-test_E[colSums(!is.na(test_E)) > 0]
colnames(subs_E)
test_E<-ie.interaction.2[ie.interaction.2$subcelltype1=="EnteroEndocrine_Cell"&ie.interaction.2$age=="old"&ie.interaction.2$compartment=="C_P",]
subs_E<-test_E[colSums(!is.na(test_E)) > 0]
colnames(subs_E)


#### removed Young_20w and jitter the dot on vln plot

exp_check<-function(IE.counts_norm.2,gene="ENSG00000019582"){
test_exp<-IE.counts_norm.2[IE.counts_norm.2$Gene==gene,]
row.names(test_exp)<-test_exp[,1]
test_exp<-test_exp[,-1]
test_exp_2<-as.data.frame(t(test_exp))
test_exp_2$id<-colnames(IE.counts_norm.2)[-1]
test_exp_2$compartment<-str_split_fixed(test_exp_2$id,"\\.",3)[,1]
test_exp_2$age<-str_split_fixed(test_exp_2$id,"\\.",3)[,2]
test_exp_2$celltype<-str_split_fixed(test_exp_2$id,"\\.",3)[,3]
test_exp_2[test_exp_2$age=="Young","age"]<-"young"
test_exp_2[test_exp_2$age=="Old","age"]<-"old"
test_exp_2[test_exp_2$compartment=="ProxCol","compartment"]<-"C_P"
test_exp_2[test_exp_2$compartment=="DistCol","compartment"]<-"C_D"
test_exp_2[test_exp_2$celltype %in% immune_c,"superC"]<-"Immune"
test_exp_2[test_exp_2$celltype %in% epithelial_c,"superC"]<-"Epithelial"
test_exp_2<-test_exp_2[test_exp_2$age!="Young_20w",]
test_exp_2$age<-factor(test_exp_2$age,levels = c("young","old"))
test_exp_2$compartment<-factor(test_exp_2$compartment,levels = c("Cecum","C_P","C_D"))
test_exp_2$celltype<-factor(test_exp_2$celltype,levels=c("Stem_cell","TA-cell","Colonocyte_Precursor","Cecum-enriched-Colonocyte","ProxCol-enriched-Colonocyte","DistCol-enriched-Colonocyte","EnteroEndocrine_Cell","Goblet_Precursor","Tuft_Cell","B-cell","Bcell-Early","CD4","CD4-Early","CD4-Cytotoxic","CD4-Naive","Bcell-Naive-Isg","CD4-TEM","CD8","Macrophage","Monocyte","Neutrophil","Plasma","Th17"))
p<-ggplot(test_exp_2[!is.na(test_exp_2$superC),], aes(x=celltype, y=get(gene), fill=age)) + geom_violin() +facet_grid(compartment~superC,scales='free_x')+theme(text = element_text(size=10),axis.text.x = element_text(angle=90, hjust=1))+ scale_fill_manual(values=c("grey", "grey0"))
q<-p+ geom_dotplot(binaxis='y', stackdir='center',dotsize=0.1,position="jitter")+ylab(paste("normalized_exp",gene,sep=" "))
print(q)
}

# CD74 ENSG00000019582
exp_check(IE.counts_norm.2,gene="ENSG00000019582")
# APP "ENSG00000142192"
exp_check(IE.counts_norm.2,gene="ENSG00000142192")
# MIF ENSG00000240972
exp_check(IE.counts_norm.2,gene="ENSG00000240972")
#ADRB2 ENSG00000169252
exp_check(IE.counts_norm.2,gene="ENSG00000169252")
#VEGFB ENSG00000173511
exp_check(IE.counts_norm.2,gene="ENSG00000173511")
#PYY ENSG00000131096
exp_check(IE.counts_norm.2,gene="ENSG00000131096")
#DPP4 ENSG00000197635
exp_check(IE.counts_norm.2,gene="ENSG00000197635")
# GCG ENSG00000115263
exp_check(IE.counts_norm.2,gene="ENSG00000115263")
#IDE ENSG00000119912
exp_check(IE.counts_norm.2,gene="ENSG00000119912")
# NPY1R ENSG00000164128

save.image(file='./data_in/CellPhoneDB.RData')

#### FN messages Figure3
#### interaction between epithelial-immune with epithelial
setwd("~/Desktop/Pro_Mine/Paper_Atlas/Immune")
load("./data_in/CellPhoneDB.RData")
Epithelial<-c("ProxCol-enriched-Colonocyte","Tuft_Cell","Stem_cell","EnteroEndocrine_Cell","TA-cell","Goblet_Precursor","Colonocyte_Precursor","Cecum-enriched-Colonocyte","DistCol-enriched-Colonocyte")

#### check EnteroEndocrine_Cell in C_P
test_E<-ie.interaction.2[ie.interaction.2$subcelltype1=="EnteroEndocrine_Cell"&ie.interaction.2$compartment=="C_P"&ie.interaction.2$subcelltype2 %in% Epithelial,]
subs_E<-test_E[colSums(!is.na(test_E)) > 0]
colnames(subs_E)
subs_E[is.na(subs_E)]<-0
subs_Ex<-as.data.frame(colSums(subs_E[1:(length(subs_E)-9)]))
subs_Ex$rn<-row.names(subs_Ex)
View(subs_Ex[order(subs_Ex[,1]),])
subs_Ex[order(subs_Ex[,1]),"rn"]
#### check EnteroEndocrine_Cell in C_D
test_E<-ie.interaction.2[ie.interaction.2$subcelltype1=="EnteroEndocrine_Cell"&ie.interaction.2$compartment=="C_D"&ie.interaction.2$subcelltype2 %in% Epithelial,]
subs_E<-test_E[colSums(!is.na(test_E)) > 0]
colnames(subs_E)
subs_E[is.na(subs_E)]<-0
subs_Ex<-as.data.frame(colSums(subs_E[1:(length(subs_E)-9)]))
subs_Ex$rn<-row.names(subs_Ex)
View(subs_Ex[order(subs_Ex[,1]),])
subs_Ex[order(subs_Ex[,1]),"rn"]
unique(ie.interaction.2$subcelltype1)
#### check Bcell-Early in C_P
test_E<-ie.interaction.2[ie.interaction.2$subcelltype1=="Bcell-Early"&ie.interaction.2$compartment=="C_P"&ie.interaction.2$subcelltype2 %in% Epithelial,]
subs_E<-test_E[colSums(!is.na(test_E)) > 0]
colnames(subs_E)
subs_E[is.na(subs_E)]<-0
subs_Ex<-as.data.frame(colSums(subs_E[1:(length(subs_E)-9)]))
subs_Ex$rn<-row.names(subs_Ex)
View(subs_Ex[order(subs_Ex[,1]),])
subs_Ex[order(subs_Ex[,1]),"rn"]
#### check Bcell-Early in C_D
test_E<-ie.interaction.2[ie.interaction.2$subcelltype1=="Bcell-Early"&ie.interaction.2$compartment=="C_D"&ie.interaction.2$subcelltype2 %in% Epithelial,]
subs_E<-test_E[colSums(!is.na(test_E)) > 0]
colnames(subs_E)
subs_E[is.na(subs_E)]<-0
subs_Ex<-as.data.frame(colSums(subs_E[1:(length(subs_E)-9)]))
subs_Ex$rn<-row.names(subs_Ex)
View(subs_Ex[order(subs_Ex[,1]),])
subs_Ex[order(subs_Ex[,1]),"rn"]
subs_Ex[order(subs_Ex[,1]),"rn"]

#### check gene expression in heatmap by group
library(stringr)
exp_check(IE.counts_norm.2,gene="ENSG00000019582")
#### genes check: GCG PYY IDE DPP4 NPY1R CD74 MIF APP
#### in EE ProxCol
genes_check<-c("ENSG00000115263","ENSG00000131096","ENSG00000119912","ENSG00000197635","ENSG00000164128","ENSG00000019582","ENSG00000240972","ENSG00000142192")
test_exp<-IE.counts_norm.2[IE.counts_norm.2$Gene %in% genes_check,]
row.names(test_exp)<-test_exp[,1]
test_exp<-test_exp[,-1]
test_exp_2<-as.data.frame(t(test_exp))
test_exp_2$id<-colnames(IE.counts_norm.2)[-1]
mean_group<-aggregate(. ~ id, test_exp_2, mean)
mean_group$compartment<-str_split_fixed(mean_group$id,"\\.",3)[,1]
mean_group$age<-str_split_fixed(mean_group$id,"\\.",3)[,2]
mean_group$celltype<-str_split_fixed(mean_group$id,"\\.",3)[,3]
mean_group[mean_group$age=="Young","age"]<-"young"
mean_group[mean_group$age=="Old","age"]<-"old"
mean_group[mean_group$compartment=="ProxCol","compartment"]<-"C_P"
mean_group[mean_group$compartment=="DistCol","compartment"]<-"C_D"
mean_group[mean_group$celltype %in% immune_c,"superC"]<-"Immune"
mean_group[mean_group$celltype %in% epithelial_c,"superC"]<-"Epithelial"
mean_group<-mean_group[mean_group$age!="Young_20w",]
mean_group$age<-factor(mean_group$age,levels = c("young","old"))
mean_group$compartment<-factor(mean_group$compartment,levels = c("Cecum","C_P","C_D"))
mean_group$celltype<-factor(mean_group$celltype,levels=c("Stem_cell","TA-cell","Colonocyte_Precursor","Cecum-enriched-Colonocyte","ProxCol-enriched-Colonocyte","DistCol-enriched-Colonocyte","EnteroEndocrine_Cell","Goblet_Precursor","Tuft_Cell","B-cell","Bcell-Early","CD4","CD4-Early","CD4-Cytotoxic","CD4-Naive","Bcell-Naive-Isg","CD4-TEM","CD8","Macrophage","Monocyte","Neutrophil","Plasma","Th17"))

EpithelialB<-c(Epithelial,"Bcell-Early")
mean_group$celltype<-str_split_fixed(mean_group$id,"\\.",3)[,3]
mean_group$celltype<-factor(mean_group$celltype,levels=c("Stem_cell","TA-cell","Colonocyte_Precursor","Cecum-enriched-Colonocyte","ProxCol-enriched-Colonocyte","DistCol-enriched-Colonocyte","EnteroEndocrine_Cell","Goblet_Precursor","Tuft_Cell","B-cell","Bcell-Early","CD4","CD4-Early","CD4-Cytotoxic","CD4-Naive","Bcell-Naive-Isg","CD4-TEM","CD8","Macrophage","Monocyte","Neutrophil","Plasma","Th17"))
mean_group$celltype_age<-paste(mean_group$celltype,mean_group$age,sep=".")
mean_group_cp_EE<-mean_group[mean_group$compartment=="C_P" & (mean_group$celltype %in% Epithelial),]
mean_group_cp_EE$celltype_age<-factor(mean_group_cp_EE$celltype_age,levels=c("Stem_cell.young","Stem_cell.old","TA-cell.young","TA-cell.old","Colonocyte_Precursor.young","Colonocyte_Precursor.old","Cecum-enriched-Colonocyte.young","Cecum-enriched-Colonocyte.old","ProxCol-enriched-Colonocyte.young","ProxCol-enriched-Colonocyte.old","DistCol-enriched-Colonocyte.young","DistCol-enriched-Colonocyte.old","EnteroEndocrine_Cell.young","EnteroEndocrine_Cell.old","Goblet_Precursor.young","Goblet_Precursor.old","Tuft_Cell.young","Tuft_Cell.old"))
library(pheatmap)
row.names(mean_group_cp_EE)<-mean_group_cp_EE$celltype_age
mean_group_cp_EE_t<-as.data.frame(t(mean_group_cp_EE[,c(2:9)]))
pheatmap(mean_group_cp_EE_t[c("ENSG00000115263","ENSG00000131096","ENSG00000119912","ENSG00000197635","ENSG00000164128"),c("Stem_cell.young","Stem_cell.old","TA-cell.young","TA-cell.old","Colonocyte_Precursor.young","Colonocyte_Precursor.old","Cecum-enriched-Colonocyte.young","Cecum-enriched-Colonocyte.old","ProxCol-enriched-Colonocyte.young","ProxCol-enriched-Colonocyte.old","DistCol-enriched-Colonocyte.young","DistCol-enriched-Colonocyte.old","EnteroEndocrine_Cell.young","EnteroEndocrine_Cell.old","Goblet_Precursor.young","Goblet_Precursor.old","Tuft_Cell.young","Tuft_Cell.old")],cluster_rows = F,cluster_cols = F,scale="row")
#### EE in distCol: /Users/jlu/Desktop/Pro_Mine/Paper_Atlas/Immune/data_out/all_0.1_out_var_HM/out_expression_check/C-D_EE.pdf
mean_group_cd_EE<-mean_group[mean_group$compartment=="C_D" & (mean_group$celltype %in% Epithelial),]
row.names(mean_group_cd_EE)<-mean_group_cd_EE$celltype_age
mean_group_cd_EE_t<-as.data.frame(t(mean_group_cd_EE[,c(2:9)]))
#### pheatmap(mean_group_cd_EE_t[c("ENSG00000115263","ENSG00000131096","ENSG00000119912","ENSG00000197635","ENSG00000164128"),c("Stem_cell.young","Stem_cell.old","TA-cell.young","TA-cell.old","Colonocyte_Precursor.young","Colonocyte_Precursor.old","Cecum-enriched-Colonocyte.young","Cecum-enriched-Colonocyte.old","ProxCol-enriched-Colonocyte.young","ProxCol-enriched-Colonocyte.old","DistCol-enriched-Colonocyte.young","DistCol-enriched-Colonocyte.old","EnteroEndocrine_Cell.young","EnteroEndocrine_Cell.old","Goblet_Precursor.young","Goblet_Precursor.old","Tuft_Cell.young","Tuft_Cell.old")],cluster_rows = F,cluster_cols = F,scale="row")
mean_group_cd_EE_t[c("ENSG00000115263","ENSG00000131096","ENSG00000119912","ENSG00000197635","ENSG00000164128"),]
pheatmap(mean_group_cd_EE_t[c("ENSG00000115263","ENSG00000131096","ENSG00000119912","ENSG00000197635","ENSG00000164128"),c("Stem_cell.young","Stem_cell.old","TA-cell.young","TA-cell.old","Colonocyte_Precursor.young","Colonocyte_Precursor.old","Cecum-enriched-Colonocyte.young","Cecum-enriched-Colonocyte.old","ProxCol-enriched-Colonocyte.young","ProxCol-enriched-Colonocyte.old","EnteroEndocrine_Cell.young","EnteroEndocrine_Cell.old","Goblet_Precursor.young","Goblet_Precursor.old","Tuft_Cell.young","Tuft_Cell.old")],cluster_rows = F,cluster_cols = F,scale="row")
#### EE in Cecum: /Users/jlu/Desktop/Pro_Mine/Paper_Atlas/Immune/data_out/all_0.1_out_var_HM/out_expression_check/Cecum_EE.pdf
mean_group_cecum_EE<-mean_group[mean_group$compartment=="Cecum" & (mean_group$celltype %in% Epithelial),]
row.names(mean_group_cecum_EE)<-mean_group_cecum_EE$celltype_age
mean_group_cecum_EE_t<-as.data.frame(t(mean_group_cecum_EE[,c(2:9)]))
pheatmap(mean_group_cecum_EE_t[c("ENSG00000115263","ENSG00000131096","ENSG00000119912","ENSG00000197635","ENSG00000164128"),c("Stem_cell.young","Stem_cell.old","TA-cell.young","TA-cell.old","Colonocyte_Precursor.young","Colonocyte_Precursor.old","Cecum-enriched-Colonocyte.young","Cecum-enriched-Colonocyte.old","ProxCol-enriched-Colonocyte.young","ProxCol-enriched-Colonocyte.old","DistCol-enriched-Colonocyte.young","DistCol-enriched-Colonocyte.old","EnteroEndocrine_Cell.young","EnteroEndocrine_Cell.old","Goblet_Precursor.young","Goblet_Precursor.old","Tuft_Cell.young","Tuft_Cell.old")],cluster_rows = F,cluster_cols = F,scale="row")
#### EB in Cecum: /Users/jlu/Desktop/Pro_Mine/Paper_Atlas/Immune/data_out/all_0.1_out_var_HM/out_expression_check/Cecum_EB.pdf
mean_group_cecum_EB<-mean_group[mean_group$compartment=="Cecum" & (mean_group$celltype %in% EpithelialB),]
row.names(mean_group_cecum_EB)<-mean_group_cecum_EB$celltype_age
mean_group_cecum_EB_t<-as.data.frame(t(mean_group_cecum_EB[,c(2:9)]))
pheatmap(mean_group_cecum_EB_t[c("ENSG00000019582","ENSG00000240972","ENSG00000142192"),c("Stem_cell.young","Stem_cell.old","TA-cell.young","TA-cell.old","Colonocyte_Precursor.young","Colonocyte_Precursor.old","Cecum-enriched-Colonocyte.young","Cecum-enriched-Colonocyte.old","ProxCol-enriched-Colonocyte.young","ProxCol-enriched-Colonocyte.old","DistCol-enriched-Colonocyte.young","DistCol-enriched-Colonocyte.old","EnteroEndocrine_Cell.young","EnteroEndocrine_Cell.old","Goblet_Precursor.young","Goblet_Precursor.old","Tuft_Cell.young","Tuft_Cell.old","Bcell-Early.young","Bcell-Early.old")],cluster_rows = F,cluster_cols = F,scale="row")

#### EB in proxCol: /Users/jlu/Desktop/Pro_Mine/Paper_Atlas/Immune/data_out/all_0.1_out_var_HM/out_expression_check/C-P_EB.pdf
mean_group_cp_EB<-mean_group[mean_group$compartment=="C_P" & (mean_group$celltype %in% EpithelialB),]
row.names(mean_group_cp_EB)<-mean_group_cp_EB$celltype_age
mean_group_cp_EB_t<-as.data.frame(t(mean_group_cp_EB[,c(2:9)]))
pheatmap(mean_group_cp_EB_t[c("ENSG00000019582","ENSG00000240972","ENSG00000142192"),c("Stem_cell.young","Stem_cell.old","TA-cell.young","TA-cell.old","Colonocyte_Precursor.young","Colonocyte_Precursor.old","Cecum-enriched-Colonocyte.young","Cecum-enriched-Colonocyte.old","ProxCol-enriched-Colonocyte.young","ProxCol-enriched-Colonocyte.old","DistCol-enriched-Colonocyte.young","DistCol-enriched-Colonocyte.old","EnteroEndocrine_Cell.young","EnteroEndocrine_Cell.old","Goblet_Precursor.young","Goblet_Precursor.old","Tuft_Cell.young","Tuft_Cell.old","Bcell-Early.young","Bcell-Early.old")],cluster_rows = F,cluster_cols = F,scale="row")
#### EB in distCol: /Users/jlu/Desktop/Pro_Mine/Paper_Atlas/Immune/data_out/all_0.1_out_var_HM/out_expression_check/C-D_EB.pdf
mean_group_cd_EB<-mean_group[mean_group$compartment=="C_D" & (mean_group$celltype %in% EpithelialB),]
row.names(mean_group_cd_EB)<-mean_group_cd_EB$celltype_age
mean_group_cd_EB_t<-as.data.frame(t(mean_group_cd_EB[,c(2:9)]))
pheatmap(mean_group_cd_EB_t[c("ENSG00000019582","ENSG00000240972","ENSG00000142192"),c("Stem_cell.young","Stem_cell.old","TA-cell.young","TA-cell.old","Colonocyte_Precursor.young","Colonocyte_Precursor.old","Cecum-enriched-Colonocyte.young","Cecum-enriched-Colonocyte.old","ProxCol-enriched-Colonocyte.young","ProxCol-enriched-Colonocyte.old","DistCol-enriched-Colonocyte.young","DistCol-enriched-Colonocyte.old","EnteroEndocrine_Cell.young","EnteroEndocrine_Cell.old","Goblet_Precursor.young","Goblet_Precursor.old","Tuft_Cell.young","Tuft_Cell.old","Bcell-Early.young","Bcell-Early.old")],cluster_rows = F,cluster_cols = F,scale="row")
pheatmap(mean_group_cd_EB_t[c("ENSG00000019582","ENSG00000240972","ENSG00000142192"),c("Stem_cell.young","Stem_cell.old","TA-cell.young","TA-cell.old","Colonocyte_Precursor.young","Colonocyte_Precursor.old","Cecum-enriched-Colonocyte.young","Cecum-enriched-Colonocyte.old","ProxCol-enriched-Colonocyte.young","ProxCol-enriched-Colonocyte.old","EnteroEndocrine_Cell.young","EnteroEndocrine_Cell.old","Goblet_Precursor.young","Goblet_Precursor.old","Tuft_Cell.young","Tuft_Cell.old","Bcell-Early.young","Bcell-Early.old")],cluster_rows = F,cluster_cols = F,scale="row")

#### Bar chart of up-stream cytokines predicted by Epithelial
up_cytokine_cecum<-read.csv("/Users/jlu/Desktop/Pro_Mine/Paper_Atlas/Epithelial/data_in/Upstream_Analysis/Cecum_YvO_Upstream_Analysis.xls.csv.cytokine.csv",skip = 2,head=T,sep=";",stringsAsFactors = F)
up_cytokine_cp<-read.csv("/Users/jlu/Desktop/Pro_Mine/Paper_Atlas/Epithelial/data_in/Upstream_Analysis/ProxCol_YvO_Upstream_Analysis.xls.csv.cytokine.csv",skip = 2,head=T,sep=";",stringsAsFactors = F)
up_cytokine_cd<-read.csv("/Users/jlu/Desktop/Pro_Mine/Paper_Atlas/Epithelial/data_in/Upstream_Analysis/DistCol_YvO_Upstream_Analysis.xls.csv.cytokine.csv",skip = 2,head=T,sep=";",stringsAsFactors = F)
up_cytokine_cecum$pos<-"Cecum"
up_cytokine_cp$pos<-"C_P"
up_cytokine_cd$pos<-"C_D"
up_cytokines<-rbind(up_cytokine_cecum,up_cytokine_cp,up_cytokine_cd)
up_cytokines$p_mlog10<-(-log10(up_cytokines$p.value.of.overlap))
-log10(0.05)=1.3


#### facet_grid coudl adjust size automatically by parametre of "space" but facet_wrap could not

ggplot(data=up_cytokines, aes(x=Upstream.Regulator, y=p_mlog10)) +
geom_bar(stat="identity")+ coord_flip()+theme_classic()+facet_grid(pos ~ ., scales = "free", space='free')


up_cytokines$pos<-factor(up_cytokines$pos,levels = c("Cecum","C_P","C_D"))
#### order by value
up_cytokines_copy<-up_cytokines
up_cytokines$Upstream.Regulator_pos<-paste(up_cytokines$pos,up_cytokines$Upstream.Regulator,sep="#")
up_cytokines$Upstream.Regulator_pos<-factor(up_cytokines$Upstream.Regulator_pos, levels = unique(up_cytokines$Upstream.Regulator_pos[order(up_cytokines$p_mlog10)]))


ggplot(data=up_cytokines, aes(x=Upstream.Regulator_pos, y=p_mlog10)) +
geom_bar(stat="identity")+ coord_flip()+theme_classic()+facet_grid(pos ~ ., scales = "free", space='free') + geom_hline(yintercept = 1.3, linetype="dashed",size=0.5)

#### remove IL22 in Cecum which is not significant
#### FN_msg4_Up_cytokines_barcharts.pdf
up_cytokines<-up_cytokines[-3,]
ggplot(data=up_cytokines, aes(x=Upstream.Regulator_pos, y=p_mlog10)) +
geom_bar(stat="identity")+ coord_flip()+theme_classic()+facet_grid(pos ~ ., scales = "free", space='free') + geom_hline(yintercept = 1.3, linetype="dashed",size=0.5)

#### input for venn plot: /Users/jlu/Desktop/Pro_Mine/Paper_Atlas/Epithelial/data_in/Upstream_Analysis/compartment_venn.svg
paste(as.character(up_cytokine_cd$Upstream.Regulator), sep="' '", collapse=", ")
paste(as.character(up_cytokine_cp$Upstream.Regulator), sep="' '", collapse=", ")
dput(as.character(up_cytokine_cp$Upstream.Regulator))
dput(as.character(up_cytokine_cd$Upstream.Regulator))
up_cytokine_cd[,"Upstream.Regulator",drop=F]
print(up_cytokine_cd[,"Upstream.Regulator",drop=F],row.names=F,colnames=F)
print(up_cytokine_cp[,"Upstream.Regulator",drop=F],row.names=F,colnames=F)
print(up_cytokine_cecum[,"Upstream.Regulator",drop=F],row.names=F,colnames=F)


