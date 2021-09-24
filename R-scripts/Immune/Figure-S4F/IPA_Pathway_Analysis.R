
# IPA pathway analysis 
# load packages
options(stringsAsFactors = FALSE)

library(ggplot2)
library(treemap)
library(pheatmap)
library(reshape2)


# read input
args = commandArgs(trailingOnly=TRUE)
IPA_input<-read.table(args[1],sep="\t",skip=2,stringsAsFactors = FALSE,header=TRUE,na.strings="NaN",quote = "")
IPA_input<-IPA_input[!substring(IPA_input[,1],1,1)=="#",]
colnames(IPA_input)[1]<-"icp"
colnames(IPA_input)[2]<-"logp"
IPA_input_Sig<-IPA_input[IPA_input$logp>-log10(0.05),]  

# IPA supercategory correlation method
ipa_catagory<-data.frame()
for (i in 1:nrow(IPA_input_Sig)) {
  tmp<-data.frame("genes"=strsplit(IPA_input_Sig$Molecules[i], ",")[[1]])
  tmp$icp<-IPA_input_Sig$icp[i]
  ipa_catagory<-rbind(ipa_catagory,tmp)
}
ipa_catagory<-unique(ipa_catagory)
ipa_catagory$value<-1
ipa_catagory_wide<-dcast(ipa_catagory,icp~genes,value.var = 'value')
ipa_catagory_wide[is.na(ipa_catagory_wide)]<-0
row.names(ipa_catagory_wide)<-ipa_catagory_wide$icp
ipa_catagory_wide<-ipa_catagory_wide[,-1]
res <- cor(t(ipa_catagory_wide),method = c("pearson"))
p1<-pheatmap(res,fontsize=5,treeheight_row = 0, treeheight_col = 0)

# IPA supercategory ipa method
IPA<-read.table("/Users/jlu/Desktop/Pro_Mine/IpaC/db_ipa//IPAPathways_Final.txt",stringsAsFactors = FALSE,header=TRUE)
IPAm<-as.matrix(IPA)
level<-3
upstr<-function(x){p<-which(IPAm==x,arr.ind=TRUE);return(unique(IPAm[p[,1],level]))}
IPAPath<-function(p){if(all(! is.na(p))){return(p[length(p)])}else{return(p[min(which(is.na(p)))-1])}}
IPN<-apply(IPA,1,IPAPath)
IPAres<-ipa_catagory_wide[row.names(ipa_catagory_wide) %in% IPN,]
IPAres$genesNum<-rowSums(IPAres)
for(i in 1:nrow(IPAres)){IPAres$Upstream[i]<-upstr(row.names(IPAres)[i])}
IPAres$Pathway<-row.names(IPAres)
p2<-treemap(
IPAres,
fontcolor.labels=c("white","black"),
fontsize.labels=c(25,20),
position.legend="none",
index = c("Upstream","Pathway"),
vSize = "genesNum",
title = "title",
inflate.labels = FALSE,
lowerbound.cex.labels = 0,
bg.labels = "transparent" ,
overlap.labels=0.5,
)

# bar plot for significant ones
IPA_input_Sig<-IPA_input_Sig[order(IPA_input_Sig$logp),]
IPA_input_Sig$icp<-factor(IPA_input_Sig$icp,levels=IPA_input_Sig$icp)
myColor <- colorRampPalette(c("blue","white", "red"))(50)
p3<-ggplot(IPA_input_Sig,aes(x=icp,y=logp))+
    scale_fill_gradientn(colours=myColor)+
    geom_bar(stat = "identity")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("Pathway category")+
    ylab("-log10 p")+
    coord_flip()+
    geom_hline(yintercept =1.30103)+
    theme(text = element_text(size=10))
print(p3)



