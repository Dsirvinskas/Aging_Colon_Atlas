# Rscript UpstreamAnalysis B-cell_age.IPA.UpstreamAnalysis.txt all 20
# IPA pathway analysis 
# load packages
options(stringsAsFactors = FALSE)

library(ggplot2)
library(treemap)
library(pheatmap)
library(reshape2)
library(stringr)


# read input
args = commandArgs(trailingOnly=TRUE)
IPA_input<-read.table(args[1],sep="\t",skip=2,stringsAsFactors = FALSE,header=TRUE,na.strings="NaN",quote = "")
colnames(IPA_input)[1]<-"ur"
colnames(IPA_input)[3]<-"type"
colnames(IPA_input)[5]<-"zscore"
colnames(IPA_input)[7]<-"pvalue"
colnames(IPA_input)[8]<-"genes_down"
IPA_input_Sig<-IPA_input[IPA_input$pvalue<0.05,]  



if(args[2]!="all"){IPA_input_Sig<-IPA_input_Sig[IPA_input_Sig$type==args[2],]}
IPA_input_Sig<-head(IPA_input_Sig,n=args[3])
IPA_input_Sig[is.na(IPA_input_Sig$zscore),"zscore"]<-0
IPA_input_Sig$genes_down_num<-lengths(strsplit(IPA_input_Sig$genes_down,","))
IPA_input_Sig$ur<-factor(IPA_input_Sig$ur,levels=unique(IPA_input_Sig$ur))
ggplot(IPA_input_Sig,aes(x=-1*log10(pvalue),y=ur)) + 
geom_point(aes(size=genes_down_num,color=zscore))+
scale_colour_gradient(low="green",high="red")+
labs(
color=expression(zscore),
size="num of target genes",
x="-log10(pvalue)" # Fold.Enrichment for KEGG | Gene Ratio for BP
# y="Pathway name",
# title="Pathway enrichment")
)+
theme_classic()+
theme(
axis.text.y = element_text(size = rel(2)),
#axis.title.x = element_text(size=rel(1.3)),
axis.title.y = element_blank()#,legend.position="bottom"
)+scale_y_discrete(limits = rev(levels(IPA_input_Sig$ur)))+scale_size_continuous(range=c(10,20)) # + xlim(0,20)
#ggsave(paste(name,"tiff",sep='.'))
ggsave(paste("URplot","pdf",sep='.'),width=15,height=10)
