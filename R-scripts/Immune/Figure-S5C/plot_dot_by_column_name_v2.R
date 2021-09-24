library(ggplot2)
library(stringr)
library(forcats)
args = commandArgs(trailingOnly=TRUE)
print(args[1])
print(args[2])
dir.create("./out")
dot_plot = function(selected_rows = NULL,
                    selected_columns = NULL,
                    filename = paste("./out/plot_",args[1],".pdf",sep=''),
                    width = 8,
                    height = 10,
                    means_path = './means.txt',
                    pvalues_path = './pvalues.txt',
                    means_separator = '\t',
                    pvalues_separator = '\t',
                    output_extension = '.pdf'
){
  all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
  all_means = read.table(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)

  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]

  if(is.na(selected_rows)){
    selected_rows = intr_pairs
    print("AAAAAA")
  } else {print("BBB"); selected_rows<-read.table(selected_rows,header=F,stringsAsFactors = F,check.names=F,,sep="\t");selected_rows<-as.character(selected_rows[,1])}

  if(is.na(selected_columns)){
    print("CCC")
    selected_columns = colnames(all_pval)
  } else {print("DDD");selected_columns<-read.table(selected_columns,header=F,stringsAsFactors = F,check.names=F,sep="\t");selected_columns<-as.character(selected_columns[,1]);print(head(selected_columns));print(!grepl("\\|",selected_columns[1]))}
  if (!grepl("\\|",selected_columns[1])) {
    print("Combination")
    pairs1 = c()
    for (i in 1:length(selected_columns))
      for (j in 1:length(selected_columns))
          pairs1 = c(pairs1,paste(selected_columns[i],selected_columns[j],sep="|"))
    selected_columns<-pairs1
  }
  
  sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]

  df_names = expand.grid(selected_rows, selected_columns)
  pval = unlist(sel_pval)
  pval[pval==0] = 0.0009
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  pr[pr==0] = 1
  plot.data = cbind(plot.data,log2(pr))
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')


##################################################################
  
  plot.data$one<-str_split_fixed(plot.data$clusters,"\\|",2)[,1]
  plot.data$two<-str_split_fixed(plot.data$clusters,"\\|",2)[,2]  
  plot.data$one_ct<-str_split_fixed(plot.data$one,"\\.",3)[,3]
  plot.data$two_ct<-str_split_fixed(plot.data$two,"\\.",3)[,3]
  plot.data$age<-str_split_fixed(plot.data$one,"\\.",3)[,2]
  plot.data$pos<-str_split_fixed(plot.data$one,"\\.",3)[,1]
  plot.data$both_ct<-paste(plot.data$one_ct,plot.data$two_ct,sep="|")
  plot.data$pair_age<-paste(plot.data$pair,plot.data$age,sep=".")
  plot.data$both_ct<-factor(plot.data$both_ct,levels=unique(plot.data$both_ct))
  plot.data$pair_age<-factor(plot.data$pair_age,levels=unique(plot.data$pair_age))
  plot.data$age<-factor(plot.data$age,levels=unique(plot.data$age))
  plot.data$pair<-factor(plot.data$pair,levels=unique(plot.data$pair))

  my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)

  saveRDS(plot.data,paste(args[1],".plot.data.RDS"))
  
  ggplot(plot.data,aes(x=both_ct,y=fct_rev(age))) +
    geom_point(aes(size=-log10(pvalue),color=mean)) +
    scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text=element_text(size=5, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1,size=10),
          axis.text.y = element_text(size=10, colour = "black"),
          axis.title=element_blank(),
          strip.text.y = element_text(angle=0),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black")) + facet_grid(fct_rev(plot.data$pair)~.)

#  ggplot(plot.data,aes(x=both_ct,y=pair_age)) +
#  geom_point(aes(size=-log10(pvalue),color=mean)) +
#  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
#  theme_bw() +
#  theme(panel.grid.minor = element_blank(),
#        panel.grid.major = element_blank(),
#        axis.text=element_text(size=5, colour = "black"),
#        axis.text.x = element_text(angle = 90, hjust = 1,size=10),
#        axis.text.y = element_text(size=10, colour = "black"),
#        axis.title=element_blank(),
#        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

  if (output_extension == '.pdf') {
      ggsave(filename, width = width, height = height, device = cairo_pdf, limitsize=F)
  }
  else {
      ggsave(filename, width = width, height = height, limitsize=F)
  }
}
# execute

dot_plot(selected_rows=args[1],selected_columns=args[2], width = 8, height = 14)

