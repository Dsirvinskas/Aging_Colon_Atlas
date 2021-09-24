library(pheatmap)
library(stringr)
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(args[1])
print(args[2])

dir.create("./out")
heatmaps_plot = function(meta_file=args[1], pvalues_file=args[2], count_filename="./out/heatmap_count.pdf", log_filename="./out/heatmap_log_count.pdf", count_network_filename="./out/count_network.txt", interaction_count_filename="./out/interaction_count.txt", count_network_separator="\t", interaction_count_separator="\t", show_rownames = T, show_colnames = T,
                         scale="none", cluster_cols = T,border_color='white', cluster_rows = T, fontsize_row=1,
                         fontsize_col = 1, main = '',treeheight_row=0, family='Arial', treeheight_col = 0,
                         col1 = "dodgerblue4", col2 = 'peachpuff', col3 = 'deeppink4', meta_sep='\t', pvalues_sep='\t', pvalue=0.05){
  #######   Network

  meta = read.csv(meta_file, comment.char = '', sep=meta_sep)

  all_intr = read.table(pvalues_file, header=T, stringsAsFactors = F, sep=pvalues_sep, comment.char = '', check.names = F)
  intr_pairs = all_intr$interacting_pair
  all_intr = all_intr[,-c(1:11)]


  split_sep = '\\|'
  join_sep = '|'

  pairs1_all = unique(meta[,2])

  pairs1 = c()
  for (i in 1:length(pairs1_all))
    for (j in 1:length(pairs1_all))
        pairs1 = c(pairs1,paste(pairs1_all[i],pairs1_all[j],sep=join_sep))

  all_count = matrix(ncol=3)
  colnames(all_count) = c('SOURCE','TARGET','count')
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]

    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]

    if(p1!=p2)
      count1 = length(unique(n2))
    else
      count1 = length(unique(n2))

    new_count = c(p1,p2,count1)
    names(new_count) = c('SOURCE','TARGET','count')
    all_count = rbind(all_count, new_count)
  }

  all_count = all_count[-1,]
  write.table(all_count, count_network_filename, sep=count_network_separator, quote=F, row.names = F)

  #######   count interactions

  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]

    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]

    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    if(p1!=p2)
      count1 = c(count1,length(unique(n2)))
    else
      count1 = c(count1,length(unique(n2)))

  }
  
  if (any(count1)>0)
  {
    count_matrix = matrix(count1, nrow=length(unique(meta[,2])), ncol=length(unique(meta[,2])))
    rownames(count_matrix)= unique(meta[,2])
    colnames(count_matrix)= unique(meta[,2])

    all_sum = rowSums(count_matrix)
    all_sum = cbind(names(all_sum), all_sum)
    #write.table(all_sum, file=interaction_count_filename, quote=F, sep=count_network_separator, row.names=F)
    write.table(count_matrix, file=interaction_count_filename, quote=F, sep=count_network_separator, row.names=T)
    col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( 1000 )

   # pheatmap(count_matrix, show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
   #          border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
   #          main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = count_filename)
    
    df<-data.frame(coln=colnames(count_matrix),leve=str_split_fixed(colnames(count_matrix),"\\.",3)[,3])
    df$leve<-as.character(df$leve)
    df$coln<-as.character(df$coln)
    orders_ie<-c("Stem_cell","TA-cell","Colonocyte_Precursor","Cecum-enriched-Colonocyte","ProxCol-enriched-Colonocyte","DistCol-enriched-Colonocyte","EnteroEndocrine_Cell","Goblet_Precursor","Goblet","Tuft_Cell","Bcell-Early","B-cell","Bcell-Naive-Isg","Plasma","CD4-Early","CD4-Naive","CD4","CD4-Cytotoxic","CD4-TEM","CD8","Th17","Monocyte","Macrophage","Neutrophil")
    dfo<-df[match(orders_ie,df$leve),] 
    dfo<-na.omit(dfo)
    count_matrix<-count_matrix[dfo$coln,dfo$coln]
    colnames(count_matrix)<-str_split_fixed(colnames(count_matrix),"\\.",3)[,3]
    row.names(count_matrix)<-str_split_fixed(row.names(count_matrix),"\\.",3)[,3]
    breaksList<-seq(0, 15, by = 1)
    col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( length(breaksList) )    
    pheatmap(count_matrix, show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
             border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
             main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap,breaks = breaksList, treeheight_col = treeheight_col, filename = count_filename)

   # pheatmap(log(count_matrix+1), show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
   #          border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
   #          main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = log_filename)
  } else {
    stop("There are no significant results using p-value of: ", pvalue, call.=FALSE)
  }
}

## execute
#heatmaps_plot(meta_file=args[1], pvalues_file=args[2],fontsize_row=2,fontsize_col=2)
heatmaps_plot(meta_file=args[1], pvalues_file=args[2],count_filename=args[3],log_filename=args[4],interaction_count_filename=args[5],count_network_filename=args[6],fontsize_row=10,fontsize_col=10,cluster_cols=F,cluster_rows=F)
