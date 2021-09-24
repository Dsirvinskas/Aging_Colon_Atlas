require(ggplot2)
require(dplyr)
require(Matrix)
require(cowplot)
require(Seurat)


data = readMM(file= "/data/tmp/Dovydas/Single_Cell_RNAseq/2020_02_05_YvO_Single_Cell/Merged/outs/filtered_feature_bc_matrix/matrix.mtx.gz") 
cellID = read.csv(file="/data/tmp/Dovydas/Single_Cell_RNAseq/2020_02_05_YvO_Single_Cell/Merged/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", header = FALSE, sep="") 
genesID = read.csv(file="/data/tmp/Dovydas/Single_Cell_RNAseq/2020_02_05_YvO_Single_Cell/Merged/outs/filtered_feature_bc_matrix/features.tsv.gz", header = FALSE, sep="") 
colnames(data) = cellID$V1 
rownames(data) = genesID$V2 

## Separating genes from hashtags
hashtag = data[c((nrow(data)-11):nrow(data)),]
data = data[c(1:(nrow(data)-12)),]

## Checking and assigning hashtags to cells
hashtag = as.matrix(hashtag) ## Fill in the blank spaces in sparce matrix
hashtag = rbind(hashtag, apply(hashtag,2,sum)) ## Sum of reads for hashtag - will further utilize cells with more than 10
row.names(hashtag)[13] = "Sum"

## Calculating the fraction of reads coming from each hashtag in every cell
hashtag1 = hashtag
for (i in 1:12){
  hashtag1[i,] = hashtag[i,]/hashtag[13,]
}

## Removal of cells that have less than 10 reads for all hashtags combined
hashtag1 = hashtag1[,hashtag1[13,] >9]   

## Assigning hashtag to cells - if one hashtag is over 65% of the total hashtag reads, the cell is assigned this hashtag
hashtag2 = hashtag1[1:12,] > 0.65
hash = data.frame(CellID=colnames(hashtag2),hashtag=as.character(sapply(1:ncol(hashtag2),function(x) names(which(hashtag2[,x]==TRUE)))))
hash = hash[hash[,2] != "character(0)",]
hash = droplevels(hash)

## Replacing hashtag with correlating sample compartment
hash = cbind(hash, "Batch" = substr(hash[1:nrow(hash),1],18,18))


info = read.csv(file="hashtags.csv", header = TRUE)
info_age = read.csv(file="hashtags_age.csv", header = TRUE)

comp = sapply(1:nrow(hash), FUN= function(x){
  a = hash[x,3]
  b = as.character(hash[x,2])
  info[a,b]
})
hash = cbind(hash, "Compartment" = comp)

samp = sapply(1:nrow(hash), FUN= function(x){
  a = hash[x,3]
  b = as.character(hash[x,2])
  info_age[a,b]
})
hash = cbind(hash, "Sample" = samp)

save.image("StartingData.RData")
## Can start running Seurat code 