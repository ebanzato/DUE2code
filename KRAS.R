
#---- Libraries -----------------------------------------------------

library(graphite)
library(ggplot2)
library(gRbase)
library(graph)
library(igraph)
library(stringr)
library(scran)

# theme_set(theme_bw())

path = 'E:/Università/PHD/DATI/KRAS/data'



#---- DATA ---------------------------------------------------------------------

cell.annotation = read.csv(paste0(path,'/GSE137912_cell.annotation.csv.gz'))
      
## Upload data

## H0 ----
# H358
H358 = as.matrix(readMM(file = paste0(path,'/GSM4094251_H358.1.matrix.mtx.gz')))
colnames(H358) = read.table(paste0(path,'/GSM4094251_H358.1.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(H358) = read.table(paste0(path,'/GSM4094251_H358.1.genes.tsv.gz'), header=F)[,1]

colnames(H358) = str_sub(colnames(H358), end=-3)
keep = cell.annotation$X[cell.annotation$Hours==0 & cell.annotation$Line=='H358']
keep = str_sub(keep, start = 8)
H358 = H358[,keep]

# H2122
H2122 = as.matrix(readMM(file = paste0(path,'/GSM4094255_H2122.0.matrix.mtx.gz')))
colnames(H2122) = read.table(paste0(path,'/GSM4094255_H2122.0.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(H2122) = read.table(paste0(path,'/GSM4094255_H2122.0.genes.tsv.gz'), header=F)[,1]

colnames(H2122) = str_sub(colnames(H2122), end=-3)
keep = cell.annotation$X[cell.annotation$Hours==0 & cell.annotation$Line=='H2122']
keep = str_sub(keep, start = 9)
H2122 = H2122[,keep]

# SW1573
SW1573 = as.matrix(readMM(file = paste0(path,'/GSM4094259_SW1573.0.matrix.mtx.gz')))
colnames(SW1573) = read.table(paste0(path,'/GSM4094259_SW1573.0.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(SW1573) = read.table(paste0(path,'/GSM4094259_SW1573.0.genes.tsv.gz'), header=F)[,1]

colnames(SW1573) = str_sub(colnames(SW1573), end=-3)
keep = cell.annotation$X[cell.annotation$Hours==0 & cell.annotation$Line=='SW1573']
keep = str_sub(keep, start = 10)
SW1573 = SW1573[,keep]

## H4 ----
# H358
H358_4 = as.matrix(readMM(file = paste0(path,'/GSM4094252_H358.2.matrix.mtx.gz')))
colnames(H358_4) = read.table(paste0(path,'/GSM4094252_H358.2.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(H358_4) = read.table(paste0(path,'/GSM4094252_H358.2.genes.tsv.gz'), header=F)[,1]

colnames(H358_4) = str_sub(colnames(H358_4), end=-3)
keep = cell.annotation$X[cell.annotation$Hours==4 & cell.annotation$Line=='H358']
keep = str_sub(keep, start = 8)
H358_4 = H358_4[,keep]

# H2122
H2122_4 = as.matrix(readMM(file = paste0(path,'/GSM4094256_H2122.4.matrix.mtx.gz')))
colnames(H2122_4) = read.table(paste0(path,'/GSM4094256_H2122.4.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(H2122_4) = read.table(paste0(path,'/GSM4094256_H2122.4.genes.tsv.gz'), header=F)[,1]

colnames(H2122_4) = str_sub(colnames(H2122_4), end=-3)
keep = cell.annotation$X[cell.annotation$Hours==4 & cell.annotation$Line=='H2122']
keep = str_sub(keep, start = 9)
H2122_4 = H2122_4[,keep]

# SW1573
SW1573_4 = as.matrix(readMM(file = paste0(path,'/GSM4094260_SW1573.4.matrix.mtx.gz')))
colnames(SW1573_4) = read.table(paste0(path,'/GSM4094260_SW1573.4.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(SW1573_4) = read.table(paste0(path,'/GSM4094260_SW1573.4.genes.tsv.gz'), header=F)[,1]

colnames(SW1573_4) = str_sub(colnames(SW1573_4), end=-3)
keep = cell.annotation$X[cell.annotation$Hours==4 & cell.annotation$Line=='SW1573']
keep = str_sub(keep, start = 10)
SW1573_4 = SW1573_4[,keep]

## H24 ----
# H358
H358_24 = as.matrix(readMM(file = paste0(path,'/GSM4094253_H358.3.matrix.mtx.gz')))
colnames(H358_24) = read.table(paste0(path,'/GSM4094253_H358.3.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(H358_24) = read.table(paste0(path,'/GSM4094253_H358.3.genes.tsv.gz'), header=F)[,1]

colnames(H358_24) = str_sub(colnames(H358_24), end=-3)
keep = cell.annotation$X[cell.annotation$Hours==24 & cell.annotation$Line=='H358']
keep = str_sub(keep, start = 8)
H358_24 = H358_24[,keep]

# H2122
H2122_24 = as.matrix(readMM(file = paste0(path,'/GSM4094257_H2122.24.matrix.mtx.gz')))
colnames(H2122_24) = read.table(paste0(path,'/GSM4094257_H2122.24.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(H2122_24) = read.table(paste0(path,'/GSM4094257_H2122.24.genes.tsv.gz'), header=F)[,1]

colnames(H2122_24) = str_sub(colnames(H2122_24), end=-3)
keep = cell.annotation$X[cell.annotation$Hours==24 & cell.annotation$Line=='H2122']
keep = str_sub(keep, start = 9)
H2122_24 = H2122_24[,keep]
rm(H2122_24)

# SW1573
SW1573_24 = as.matrix(readMM(file = paste0(path,'/GSM4094261_SW1573.24.matrix.mtx.gz')))
colnames(SW1573_24) = read.table(paste0(path,'/GSM4094261_SW1573.24.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(SW1573_24) = read.table(paste0(path,'/GSM4094261_SW1573.24.genes.tsv.gz'), header=F)[,1]

colnames(SW1573_24) = str_sub(colnames(SW1573_24), end=-3)
keep = cell.annotation$X[cell.annotation$Hours==24 & cell.annotation$Line=='SW1573']
keep = str_sub(keep, start = 10)
SW1573_24 = SW1573_24[,keep]


## H72 ----
# H358
H358_72 = as.matrix(readMM(file = paste0(path,'/GSM4094254_H358.5.matrix.mtx.gz')))
colnames(H358_72) = read.table(paste0(path,'/GSM4094254_H358.5.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(H358_72) = read.table(paste0(path,'/GSM4094254_H358.5.genes.tsv.gz'), header=F)[,1]

colnames(H358_72) = str_sub(colnames(H358_72), end=-3)
keep = cell.annotation$X[cell.annotation$Hours==72 & cell.annotation$Line=='H358']
keep = str_sub(keep, start = 8)
H358_72 = H358_72[,keep]

# H2122
H2122_72 = as.matrix(readMM(file = paste0(path,'/GSM4094258_H2122.72.matrix.mtx.gz')))
colnames(H2122_72) = read.table(paste0(path,'/GSM4094258_H2122.72.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(H2122_72) = read.table(paste0(path,'/GSM4094258_H2122.72.genes.tsv.gz'), header=F)[,1]

colnames(H2122_72) = str_sub(colnames(H2122_72), end=-3)
keep = cell.annotation$X[cell.annotation$Hours==72 & cell.annotation$Line=='H2122']
keep = str_sub(keep, start = 9)
H2122_72 = H2122_72[,keep]

# SW1573
SW1573_72 = as.matrix(readMM(file = paste0(path,'/GSM4094262_SW1573.72.matrix.mtx.gz')))
colnames(SW1573_72) = read.table(paste0(path,'/GSM4094262_SW1573.72.barcodes.tsv.gz'), quote="\"", comment.char="")[,1]
rownames(SW1573_72) = read.table(paste0(path,'/GSM4094262_SW1573.72.genes.tsv.gz'), header=F)[,1]

colnames(SW1573_72) = str_sub(colnames(SW1573_72), end=-3)
keep = cell.annotation$X[cell.annotation$Hours==72 & cell.annotation$Line=='SW1573']
keep = str_sub(keep, start = 10)
SW1573_72 = SW1573_72[,keep]





# Merge -----

# table(rownames(H358) %in% rownames(H2122))

rows_keep = intersect(rownames(H358), rownames(SW1573))
rows_keep = intersect(rows_keep, rownames(H2122))

counts_t0 = cbind(H358[rows_keep,], SW1573[rows_keep,], H2122[rows_keep,])
cellgr_t0 = rep(c('H358','SW1573','H2122'), times = c(ncol(H358),ncol(SW1573),ncol(H2122)))
rm(H358,SW1573,H2122); gc()

counts_t4 = cbind(H358_4[rows_keep,], SW1573_4[rows_keep,], H2122_4[rows_keep,])
cellgr_t4 = rep(c('H358','SW1573','H2122'), times = c(ncol(H358_4),ncol(SW1573_4),ncol(H2122_4)))
rm(H358_4,SW1573_4,H2122_4); gc()

counts_t24 = cbind(H358_24[rows_keep,], SW1573_24[rows_keep,])
cellgr_t24 = rep(c('H358','SW1573'), times = c(ncol(H358_24),ncol(SW1573_24)))
rm(H358_24,SW1573_24); gc()

counts_t72 = cbind(H358_72[rows_keep,], SW1573_72[rows_keep,], H2122_72[rows_keep,])
cellgr_t72 = rep(c('H358','SW1573','H2122'), times = c(ncol(H358_72),ncol(SW1573_72),ncol(H2122_72)))
rm(H358_72,SW1573_72,H2122_72); gc()

rm(rows_keep)


#---- Filtering out genes with zero counts ----------------------------------

# pseudo time
tt = rep(c('t0','t4','t24','t72'),times=c(ncol(counts_t0),ncol(counts_t4),ncol(counts_t24),ncol(counts_t72)))

# merge the datasets
data_allt = cbind(counts_t0,counts_t4,counts_t24,counts_t72)

data_allt = data_allt[which(rowSums(data_allt)>0), ]



#---- Filtering out genes expressed in less than 5% of samples ----------------------------

rgenes = rownames(data_allt)[which(rowSums(data_allt>0) <= ncol(data_allt)*0.05)]

data_allt = data_allt[!(rownames(data_allt) %in% rgenes), ]


#---- Size factor ------------------------------------------------------

sf = librarySizeFactors(data_allt)





#---- ENTREZID ----------------------------------------------------------

# Genes in the pathways are encoded with ENTREZIDs, so we convert the gene names in our dataset accordingly.
# We download the mapping table from biomaRt.
# (For genes with multiple matches to the same ENTREZID, we discard the second occurrence)

genes_keep = rownames(data_allt)
 
# library(biomaRt)
# ensembl = useEnsembl(biomart = "genes", mirror='asia')
# ensembl = useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
# names_converter = getBM(filters = "ensembl_gene_id", values=genes_keep,
#                         attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
#                         mart=ensembl)
# 
# write.csv(names_converter, file=paste0(path,'/names_converter_norm_finale.csv'))

names_converter = read.csv(paste0(path,'/names_converter_norm_finale.csv'), he=T)[,-1]

# Keep only genes that have a unique Ensembl-ENTREZID correspondence, without duplicates
names_converter = names_converter[-which(is.na(names_converter$entrezgene_id)),] 
names_converter = names_converter[-which(duplicated(names_converter$ensembl_gene_id)),] 
names_converter = names_converter[-which(duplicated(names_converter$entrezgene_id)),] 

# Only genes from the dataset with a valid Ensembl–ENTREZID mapping are retained
# table(genes_keep %in% names_converter$ensembl_gene_id) 
data_allt = data_allt[genes_keep %in% names_converter$ensembl_gene_id,]

# Sort the names according to the order of the dataset
names_converter = names_converter[order(match(names_converter$ensembl_gene_id, rownames(data_allt))), ]

rownames(data_allt) = names_converter$entrezgene_id


#---- PATHWAYS -----------------------------------------------------------------------------------

# Download the pathways
pathways = pathways("hsapiens", "kegg")[] 
graphs = lapply(pathways, function(p) pathwayGraph(p)) 

for(i in 1:length(graphs)) graph::nodes(graphs[[i]]) = gsub("ENTREZID:","", graph::nodes(graphs[[i]]))

# Nodes not present in the dataset are removed
new_nodes = lapply(graphs, function(x){nodes(x)[which(nodes(x) %in% names_converter$entrezgene_id)]})
graphs = mapply(graphs, new_nodes, FUN = function(x,y){ x = graph::subGraph(y,x)})

# Only datasets including KRAS (3845) are selected
selector = function(x){ if(3845 %in% nodes(x)){ x } }
graphs = lapply(graphs, selector)

graphs = graphs[-which(sapply(graphs, is.null))]

# Only genes present in the pathways are retained in the dataset 
p_nodes = unlist(sapply(graphs, function(x){nodes(x)}), use.names = FALSE)
p_nodes = unique(p_nodes)

data_allt = data_allt[rownames(data_allt) %in% p_nodes,]

rm(p_nodes, new_nodes, pathways, i, selector, genes_keep, path)
rm(counts_t0, counts_t4, counts_t24, counts_t72)
rm(rgenes,rgenes0,rgenes4,rgenes24,rgenes72)
rm(cell.annotation)



#---- Theta -----------------------

library(edgeR)
d <- DGEList(counts=data_allt)
d1 <- estimateDisp(d, robust=TRUE)


rm(d)


save.image(file='D:\\Assegno\\Davide\\Two Stage\\Paper_twostage\\SCRIPT\\REAL DATA FINAL\\KRAS_dataset.RData')
