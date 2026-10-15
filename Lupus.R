# Load data

library(muscData)
library(scater)
library(scrapper)
library(readr)

sce <- Kang18_8vs8()

# Filter out doublets and visualize the data.

sce <- sce[,sce$multiplets == "singlet"]
plotTSNE(sce, color_by = "cell")
plotTSNE(sce, color_by = "stim")

# B cell
bcell <- sce[,which(sce$cell == "B cells")]

bcell <- scrapper::normalizeRnaCounts.se(bcell)
bcell <- runPCA(bcell)
bcell <- runTSNE(bcell)
bcell$ind <- factor(bcell$ind)
table(bcell$ind)

plotTSNE(bcell, color_by = "stim")

## Differential expression
library(glmGamPoi)

bcell <- bcell[rowSums(counts(bcell)) > 5,]

reduced_sce <- pseudobulk(bcell, group_by = vars(ind, condition = stim), n_cells = n())

fit <- glm_gp(reduced_sce, design = ~ condition + ind, size_factor = "ratio", verbose = TRUE)

res <- test_de(fit, contrast = cond(condition = "stim") - cond(condition = "ctrl"))
head(res)
table(res$adj_pval <= 0.05)

degs <- res[res$adj_pval <= 0.05, "name"]
head(degs)

write.csv(res, file="Bcells_de_results.csv")

## Pathway selection
library(AnnotationDbi)
library(org.Hs.eg.db)

genes_keep <- rownames(bcell)

names_converter <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = unique(genes_keep),
  keytype  = "SYMBOL",
  columns  = c("ENSEMBL", "ENTREZID")
)

names(names_converter) <- c(
  "external_gene_name", "ensembl_gene_id", "entrezgene_id"
)

dim(names_converter)
head(names_converter)

library(tidyverse)
df <- inner_join(data.frame(external_gene_name = genes_keep), names_converter)

TS_counts = assay(bcell)[df$external_gene_name,]
rownames(TS_counts) = df$entrezgene_id

# Download and select pathways.

library(graphite)
pathways = pathways("hsapiens", "kegg")[] # !!! 'hsapiens' if human
graphs = lapply(pathways, function(p) pathwayGraph(p)) 

for(i in 1:length(graphs)) graph::nodes(graphs[[i]]) = gsub("ENTREZID:","", graph::nodes(graphs[[i]]))

# Nodes not in the dataset are removed
new_nodes = lapply(graphs, function(x){graph::nodes(x)[which(graph::nodes(x) %in% df$entrezgene_id)]})
graphs = mapply(graphs, new_nodes, FUN = function(x,y){ x = graph::subGraph(y,x)})

# Remove empty pathways
graphs = graphs[-which(sapply(graphs, function(x) length(graph::nodes(x))==0))]

# Only genes present in the pathways are retained in the dataset 
p_nodes = unlist(sapply(graphs, function(x){graph::nodes(x)}), use.names = FALSE)
p_nodes = unique(p_nodes)

TS_counts = TS_counts[rownames(TS_counts) %in% p_nodes, ]
# length(p_nodes); dim(TS_counts)

#table(names_converter$entrezgene_id %in% rownames(TS_counts))
df = df[df$entrezgene_id %in% rownames(TS_counts),]
dim(TS_counts)
dim(df)

# Create graphs.

library(gRbase)
library(igraph)
# Graphs moralization
m.graphs = lapply(graphs, function(x) gRbase::moralize(graph_from_graphnel(x)))

# Get the adjacency matrices
a.graph = lapply(m.graphs, function(x) as.matrix(as_adjacency_matrix(x)))


# Convert gene identifiers to common gene names
# head(names_converter$entrezgene_id); rownames(TS_counts)[1:6]
df$external_gene_name = gsub("-", "", df$external_gene_name)

rownames(TS_counts) = df$external_gene_name

# Change gene names also in the graphs
a.graph = lapply(a.graph, function(x){
  sub.names = df[df$entrezgene_id %in% rownames(x),]
  sub.names = sub.names[match(rownames(x), sub.names$entrezgene_id),]
  rownames(x) = colnames(x) = sub.names$external_gene_name
  x
})

# Filter pathways.

ll = sapply(a.graph, ncol)#; unname(ll)
a.graph = a.graph[ll >= 10]
length(a.graph)

## Two stage approach
library(DUE2)

SF <- librarySizeFactors(bcell)
head(SF)

for(j in seq_along(a.graph)){

  print(j)

  # graph
  g = a.graph[[j]]
  # dataset
  data.ts = as.matrix(t(TS_counts[colnames(g),]))

  res = tryCatch(
    {
      DUE2(
        graph = g,
        data = data.ts,
        group = as.factor(bcell$stim),
        cell.group = as.factor(bcell$ind),
        glm.family = "nb",
        method.FDR = "BH",
        method.FWER = "holm",
        alpha = 0.05,
        sf = SF,
        progressbar = TRUE
      )
    },
    error = function(e) NA
  )

  saveRDS(res, file = paste0("res/kang_paired", j, ".rds"))
}

names(res_all) = names(a.graph)

# Remove NA's and sort.

library(ggplot2)
theme_set(theme_bw())
my.theme = theme(axis.text=element_text(size=20), plot.title=element_text(size=22),
                 axis.title=element_text(size=20), legend.text=element_text(size=20),
                 legend.title=element_text(size=22))
table(is.na(res_all))

res_all = res_all[!is.na(res_all)]

mean_dis = sort(sapply(res_all, function(x) mean(x[[1]][,1])), decreasing = TRUE)
df = data.frame('pathway'=names(mean_dis),'value'=mean_dis)
df = df[order(df$value, decreasing = TRUE),]
df$x = 1:nrow(df)

p <- ggplot(df, aes(x = x, y = value)) +
  geom_point(size = 1.5, col='slateblue4') +
  labs(x = "", y = "") + ggtitle('Proportion of rejected genes in stage 1') +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  ) + my.theme

p

ggsave("Bcell_ranking.pdf", plot = p, width = 8, height = 6)

pathway_final <- data.frame(pathway=names(mean_dis), prop=unname(round(mean_dis, 4)))
head(pathway_final)

write_csv(pathway_final, file = 'Bcell_pathways_paired.csv')

# Compare with DE results.

path_genes <- unique(unlist(lapply(res_all, function(x) rownames(x[[1]]))))
disreg <- unique(unlist(lapply(res_all, function(x) rownames(x[[1]])[which(res_all[[1]][[1]][,1]>0)])))

length(disreg)
length(degs)
length(path_genes)
length(intersect(degs, path_genes))
length(intersect(degs, disreg))
intersect(degs, disreg)

# Plot.

set.seed(13)
path_names <- c("Influenza A", "Herpes simplex virus 1 infection", "Hepatitis C", "Measles",  "Coronavirus disease - COVID-19", "Epstein-Barr virus infection")
pathways = lapply(res_all[path_names], function(x) x[[1]][,1])

all_genes = sort(unique(unlist(lapply(pathways, names))))

# matrice pathways x geni (NA = gene assente)
mat = matrix(NA,
             nrow = length(pathways),
             ncol = length(all_genes),
             dimnames = list(names(pathways), all_genes))

for (pw in names(pathways)) {
  mat[pw, names(pathways[[pw]])] = pathways[[pw]]
}

# percentuale di 1 per pathway
perc_ones = rowSums(mat == 1, na.rm = TRUE) / rowSums(!is.na(mat))

# selezione top pathways
top_pw = names(sort(perc_ones, decreasing = TRUE))
mat_pw = mat[top_pw, , drop = FALSE]

# frequenza dei geni = 1 nei top pathways
gene_freq = colSums(mat_pw == 1, na.rm = TRUE)

# selezione top 15 geni
top_genes = names(gene_freq[gene_freq>=2])
# top_genes = names(gene_freq)[1:60]
mat_pw = mat_pw[, top_genes, drop = FALSE]

# codifica colori
# 0 = assente, 1 = presente ma 0, 2 = presente e 1
mat_col = matrix(0,
                 nrow = nrow(mat_pw),
                 ncol = ncol(mat_pw),
                 dimnames = dimnames(mat_pw))

mat_col[!is.na(mat_pw) & mat_pw == 0] = 1
mat_col[!is.na(mat_pw) & mat_pw == 1] = 2

library(pheatmap)

ph = pheatmap(t(mat_col),
              cluster_rows = FALSE,
              cluster_cols = FALSE,
              color = c("grey90", "skyblue3", "darkblue"),
              legend_breaks = c(0, 1, 2),
              legend_labels = c("Absent", "Not significant", "Significant"),
              border_color = "white",
              silent = TRUE, fontsize = 8)


pdf("Bcell_heatmap.pdf")
ph
dev.off()
