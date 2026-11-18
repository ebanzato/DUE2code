
#---- Libraries -----------------------------------------------------

library(parallel)
library(Matrix)
library(gRbase)
library(igraph)
library(edgeR)

# DUE2 is available on GitHub. To install it:
# install.packages('remotes')
# remotes::install_github('ebanzato/DUE2')
library(DUE2)


#---- Paths ---------------------------------------------------------
# Set these paths before running the script.
# 'data_file' is the .RData file created by KRAS.R (its 'out_file').
# 'res_file'  is where the simulation results are saved.

data_file = 'path/to/output/KRAS_dataset.RData'
res_file  = 'path/to/results/H0_kras.RData'


#---- Data -------------------------------------------------------------------

load(data_file)

gRAS  = graphs$`Ras signaling pathway`

gRAS = graph_from_graphnel(gRAS)
gRAS = subgraph(gRAS, V(gRAS)[!igraph::degree(gRAS)==0])

ugRAS = moralize(gRAS)
ugRAS = as.matrix(as_adjacency_matrix(ugRAS))

gene_names = names_converter[names_converter$entrezgene_id %in% colnames(ugRAS),2:3]
gene_names = gene_names[order(match(gene_names$entrezgene_id, rownames(ugRAS))), ]

df_ras = t(data_allt[rownames(data_allt) %in% colnames(ugRAS),])
df_ras = as.data.frame(df_ras[,order(match(colnames(df_ras), rownames(ugRAS)))])

# Use gene symbols as names
colnames(df_ras) = gene_names$hgnc_symbol
colnames(ugRAS) = rownames(ugRAS) = gene_names$hgnc_symbol

# Select only t=0
df_0 = df_ras[tt=='t0',]#; df_0 = df_0[-1,]
sf_0 = sf[tt=='t0']#; sf_0 = sf_0[-1]
cg_0 = cellgr_t0#[-1]

dim(df_0)
length(sf_0)
length(cg_0)

# rm(list = ls()[!ls() %in% c('df_0','sf_0','cg_0','ugRAS')])


#---- TS n ------------------------

# Define the condition (group labels, randomly permuted at each simulation)
cond = c(1, rep(c(0,1), each=nrow(df_0)/2))

nSim = 1000

set.seed(2706)

cl = makeCluster(5)
clusterExport(cl, ls())
clusterEvalQ(cl, c(library('utils'), library('DUE2')))

start = Sys.time()

res_h0 = parLapply(cl, 1:nSim, function(sim){

  cond = sample(cond, length(cond), replace=FALSE)

  res_nb = DUE2(graph=ugRAS,
                   data = df_0,
                   group = cond,
                   glm.family='nb',
                   method.FDR='BH',
                   method.FWER='holm',
                   alpha=0.05,
                   cell.group= cg_0,
                   sf = sf_0,
                   progressbar=TRUE)

  res_p = DUE2(graph=ugRAS,
                   data = df_0,
                   group = cond,
                   glm.family='poisson',
                   method.FDR='BH',
                   method.FWER='holm',
                   alpha=0.05,
                   cell.group= cg_0,
                   sf = sf_0,
                   progressbar=TRUE)

  return(list('nb'=res_nb, 'pois'=res_p))

})

stopCluster(cl)

Sys.time() - start

save(res_h0, file = res_file)

# To skip the simulation, load the saved results instead
# load(res_file)


# Results
res_nb = lapply(res_h0, function(x) x$nb)
res_p = lapply(res_h0,  function(x) x$pois)

mean(sapply(res_p,  function(x) ifelse(sum(x$g.diff[,1])>0,1,0)))
mean(sapply(res_nb, function(x) ifelse(sum(x$g.diff[,1])>0,1,0)))



#---- PLOT --------

pv_nb = t(sapply(res_nb, function(x) x$delta_pvalues[,1]))
pv_p  = t(sapply(res_p,  function(x) x$delta_pvalues[,1]))

#dim(pv_p)


fdr_emp = sapply(c(0.001, 0.01, 0.05, 0.1), function(a){
  nb = apply(pv_nb, 1, function(i) sum(p.adjust(i, method='BH')<a)>0)
  p  = apply(pv_p,  1, function(i) sum(p.adjust(i, method='BH')<a)>0)
  c('nb'=mean(nb),'pois'=mean(p))
})

dataplot = data.frame(
  alpha = rep(c(0.001, 0.01, 0.05, 0.10), 2),
  method = rep(c("Poisson", "NB"), each = 4),
  fdr_emp = c(fdr_emp[2,],fdr_emp[1,])
)








dev.new()
{
  par(mfrow=c(1,3),
      cex.lab = 2,        # axis.title
      cex.axis = 1.7,     # axis.text
      cex.main = 2.1,     # plot.title
      las = 1.5,          # axis labels orientation
      mar = c(5, 5, 3, 3))

cols = c("#1F77B4", "#FF7F0E")
methods = unique(dataplot$method)

xlim = c(0,0.1)
ylim = c(0,1)

# wider right margin
# par(mar = c(5, 5, 2, 12))#, las = 1, cex.lab = 1.3, cex.axis = 1.2)

# empty plot
plot(NA, xlim = xlim, ylim = ylim, bty = "l",
     xlab = "α", ylab = "Empirical FDR")

# diagonal line
abline(a = 0, b = 1,lty = 2, lwd = 2)

# lines and points for each method
for(i in 1:2) {
  m = methods[i]
  sub = dataplot[dataplot$method == m, ]
  sub = sub[order(sub$alpha), ]
  lines(sub$alpha, sub$fdr_emp, col = cols[i],lwd = 2)
  points(sub$alpha,sub$fdr_emp,col = cols[i], pch = 16, cex=2)
}

## legend
legend(x = 0.075,
       y = 0.75,
       legend = methods,
       col = cols,
       lwd = 2,
       pch = 16,
       cex=2, seg.len = 1,
       title = "Method",
       xpd = TRUE,
       bty = "n")


hist(pv_nb, col="gray90", prob=T, xlab='p-value',main='Negative binomial',ylim=c(0,3.3),cex=2)
hist(pv_p, col="gray90", prob=T, xlab='p-value',main='Poisson',ylim=c(0,3.3),cex=2)

}
