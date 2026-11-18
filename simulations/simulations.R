#---- Libraries -----------------------------------------

library(ggsci)
library(ggpubr)
library(igraph)
library(ggplot2)
library(spsUtil)
library(parallel)
library(reshape2)
library(ModPGMInference)
library(DUE2)

# Functions for data simulation and summary measures
source('path/to/sim_functions.R')

theme_set(theme_bw())
my.theme = theme(axis.text=element_text(size=20), 
                 plot.title=element_text(size=22),
                 axis.title=element_text(size=20), 
                 legend.text=element_text(size=20),
                 legend.title=element_text(size=22),  
                 strip.text = element_text(size=20))


#---- Scenario ------------------------------------------

# Graph: subgraph of the KEGG Ras signaling pathway (29 nodes, 43 edges)

p = 29

# Edge list
edges = matrix(c( 1, 2,   1,10,   1,16,
                  2,10,   2,16,
                  3, 5,   3,11,   3,13,   3,14,   3,15,   3,16,   3,17,   3,19,
                  4,16,
                  5,17,
                  6, 7,   6, 8,   6, 9,   6,20,   6,21,   6,22,
                  7, 8,   7, 9,   7,20,   7,21,   7,22,
                  8, 9,   8,13,
                  9,13,
                 10,23,
                 11,14,  11,15,
                 12,16,  12,24,  12,25,  12,26,  12,27,
                 14,15,  14,18,  14,28,  14,29,
                 18,28,  18,29),
               ncol = 2, byrow = TRUE)

# Adjacency matrix
g = matrix(0, p, p)
g[edges] = 1
g[edges[, c(2,1)]] = 1
colnames(g) = rownames(g) = paste0('v',1:p)

e = nrow(edges)

pg = graph_from_adjacency_matrix(g, mode='undirected')

# PLOT
pg = g; colnames(pg) = rownames(pg) = 1:p
pg = graph_from_adjacency_matrix(pg, mode='undirected')

V(pg)$color = 'white'
V(pg)$frame.color = rep('gray59', 29); V(pg)$frame.color[c(3,8,23)] = 'tomato'
V(pg)$frame.width = rep(1, 29); V(pg)$frame.width[c(3,8,23)] = 3
E(pg)$color = rep('gray59',43); E(pg)$color[c(4,35)] = 'tomato'
E(pg)$width = rep(1,43); E(pg)$width[c(4,35)] = 3

{x11()
  plot(pg, vertex.size=12, vertex.label.color='black', 
       vertex.label.cex=1.5, layout=layout_with_kk)
}


#---- Set parameters -------------------------

set.seed(1)
th1 = runif(p, 0.1, 0.5)
names(th1) = colnames(g)
th1[c('v3','v8','v23')] = 0.1

set.seed(2)
phi1 = g; phi1[lower.tri(phi1)] = 0
phi1[phi1==1] = runif(length(phi1[phi1==1]),-0.3,-0.1) 
phi1 = phi1 + t(phi1)
phi1['v12','v25'] = phi1['v25','v12'] = -0.1 
phi1['v2','v10'] = phi1['v10','v2'] = -0.1 

#sort(igraph::degree(pg))

## DIFFERENCES
thx2 = thx3 = thx4 = th1
thx2[c('v3','v8','v23')] = 0.1 * 2
thx3[c('v3','v8','v23')] = 0.1 * 3
thx4[c('v3','v8','v23')] = 0.1 * 4

phix2 = phix3 = phix4 = phi1
phix2['v12','v25'] = phix2['v25','v12'] = -0.1  * 2
phix3['v12','v25'] = phix3['v25','v12'] = -0.1  * 3
phix4['v12','v25'] = phix4['v25','v12'] = -0.1  * 4
phix2['v2','v10'] = phix2['v10','v2'] = -0.1  * 2
phix3['v2','v10'] = phix3['v10','v2'] = -0.1  * 3
phix4['v2','v10'] = phix4['v10','v2'] = -0.1  * 4

true.diff = cbind('th' = ifelse(th1-thx2==0,0,1),
                  ifelse(phi1-phix2==0,0,1))


#---- Simulations --------------------------------------

R = rep(30,p)
nSim = 1000
N = c(500,1000,1500,2000,3000,5000)

glm.family = 'poisson'
alpha = 0.05
method.FDR = 'BH'
method.FWER = 'holm'

ncores = 7

set.seed(1301)

#### k = 2  -----
param_list = list('th1'=th1,'th2'=thx2,'phi1'=phi1,'phi2'=phix2,'R'=R)

cl = makeCluster(ncores)
clusterExport(cl, ls())
clusterEvalQ(cl, c(library('spsUtil'), library('ModPGMInference'), library('DUE2')))

start = Sys.time()
res = lapply(N, function(n){
  
  print(paste0('2, n=',n))
  
  res_sqrt = parLapply(cl, 1:nSim, function(sim){
    
    # Groups
    group = rep(c(0,1), each=n)
    
    # Sim data
    data = simdata_tpgm(n, param_list)
    colnames(data) = colnames(g)
    
    # Two stage
    res_p = DUE2(graph = g, data = data, group = group,
                 glm.family = glm.family, alpha = alpha,
                 method.FDR = method.FDR, method.FWER = method.FWER)
    
    # OUT
    return(res_p$g.diff)
    
  })
  
  res_sqrt
})
names(res) = paste0('n', N)
Sys.time() - start

stopCluster(cl)


save(res, file='path/to/results/TS_tpgm_2.RData')

#### k = 3  -----
param_list = list('th1'=th1,'th2'=thx3,'phi1'=phi1,'phi2'=phix3,'R'=R)

cl = makeCluster(ncores)
clusterExport(cl, ls())
clusterEvalQ(cl, c(library('spsUtil'), library('ModPGMInference'), library('DUE2')))

start = Sys.time()
res = lapply(N, function(n){
  
  print(paste0('3, n=',n))
  
  res_sqrt = parLapply(cl, 1:nSim, function(sim){
    
    # Groups
    group = rep(c(0,1), each=n)
    
    # Sim data
    data = simdata_tpgm(n, param_list)
    colnames(data) = colnames(g)
    
    # Two stage
    res_p = DUE2(graph = g, data = data, group = group,
                 glm.family = glm.family, alpha = alpha,
                 method.FDR = method.FDR, method.FWER = method.FWER)
    
    # OUT
    return(res_p$g.diff)
    
  })
  
  res_sqrt
})
names(res) = paste0('n', N)
Sys.time() - start

stopCluster(cl)


save(res, file='path/to/results/TS_tpgm_3.RData')



#### k = 4  -----
param_list = list('th1'=th1,'th2'=thx4,'phi1'=phi1,'phi2'=phix4,'R'=R)

cl = makeCluster(ncores)
clusterExport(cl, ls())
clusterEvalQ(cl, c(library('spsUtil'), library('ModPGMInference'), library('DUE2')))

start = Sys.time()
res = lapply(N, function(n){
  
  print(paste0('4, n=',n))
  
  res_sqrt = parLapply(cl, 1:nSim, function(sim){
    
    # Groups
    group = rep(c(0,1), each=n)
    
    # Sim data
    data = simdata_tpgm(n, param_list)
    colnames(data) = colnames(g)
    
    # Two stage
    res_p = DUE2(graph = g, data = data, group = group,
                 glm.family = glm.family, alpha = alpha,
                 method.FDR = method.FDR, method.FWER = method.FWER)
    
    # OUT
    return(res_p$g.diff)
    
  })
  
  res_sqrt
})
names(res) = paste0('n', N)
Sys.time() - start

stopCluster(cl)


save(res, file='path/to/results/TS_tpgm_4.RData')




#---- FIGURE ---------------------------------------------------

tn_lrt = ifelse(rowSums(true.diff)==0, 0, 1) 
lrt = names(tn_lrt[which(tn_lrt>0)])


#### results for k = 2 -----------------------------------------------------

load('path/to/results/TS_tpgm_2.RData')

ofdr_k_2 =  sapply(res, function(x) OFDR(x, lrt))
pw_k_2 = sapply(res,  function(x) power.1stage(x, lrt))

fwer_k_2 = t(sapply(res, function(dn) fwer.2stage(dn, true.diff)))
table(fwer_k_2>0.05)

pw2_k = lapply(res, function(dn) power.2stage(dn, true.diff))

pw2_k_2 = sapply(pw2_k, function(dn) dn[which(true.diff==1, arr.ind=T)])
rownames(pw2_k_2) = rownames(which(true.diff==1, arr.ind=T))



#### results for k = 3 -----------------------------------------------------

load('path/to/results/TS_tpgm_3.RData')

ofdr_k_3 =  sapply(res, function(x) OFDR(x, lrt))
pw_k_3 = sapply(res,  function(x) power.1stage(x, lrt))

fwer_k_3 = t(sapply(res, function(dn) fwer.2stage(dn, true.diff)))
table(fwer_k_3>0.05)

pw2_k = lapply(res, function(dn) power.2stage(dn, true.diff))

pw2_k_3 = sapply(pw2_k, function(dn) dn[which(true.diff==1, arr.ind=T)])
rownames(pw2_k_3) = rownames(which(true.diff==1, arr.ind=T))



#### results for k = 4 -----------------------------------------------------

load('path/to/results/TS_tpgm_4.RData')

ofdr_k_4 =  sapply(res, function(x) OFDR(x, lrt))
pw_k_4 = sapply(res,  function(x) power.1stage(x, lrt))

fwer_k_4 = t(sapply(res, function(dn) fwer.2stage(dn, true.diff)))
table(fwer_k_4>0.05)

pw2_k = lapply(res, function(dn) power.2stage(dn, true.diff))

pw2_k_4 = sapply(pw2_k, function(dn) dn[which(true.diff==1, arr.ind=T)])
rownames(pw2_k_4) = rownames(which(true.diff==1, arr.ind=T))


#------ Summary figures ---------------------------------------------------

#### OFDR ----
dataplot = cbind('2'=ofdr_k_2, '3'=ofdr_k_3, '4'=ofdr_k_4)
rownames(dataplot) = N

dataplot = melt(dataplot)
colnames(dataplot) = c('N', 'k', 'OFDR')
dataplot$k = as.factor(dataplot$k)

ofdr = ggplot(dataplot, aes(x=N, y=OFDR, col=k)) +
  geom_line(linewidth=1) + geom_point() + ylim(0,0.1) +
  ylab('') + ggtitle('OFDR') + my.theme + scale_color_aaas()+
  geom_hline(yintercept = 0.05, linetype = "dashed")


#### FWER ----
rownames(fwer_k_2) = rownames(fwer_k_3) = rownames(fwer_k_4) = N
dataplot = rbind(cbind(melt(fwer_k_2)[,-2], 'k'='2'),
                 cbind(melt(fwer_k_3)[,-2], 'k'='3'),
                 cbind(melt(fwer_k_4)[,-2], 'k'='4'))
colnames(dataplot) = c('N', 'FWER', 'k')
dataplot$k = as.factor(dataplot$k)
dataplot$N = as.factor(dataplot$N)

fwer = ggplot(dataplot, aes(x=N, y=FWER, col=k)) +
  geom_boxplot(linewidth=1) + ylim(0,0.1) +
  ylab('') + ggtitle('FWER') + my.theme + scale_color_aaas()+
  geom_hline(yintercept = 0.05, linetype = "dashed")

{
x11()
ggarrange(ofdr, fwer, common.legend = T, legend='right')
}


#### Power ----
rownames(pw_k_2) = rownames(pw_k_3) = rownames(pw_k_4) = substr(rownames(pw_k_4), 2,3)
colnames(pw_k_2) = colnames(pw_k_3) = colnames(pw_k_4) = N
rownames(pw2_k_2) = rownames(pw2_k_3) = rownames(pw2_k_4) = substr(rownames(pw2_k_4), 2,3)
colnames(pw2_k_2) = colnames(pw2_k_3) = colnames(pw2_k_4) = N


# nodes
dp1 = cbind(rbind(cbind(melt(pw_k_2[c(2:3,6),]), 'k'='k=2'),
           cbind(melt(pw_k_3[c(2:3,6),]), 'k'='k=3'),
           cbind(melt(pw_k_4[c(2:3,6),]), 'k'='k=4')),'Stage'='Stage 1')
dp2 = cbind(rbind(cbind(melt(pw2_k_2[c(1:3),]), 'k'='k=2'),
           cbind(melt(pw2_k_3[c(1:3),]), 'k'='k=3'),
           cbind(melt(pw2_k_4[c(1:3),]), 'k'='k=4')),'Stage'='Stage 2')
dp = rbind(dp1, dp2)

colnames(dp) = c('Node', 'N', 'Power', 'k', 'Stage')
dp$Node = as.factor(dp$Node)

pw_n = ggplot(dp, aes(x=N, y=Power, col=Node)) +
  geom_line(linewidth=1) + geom_point() + #ylim(0,0.1) +
  ylab('') + my.theme + scale_color_d3() +
  facet_wrap(~ Stage + k)


{x11()
pw_n
}


# edges
dp1 = cbind(rbind(cbind(melt(pw_k_2[-c(2:3,6),]), 'k'='k=2'),
                  cbind(melt(pw_k_3[-c(2:3,6),]), 'k'='k=3'),
                  cbind(melt(pw_k_4[-c(2:3,6),]), 'k'='k=4')),'Stage'='Stage 1')
dp2 = cbind(rbind(cbind(melt(pw2_k_2[-c(1:3),]), 'k'='k=2'),
                  cbind(melt(pw2_k_3[-c(1:3),]), 'k'='k=3'),
                  cbind(melt(pw2_k_4[-c(1:3),]), 'k'='k=4')),'Stage'='Stage 2')
dp = rbind(dp1, dp2)

colnames(dp) = c('Node', 'N', 'Power', 'k', 'Stage')
dp$Node = as.factor(dp$Node)

pw_e = ggplot(dp, aes(x=N, y=Power, col=Node)) +
  geom_line(linewidth=1) + geom_point() + #ylim(0,0.1) +
  ylab('') + my.theme + scale_color_nejm() +
  facet_wrap(~ Stage + k)


{x11()
pw_e
  }
