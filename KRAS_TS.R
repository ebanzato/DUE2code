
#---- Libraries -----------------------------------------------------

library(parallel)
library(Matrix)
library(graphite)
library(ggplot2)
library(tidyr)
library(gRbase)
library(igraph)
library(reshape2)
library(ggsci)
library(ggpubr)
library(stringr)
library(DUE2)
library(MASS)
library(xtable)

theme_set(theme_bw())


#---- Data ----------------------------------------------------------

load('D:/Assegno/Davide/Two Stage/PAPER/KRAS_dataset.RData')

gRAS  = graphs$`Ras signaling pathway`

gRAS = graph_from_graphnel(gRAS)
gRAS = subgraph(gRAS, V(gRAS)[!igraph::degree(gRAS)==0])

ugRAS = moralize(gRAS)
ugRAS = as.matrix(as_adjacency_matrix(ugRAS))

gene_names = names_converter[names_converter$entrezgene_id %in% colnames(ugRAS),2:3]
gene_names = gene_names[order(match(gene_names$entrezgene_id, rownames(ugRAS))), ]

df_ras = t(data_allt[rownames(data_allt) %in% colnames(ugRAS),])
df_ras = as.data.frame(df_ras[,order(match(colnames(df_ras), rownames(ugRAS)))])


# Aggiungo il nome comune
colnames(df_ras) = gene_names$hgnc_symbol
colnames(ugRAS) = rownames(ugRAS) = gene_names$hgnc_symbol




#---- Definisco i dataset e i sf ----

df_0 = df_ras[tt=='t0',] 
df_4 = df_ras[tt=='t4',] 
df_24 = df_ras[tt=='t24',] 
df_72 = df_ras[tt=='t72',] 


# Definisco i size factors

sf_0 = sf[tt=='t0']
sf_4 = sf[tt=='t4']
sf_24 = sf[tt=='t24']
sf_72 = sf[tt=='t72']


#---- Bilanciamo -----------------------------------------------

# tutti grandi quanto t24: 1006
set.seed(1301)
s0 = sample(1:nrow(df_0), nrow(df_24), replace=F)
s4 = sample(1:nrow(df_4), nrow(df_24), replace=F)
s72 = sample(1:nrow(df_72), nrow(df_24), replace=F)

# seleziono
df_0 = df_0[s0, ]; group0 = cellgr_t0[s0]; sf_0 = sf_0[s0]
df_4 = df_4[s4, ]; group4 = cellgr_t4[s4]; sf_4 = sf_4[s4]
df_72 = df_72[s72, ]; group72 = cellgr_t72[s72]; sf_72 = sf_72[s72]
group24 = cellgr_t24



#---- TS SEQUENCE ------------------------

pstime0_4 = as.factor(rep(c('t0','t4'), times=c(nrow(df_0),nrow(df_4))))
pstime4_24 = as.factor(rep(c('t4','t24'), times=c(nrow(df_4),nrow(df_24))))
pstime24_72 = as.factor(rep(c('t24','t72'), times=c(nrow(df_24),nrow(df_72))))


#----- t0 vs t4 ----
res_04 = DUE2(graph=ugRAS,
                data = rbind(df_0,df_4),
                group = pstime0_4,
                glm.family='nb',
                method.FDR='BH',
                method.FWER='holm',
                alpha=0.05,
                cell.group=c(group0,group4),
                sf = c(sf_0,sf_4),
                progressbar=TRUE)

rownames(res_04) = colnames(res_04)[-c(1:2)]


#----- t4 vs t24 ----
res_424 = DUE2(graph=ugRAS,
                 data = rbind(df_4,df_24),
                 group = pstime4_24,
                 glm.family='nb',
                 method.FDR='BH',
                 method.FWER='holm',
                 alpha=0.05,
                 cell.group=c(group4,group24),
                 sf = c(sf_4,sf_24),
                 progressbar=TRUE)

rownames(res_424) = colnames(res_424)[-c(1:2)]



#----- t24 vs t72 ----
res_2472 = DUE2(graph=ugRAS,
                 data = rbind(df_24,df_72),
                 group = pstime24_72,
                 glm.family='nb',
                 method.FDR='BH',
                 method.FWER='holm',
                 alpha=0.05,
                 cell.group=c(group24,group72),
                 sf = c(sf_24,sf_72),
                 progressbar=TRUE)

rownames(res_2472) = colnames(res_2472)[-c(1:2)]



#------ save/load ----
# save(list=c('res_04','res_424','res_2472'),
#      file = 'D:/Assegno/Davide/Two Stage/PAPER/KRAS_results_TS.RData')
load('D:/Assegno/Davide/Two Stage/PAPER/KRAS_results_TS.RData')


#------ PLOT ---------------------------------------------

plot.ts = function(results, coef.b, cond2, main.ts, ll){
  
  x11()
  par(mfrow=c(1,3))
  
  ## PLOT COND 1
  
  col_lab = ifelse(results[,1]==1, 'black', 'gray60')
  
  # Nodes
  col_nodes = ifelse(coef.b[,1]>0, 'red', 'blue')
  
  rgb.val.n = col2rgb(col_nodes)
  
  n.norm = abs(coef.b[,1])/max(abs(coef.b[,1]))
  
  t.col.n = sapply(1:ncol(rgb.val.n), function(i) {
    rgb(rgb.val.n[1,i], rgb.val.n[2,i], rgb.val.n[3,i],
        max = 255,
        alpha = n.norm[i]*255)}) 
  
  
  # Edges
  gRES = coef.b[,-1]
  gRES = (gRES + t(gRES))/2
  gRES = graph_from_adjacency_matrix(gRES, mode='undirected', weighted = T)
  
  E(gRES)$color = ifelse(E(gRES)$weight>0, 'red',
                         ifelse(E(gRES)$weight<0, 'blue', 'white'))
  
  E(gRES)$weight = abs(E(gRES)$weight)
  
  w.norm = E(gRES)$weight / max(E(gRES)$weight)
  
  rgb.val = col2rgb(E(gRES)$color)
  
  t.col.e = sapply(1:ncol(rgb.val), function(i) {
    rgb(rgb.val[1,i], rgb.val[2,i], rgb.val[3,i],
        max = 255,
        alpha = w.norm[i]*255)}) 


  ## PLOT COND 2
  
  # Nodes
  col_nodes2 = ifelse(cond2[,1]>0, 'red', 'blue')
  
  rgb.val.n2 = col2rgb(col_nodes2)
  
  n.norm2 = abs(cond2[,1])/max(abs(cond2[,1]))
  
  t.col.n2 = sapply(1:ncol(rgb.val.n2), function(i) {
    rgb(rgb.val.n2[1,i], rgb.val.n2[2,i], rgb.val.n2[3,i],
        max = 255,
        alpha = n.norm2[i]*255)}) 
  
  # Edges
  gRES2 = cond2[,-1]
  gRES2 = (gRES2 + t(gRES2))/2
  gRES2 = graph_from_adjacency_matrix(gRES2, mode='undirected', weighted = T)
  
  E(gRES2)$color = ifelse(E(gRES2)$weight>0, 'red',
                          ifelse(E(gRES2)$weight<0, 'blue', 'white'))
  E(gRES2)$weight = abs(E(gRES2)$weight)
  
  w.norm2 = E(gRES2)$weight / max(E(gRES2)$weight)
  
  rgb.val2 = col2rgb(E(gRES2)$color)
  
  t.col.e2 = sapply(1:ncol(rgb.val2), function(i) {
    rgb(rgb.val2[1,i], rgb.val2[2,i], rgb.val2[3,i],
        max = 255,
        alpha = w.norm2[i]*255)}) 
  
  
  ## LOCATION
  
  # PLOT
  color_1s = ifelse(results[,1]==1, 'black', 'gray60')
  color_2s_n = ifelse(results[,2]==1, 'tomato', 'gray90')
  color_2s_e = 'gray90'
  
  g4 = results[,-c(1:2)] + t(results[,-c(1:2)]) + 1
  
  wRAS = as.matrix(ugRAS * g4)
  wRAS = graph_from_adjacency_matrix(wRAS, weighted = T, mode='undirected')
  E(wRAS)$color = ifelse(E(wRAS)$weight==1, 'gray90',
                         ifelse(E(wRAS)$weight==2, 'orange', 'tomato'))
  
  
  # Differences
  plot(wRAS, vertex.color='white', layout=ll, 
       vertex.size=9, vertex.label.cex=0.6, 
       vertex.label.color = color_1s, # Stage 1
       vertex.frame.color = color_2s_n, # Stage 2: nodes
       edge.width=ifelse(E(wRAS)$color=='gray90',1,2), 
       vertex.frame.width = ifelse(color_2s_n=='gray90',1,2), 
       main=main.ts[1])
  
  # Cond 1
  plot(gRES, vertex.color='white', layout=ll, 
       vertex.size=9, vertex.label.cex=0.6, 
       vertex.label.color = col_lab, # Stage 1
       vertex.frame.color = t.col.n, # Stage 2: nodes
       edge.width = 10*E(gRES)$weight, 
       vertex.frame.width = 0.5 * abs(coef.b[,1]), 
       edge.color = t.col.e,
       main=main.ts[2])
  
  # Cond 2
  plot(gRES2, vertex.color='white', layout=ll, 
       vertex.size=9, vertex.label.cex=0.6, 
       vertex.label.color = col_lab, # Stage 1
       vertex.frame.color = t.col.n2, # Stage 2: nodes
       edge.width = 10*E(gRES2)$weight, 
       vertex.frame.width = 0.5 * abs(cond2[,1]), 
       edge.color = t.col.e2,
       main=main.ts[3])

}


#--------- t0 t4 ----

results_4 = res_04$g.diff
results_4[is.na(results_4)] = 0

coef_4 = res_04$base_coef
coef_4[is.na(coef_4)] = 0

delta_4 = res_04$delta_coef
delta_4[is.na(delta_4)] = 0
delta_4 = delta_4 * results_4[,-1]

cond2_4 = coef_4 + delta_4

#sort(abs(coef_4[,-1]),decreasing=T)[1:6]


#---- LAYOUT
gRES = coef_4[,-1]
gRES = (gRES + t(gRES))/2
gRES = graph_from_adjacency_matrix(gRES, mode='undirected', weighted = T)
E(gRES)$weight = abs(E(gRES)$weight)
set.seed(1301); ll = layout_with_fr(gRES)


plot.ts(results = results_4, coef.b = coef_4, cond2 = cond2_4,
        main.ts = c('Tests T0 vs T4', 'T0', 'T4'), ll = ll)




#--------- t4 t24 ----

results_24 = res_424$g.diff
results_24[is.na(results_24)] = 0

coef_24 = res_424$base_coef
coef_24[is.na(coef_24)] = 0

delta_24 = res_424$delta_coef
delta_24[is.na(delta_24)] = 0
delta_24 = delta_24 * results_24[,-1]

cond2_24 = coef_24 + delta_24

plot.ts(results = results_24, coef.b = coef_24, cond2 = cond2_24,
        main.ts = c('Tests T4 vs T24', 'T4', 'T24'), ll)



#--------- t24 t72 ----

results_72 = res_2472$g.diff
results_72[is.na(results_72)] = 0

coef_72 = res_2472$base_coef
coef_72[is.na(coef_72)] = 0

delta_72 = res_2472$delta_coef
delta_72[is.na(delta_72)] = 0
delta_72 = delta_72 * results_72[,-1]

cond2_72 = coef_72 + delta_72


plot.ts(results = results_72, coef.b = coef_72, cond2 = cond2_72,
        main.ts = c('Tests T24 vs T72', 'T24', 'T72'), ll=ll)


#------ TABLE differenze ---------------------------

#--------- t0 t4 ----

results = res_04$g.diff
results[is.na(results)] = 0

coef = res_04$base_coef
delta = res_04$delta_coef

# Intercette
i_diff = results[,2]
table(i_diff)

int_t0 = coef[i_diff==1,1]
int_d = delta[i_diff==1,1]
# table(names(int_t0) == names(int_d))
int_t4 = int_t0 + int_d

tab_int_04 = round(cbind('t_0'=int_t0, 't_4'=int_t4, 'delta'=int_d,
                         'exp(t_0)'=exp(int_t0), 'exp(t_4)'=exp(int_t4), 'exp(delta)'=exp(int_d)),3)
# tab_mu_04  = round(exp(tab_int_04),3)


# Associazioni
a_diff = results[,-c(1:2)]
id_diff= which(a_diff==1, arr.ind = T)

assoc_04 = data.frame(cbind('Gene1'=rownames(a_diff)[id_diff[,'row']],
                            'Gene2'=colnames(a_diff)[id_diff[,'col']]),
                      round(cbind('t0'=coef[,-1][id_diff],
                                  't4'=coef[,-1][id_diff]+delta[,-1][id_diff],
                                  'delta'=delta[,-1][id_diff],
                                  'exp(t0)'=exp(coef[,-1][id_diff]),
                                  'exp(t4)'=exp(coef[,-1][id_diff]+delta[,-1][id_diff]),
                                  'exp(delta)'=exp(delta[,-1][id_diff])),3))





#--------- t4 t24 ----

results = res_424$g.diff
results[is.na(results)] = 0

coef = res_424$base_coef
delta = res_424$delta_coef

# Intercette
i_diff = results[,2]
table(i_diff)

int_t0 = coef[i_diff==1,1]
int_d = delta[i_diff==1,1]
# table(names(int_t0) == names(int_d))
int_t4 = int_t0 + int_d

tab_int_424 = round(cbind('t_4'=int_t0, 't_24'=int_t4, 'delta'=int_d,
                          'exp(t_4)'=exp(int_t0), 'exp(t_24)'=exp(int_t4), 'exp(delta)'=exp(int_d)),3)
# tab_mu_424  = round(exp(tab_int_424),3)


# Associazioni
a_diff = results[,-c(1:2)]
id_diff= which(a_diff==1, arr.ind = T)

assoc_424 = data.frame(cbind('Gene1'=rownames(a_diff)[id_diff[,'row']],
                            'Gene2'=colnames(a_diff)[id_diff[,'col']]),
                       round(cbind('t4'=coef[,-1][id_diff],
                                   't24'=coef[,-1][id_diff]+delta[,-1][id_diff],
                                   'delta'=delta[,-1][id_diff],
                                   'exp(t4)'=exp(coef[,-1][id_diff]),
                                   'exp(t24)'=exp(coef[,-1][id_diff]+delta[,-1][id_diff]),
                                   'exp(delta)'=exp(delta[,-1][id_diff])),3))

# ass_e_424 = assoc_424
# ass_e_424[,3:5] = exp(ass_e_424[,3:5])



#--------- t24 t72 ----

results = res_2472$g.diff
results[is.na(results)] = 0

coef = res_2472$base_coef
delta = res_2472$delta_coef

# Intercette
i_diff = results[,2]
table(i_diff)

int_t0 = coef[i_diff==1,1]
int_d = delta[i_diff==1,1]
# table(names(int_t0) == names(int_d))
int_t4 = int_t0 + int_d

tab_int_2472 = round(cbind('t_24'=int_t0, 't_72'=int_t4, 'delta'=int_d,
                           'exp(t_24)'=exp(int_t0), 'exp(t_72)'=exp(int_t4), 'exp(delta)'=exp(int_d)),3)
# tab_mu_2472  = round(exp(tab_int_2472),3)


# Associazioni
a_diff = results[,-c(1:2)]
id_diff= which(a_diff==1, arr.ind = T)

assoc_2472 = data.frame(cbind('Gene1'=rownames(a_diff)[id_diff[,'row']],
                             'Gene2'=colnames(a_diff)[id_diff[,'col']]),
                        round(cbind('t24'=coef[,-1][id_diff],
                                    't72'=coef[,-1][id_diff]+delta[,-1][id_diff],
                                    'delta'=delta[,-1][id_diff],
                                    'exp(t24)'=exp(coef[,-1][id_diff]),
                                    'exp(t72)'=exp(coef[,-1][id_diff]+delta[,-1][id_diff]),
                                    'exp(delta)'=exp(delta[,-1][id_diff])),3))

# ass_e_2472 = assoc_2472
# ass_e_2472[,3:5] = exp(ass_e_2472[,3:5])



xtable(tab_int_04)
xtable(tab_int_424)
xtable(tab_int_2472)

xtable(assoc_04)
xtable(assoc_424)
xtable(assoc_2472)








#------ PLOT differenze ----------------------------------


#--- FUNCTIONS ----
plot.due2 = function(results, gpath, ll=layout_in_circle, main= 'DUE2'){
  
res = results$g.diff
res[is.na(res)] = 0

delta = results$delta_coef
delta[is.na(delta)] = 0
delta = delta * res[,-1]

delta[which(res[,-1]==0, arr.ind = T)] = 0
delta[,-1] = (delta[,-1]+t(delta[,-1]))/2



## LOCATION

# PLOT
color_1s = ifelse(res[,1]==1, 'black', 'gray60')
color_2s_n = ifelse(res[,2]==1, 'forestgreen', 'gray90')
color_2s_e = 'gray90'

gres = res[,-c(1:2)] + t(res[,-c(1:2)]) + 1

g = as.matrix(gpath * gres)
g = graph_from_adjacency_matrix(g, weighted = T, mode='undirected')
E(g)$color = ifelse(E(g)$weight==1, 'gray90',
                       ifelse(E(g)$weight==2, 'chartreuse2', 'forestgreen'))


# Differences
plot(g, vertex.color='white', layout=ll, 
     vertex.size=9, vertex.label.cex=0.6, 
     vertex.label.color = color_1s, # Stage 1
     vertex.frame.color = color_2s_n, # Stage 2: nodes
     edge.width=ifelse(E(g)$color=='gray90',1,2), 
     vertex.frame.width = ifelse(color_2s_n=='gray90',1,2), 
     main=main)
}

plot.coef = function(results, ll=layout_in_circle, main=expression(Coefficients ~ Delta)){
  
  res = results$g.diff
  res[is.na(res)] = 0
  
  delta = results$delta_coef
  delta[is.na(delta)] = 0
  delta = delta * res[,-1]
  
  delta[which(res[,-1]==0, arr.ind = T)] = 0
  delta[,-1] = (delta[,-1]+t(delta[,-1]))/2
  
# PLOT
color_nlab = res[,-1]
color_nlab[,-1] = color_nlab[,-1] + t(color_nlab[,-1])
color_lab = ifelse(rowSums(color_nlab)==0, 'gray60', 'black')

v_color = ifelse(delta[,1]>0, 'tomato1', 
                 ifelse(delta[,1]<0, 'cornflowerblue', 'gray85'))

g = graph_from_adjacency_matrix(delta[,-1], weighted = T, mode='undirected')
E(g)$color = ifelse(E(g)$weight>0, 'tomato1',
                    ifelse(E(g)$weight<0, 'cornflowerblue', 'gray8'))



plot(g, vertex.color='white', layout=ll, 
     vertex.size=9, vertex.label.cex=0.6, 
     vertex.label.color = color_lab, 
     vertex.frame.color = v_color, 
     vertex.frame.width = ifelse(v_color=='gray85', 1.5, abs(delta[,1])*4),
     edge.width=abs(E(g)$weight)*30, main=main)


}




#--- PLOT FINALI ----

#---- LAYOUT
coef_4 = res_04$base_coef
coef_4[is.na(coef_4)] = 0
gRES = coef_4[,-1]
gRES = (gRES + t(gRES))/2
gRES = graph_from_adjacency_matrix(gRES, mode='undirected', weighted = T)
E(gRES)$weight = abs(E(gRES)$weight)
set.seed(1301); ll = layout_with_fr(gRES)

x11()
par(mfrow=c(1,2))

# t0 t4
plot.due2(res_04, ugRAS, ll=ll, main=expression(DUE2 - T[0] ~ vs ~ T[4]))
plot.coef(res_04, ll=ll)

# t4 t24
plot.due2(res_424, ugRAS, ll=ll, main=expression(DUE2 - T[4] ~ vs ~ T[24]))
plot.coef(res_424, ll=ll)

# t24 t72
plot.due2(res_2472, ugRAS, ll=ll, main=expression(DUE2 - T[24] ~ vs ~ T[72]))
plot.coef(res_2472, ll=ll)





