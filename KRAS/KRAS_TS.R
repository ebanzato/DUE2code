
#---- Libraries -----------------------------------------------------

library(graphite)
library(gRbase)
library(igraph)
library(xtable)

# DUE2 is available on GitHub. To install it:
# install.packages('remotes')
# remotes::install_github('ebanzato/DUE2')
library(DUE2)


#---- Paths ---------------------------------------------------------
# Set these paths before running the script.
# 'data_file' is the .RData file created by KRAS.R (its 'out_file').
# 'res_file'  is where the DUE2 results are saved.

data_file = 'path/to/output/KRAS_dataset.RData'
res_file  = 'path/to/results/KRAS_results_TS.RData'


#---- Data ----------------------------------------------------------

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




#---- Define datasets and size factors ----

df_0 = df_ras[tt=='t0',]
df_4 = df_ras[tt=='t4',]
df_24 = df_ras[tt=='t24',]
df_72 = df_ras[tt=='t72',]


# Size factors

sf_0 = sf[tt=='t0']
sf_4 = sf[tt=='t4']
sf_24 = sf[tt=='t24']
sf_72 = sf[tt=='t72']


#---- Balancing -----------------------------------------------

# all time points are subsampled to the size of t24 (1006 cells)
set.seed(1301)
s0 = sample(1:nrow(df_0), nrow(df_24), replace=F)
s4 = sample(1:nrow(df_4), nrow(df_24), replace=F)
s72 = sample(1:nrow(df_72), nrow(df_24), replace=F)

# subsample
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




#------ save/load ----
save(list=c('res_04','res_424','res_2472'), file = res_file)

# To skip the analysis, load the saved results instead
# load(res_file)



#------ TABLE of differences ---------------------------

#--------- t0 t4 ----

results = res_04$g.diff
results[is.na(results)] = 0

coef = res_04$base_coef
delta = res_04$delta_coef

# Intercepts
i_diff = results[,2]
table(i_diff)

int_t0 = coef[i_diff==1,1]
int_d = delta[i_diff==1,1]
# table(names(int_t0) == names(int_d))
int_t4 = int_t0 + int_d

tab_int_04 = round(cbind('t_0'=int_t0, 't_4'=int_t4, 'delta'=int_d,
                         'exp(t_0)'=exp(int_t0), 'exp(t_4)'=exp(int_t4), 'exp(delta)'=exp(int_d)),3)
# tab_mu_04  = round(exp(tab_int_04),3)


# Associations
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

# Intercepts
i_diff = results[,2]
table(i_diff)

int_t0 = coef[i_diff==1,1]
int_d = delta[i_diff==1,1]
# table(names(int_t0) == names(int_d))
int_t4 = int_t0 + int_d

tab_int_424 = round(cbind('t_4'=int_t0, 't_24'=int_t4, 'delta'=int_d,
                          'exp(t_4)'=exp(int_t0), 'exp(t_24)'=exp(int_t4), 'exp(delta)'=exp(int_d)),3)
# tab_mu_424  = round(exp(tab_int_424),3)


# Associations
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

# Intercepts
i_diff = results[,2]
table(i_diff)

int_t0 = coef[i_diff==1,1]
int_d = delta[i_diff==1,1]
# table(names(int_t0) == names(int_d))
int_t4 = int_t0 + int_d

tab_int_2472 = round(cbind('t_24'=int_t0, 't_72'=int_t4, 'delta'=int_d,
                           'exp(t_24)'=exp(int_t0), 'exp(t_72)'=exp(int_t4), 'exp(delta)'=exp(int_d)),3)
# tab_mu_2472  = round(exp(tab_int_2472),3)


# Associations
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
