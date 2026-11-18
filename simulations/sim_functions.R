
## Functions for the simulation study ##
#
# - simdata_tpgm(): simulates data from two Truncated Poisson Graphical Models
# - OFDR(), power.1stage(): summary measures for stage 1 (screening)
# - fwer.2stage(), power.2stage(): summary measures for stage 2 (confirmatory)
#
# The two-stage test itself is performed with DUE2::DUE2().


#---- Libraries -------------------------------------------------------------

library(ModPGMInference)   # ModPGMSampler()
library(spsUtil)           # quiet()


#---- Data simulation -------------------------------------------------------

# Simulates n observations from a TPGM with parameters (th1, phi1) (group 0)
# and n observations from a TPGM with parameters (th2, phi2) (group 1).
#
# param: list with elements
#   th1, th2   : node parameters (vectors of length p)
#   phi1, phi2 : edge parameters (symmetric p x p matrices)
#   R          : truncation values (vector of length p)

simdata_tpgm = function(n, param){

  # Param
  R = param$R
  p = length(param$th1)
  th1 = param$th1
  th2 = param$th2
  phi1 = param$phi1
  phi2 = param$phi2

  # Sample
  data1 = data.frame(quiet(ModPGMSampler(psi=th1, true_graph=phi1, model="TPGM", D=R, nSample=n, burn_in=5000, thin=250)))
  data2 = data.frame(quiet(ModPGMSampler(psi=th2, true_graph=phi2, model="TPGM", D=R, nSample=n, burn_in=5000, thin=250)))
  data  = rbind(data1, data2)
  colnames(data) = 1:p

  return(data)
}



#---- Summary measures ------------------------------------------------------
#
# 'result' / 'results' is a list (one element per simulation) of the 'g.diff'
# matrices returned by DUE2().

## STAGE 1

# OFDR: overall FDR of stage 1, i.e. (# false rejected nodes)/(# rejected nodes),
#       averaged over the simulations
# power.1stage: for each truly differential node, proportion of simulations
#               in which it is rejected

OFDR = function(result, true.nodes){

  simFDR = sapply(result, function(sim){

    rownames(sim) = colnames(sim)[-c(1:2)]

    totalRej = sum(sim[,'LRT'])
    falseRej = sum(sim[!(rownames(sim) %in% true.nodes),'LRT'])

    # single simulation
    if(totalRej > 0){
      fdr = falseRej/totalRej
    }else{
      fdr = 0
    }
    return(fdr)
  })

  ofdr = mean(simFDR)
  return(ofdr)
}


power.1stage = function(result, true.nodes){
  df.res = sapply(result, function(sim) sim[,'LRT'])
  rowm = rowMeans(df.res)
  names(rowm) = colnames(result[[1]])[-c(1:2)]
  rowm[true.nodes]
}


## STAGE 2

# fwer.2stage: node-wise FWER, i.e. for each node, proportion of simulations
#              with at least one false rejection
# power.2stage: for each true difference, proportion of simulations in which
#               it is rejected
#
# true.param: p x (p+1) matrix of 0/1; first column = node (mean) differences,
#             other columns = edge (association) differences

fwer.2stage = function(results, true.param){
  res = rowMeans(sapply(1:length(results), function(sim){
    r1 = results[[sim]][,-1]
    g.mat = ifelse(is.na(r1),0,r1)
    g.mat[which(true.param==1, arr.ind=T)] = 0
    out = ifelse(rowSums(g.mat)>0,1,0)
    return(out)
  }))
  return(res)
}


power.2stage = function(results, true.param){

  res = lapply(1:length(results), function(sim){
    df = results[[sim]][,-1]
    ifelse(is.na(df),0,df)
  })

  sumres = Reduce(`+`, res)

  out = (sumres * true.param)/length(results)
  return(out)
}
