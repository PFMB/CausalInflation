rm(list = ls())
library(parallel)
library(BaBooN) # Rubin's Rule

# insert path
setwd("/cluster/home/phibauma/CausalInflation/")

# load 5 imputed data.frames that are analyzed
load("causalinfl.RData")

# initiate cluster
no_cores <- 5
cl <- makeCluster(no_cores, outfile = "")
clusterSetRNGStream(cl = cl, iseed = 1)
clusterEvalQ(cl,library(ltmle))

estimation_ltmle <- function(dat, path = path) {
  
  # load all customized (besides defaults) learner/screeners
  source("LearnerLibrary.R")
  load("SelectedLearners.RData")
  # gbm failes frequently and was not selected often during previous analyses
  SL.Est_Data <- SL.Est_Data[-c(10,21,32,43,54)] # exclude GBM
  
  # load nodes, g-/Q-formulas and interventions
  load("NodesFormInterv.RData")
  nod <- ltmle_prep$nodes
  int <- ltmle_prep$invterventions
  
  # -------- ScreenLearn ----- Static -------- #
  
  ScreenLearnSta <-  tryCatch({ltmle(dat,
                                Lnodes = nod$L, Anodes = nod$A, Ynodes = nod$Y,
                                estimate.time = FALSE,
                                Yrange = c(-3, 86),
                                variance.method = "tmle",
                                abar = list(treament = int$stat1, control = int$stat0),
                                SL.library = SL.Est_Data
  )}, error = function(cond) return(NA))
  
  # -------- ScreenLearn ----- Dynamic -------- #
  
  ScreenLearnDyn <- tryCatch({ltmle(dat,
                              Lnodes = nod$L, Anodes = nod$A, Ynodes = nod$Y,
                              estimate.time = FALSE,
                              Yrange = c(-3, 86),
                              variance.method = "tmle",
                              abar = list(treament = int$dyn_intv, control = int$stat0),
                              SL.library = SL.Est_Data
  )}, error = function(cond) return(NA))
  
  # -------- EconDAG ----- Static -------- #
  
  EconDAGSta <- tryCatch({ltmle(dat,
                               Lnodes = nod$L, Anodes = nod$A, Ynodes = nod$Y,
                               Qform = ltmle_prep$form_econ$Q_econ,
                               gform = ltmle_prep$form_econ$g_econ,
                               estimate.time = FALSE,
                               Yrange = c(-3, 86),
                               variance.method = "tmle",
                               abar = list(treament = int$stat1, control = int$stat0),
                               SL.library = SL.Est_Theory
  )}, error = function(cond) return(NA))
  
  # -------- EconDAG ----- Dynamic -------- #
  
  # estimation
  EconDAGDyn <- tryCatch({ltmle(dat,
                              Lnodes = nod$L, Anodes = nod$A, Ynodes = nod$Y,
                              Qform = ltmle_prep$form_econ$Q_econ,
                              gform = ltmle_prep$form_econ$g_econ,
                              estimate.time = FALSE,
                              Yrange = c(-3, 86),
                              variance.method = "tmle",
                              abar = list(treament = int$dyn_intv, control = int$stat0),
                              SL.library = SL.Est_Theory
  )}, error = function(cond) return(NA))
  
  # -------- PlainDAG ----- Static -------- #
  
  # estimation
  PlainDAGSta <- tryCatch({ltmle(dat,
                               Lnodes = nod$L, Anodes = nod$A, Ynodes = nod$Y,
                               Qform = ltmle_prep$form_base$Q_base,
                               gform = ltmle_prep$form_base$g_base,
                               estimate.time = FALSE,
                               Yrange = c(-3, 86),
                               variance.method = "tmle",
                               abar = list(treament = int$stat1, control = int$stat0),
                               SL.library =  SL.Est_Theory
  )}, error = function(cond) return(NA))
  
  # -------- PlainDAG ----- Dynamic -------- #
  
  # estimation
  PlainDAGDyn <- tryCatch({ltmle(dat,
                              Lnodes = nod$L, Anodes = nod$A, Ynodes = nod$Y,
                              Qform = ltmle_prep$form_base$Q_base,
                              gform = ltmle_prep$form_base$g_base,
                              estimate.time = FALSE,
                              Yrange = c(-3, 86),
                              variance.method = "tmle",
                              abar = list(treament = int$dyn_intv, control = int$stat0),
                              SL.library =  SL.Est_Theory
  )}, error = function(cond) return(NA))
  
  list(ScreenLearnSta = ScreenLearnSta, ScreenLearnDyn = ScreenLearnDyn,
       EconDAGSta = EconDAGSta, EconDAGDyn = EconDAGDyn,
       PlainDAGSta = PlainDAGSta, PlainDAGDyn = PlainDAGDyn)
}

# run estimation
start_time <- Sys.time()
res <- parLapply(cl, infl, estimation_ltmle)
end_time <- Sys.time()
stopCluster(cl)

attributes(res)$sessInfo <- sessionInfo()
attributes(res)$time <- end_time - start_time

cat("Elapsed time:", attributes(res)$time, "\n")

prep_res <- function(res){
  
  # extract results from all ltmle outputs and put in result matrix
  est <- unlist(res["estimate",])
  std <- unlist(res["std.dev",])
  CI <- do.call("rbind",res["CI",])
  pvalues <- unlist(res["pvalue",])
  cbind(est,std,CI,pvalues)
  
}

# LTMLE: combine/average the results from above with Rubins Rule and p2s
ltmle_res <- lapply(res, function(res_s, est_meth = "tmle") {
  sapply(res_s, function(res_x) unclass(summary(res_x, estimator = est_meth)$effect.measures$ATE))
})
ltmle_res <- lapply(ltmle_res, prep_res)
ltmle_res <- do.call("cbind",ltmle_res)
(ltmle_res <- sapply(rownames(ltmle_res), function(row){
  MI.inference(thetahat = ltmle_res[row, colnames(ltmle_res) == "est"], 
               varhat.thetahat = ltmle_res[row, colnames(ltmle_res) == "std"]^2)
}))

# IPTW: combine/average the results from above with Rubins Rule and p2s
iptw_res <- lapply(res, function(res_s, est_meth = "iptw") {
  sapply(res_s, function(res_x) unclass(summary(res_x, estimator = est_meth)$effect.measures$ATE))
})
iptw_res <- lapply(iptw_res, prep_res)
iptw_res <- do.call("cbind",iptw_res)
(iptw_res <- sapply(rownames(iptw_res), function(row){
  MI.inference(thetahat = iptw_res[row, colnames(iptw_res) == "est"],
               varhat.thetahat = iptw_res[row, colnames(iptw_res) == "std"]^2)
}))


