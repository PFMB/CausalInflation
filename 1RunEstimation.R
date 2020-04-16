rm(list = ls())
library(parallel)
set.seed(1)

# insert path
setwd("/opt/rdata/phibauma/CausalInflation/GitHubVersion")

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
  
  # time in milliseconds for saving
  options(digits.secs = 6)
  
  # load nodes, g-/Q-formulas and interventions
  load("NodesFormInterv.RData")
  nod <- ltmle_prep$nodes
  int <- ltmle_prep$invterventions
  
  # -------- FULL DATA----- STATIC INTERVENTION -------- #
  
  data_stat <-  tryCatch({ltmle(dat,
                                Lnodes = nod$L, Anodes = nod$A, Ynodes = nod$Y,
                                estimate.time = FALSE,
                                Yrange = c(-3, 86),
                                variance.method = "tmle",
                                abar = list(treament = int$stat1, control = int$stat0),
                                SL.library = SL.Est_Data
  )}, error = function(cond) return(NA))
  
  cat("data_stat estimation is over.")
  
  # -------- FULL DATA ----- DYNAMIC INTERVENTION -------- #
  
  data_dyn <- tryCatch({ltmle(dat,
                              Lnodes = nod$L, Anodes = nod$A, Ynodes = nod$Y,
                              estimate.time = FALSE,
                              Yrange = c(-3, 86),
                              variance.method = "tmle",
                              abar = list(treament = int$dyn_intv, control = int$stat0),
                              SL.library = SL.Est_Data
  )}, error = function(cond) return(NA))
  
  cat("data_dyn estimation is over.")
  
  # -------- ECON ----- STATIC INTERVENTION -------- #
  
  econ_stat <- tryCatch({ltmle(dat,
                               Lnodes = nod$L, Anodes = nod$A, Ynodes = nod$Y,
                               Qform = ltmle_prep$form_econ$Q_econ,
                               gform = ltmle_prep$form_econ$g_econ,
                               estimate.time = FALSE,
                               Yrange = c(-3, 86),
                               variance.method = "tmle",
                               abar = list(treament = int$stat1, control = int$stat0),
                               SL.library = SL.Est_Theory
  )}, error = function(cond) return(NA))
  
  cat("econ_stat estimation is over.\n")
  
  # -------- ECON ----- DYNAMIC INTERVENTION -------- #
  
  # estimation
  econ_dyn <- tryCatch({ltmle(dat,
                              Lnodes = nod$L, Anodes = nod$A, Ynodes = nod$Y,
                              Qform = ltmle_prep$form_econ$Q_econ,
                              gform = ltmle_prep$form_econ$g_econ,
                              estimate.time = FALSE,
                              Yrange = c(-3, 86),
                              variance.method = "tmle",
                              abar = list(treament = int$dyn_intv, control = int$stat0),
                              SL.library = SL.Est_Theory
  )}, error = function(cond) return(NA))
  
  cat("econ_dyn estimation is over.\n")
  
  # -------- BASE ----- STATIC INTERVENTION -------- #
  
  # estimation
  base_stat <- tryCatch({ltmle(dat,
                               Lnodes = nod$L, Anodes = nod$A, Ynodes = nod$Y,
                               Qform = ltmle_prep$form_base$Q_base,
                               gform = ltmle_prep$form_base$g_base,
                               estimate.time = FALSE,
                               Yrange = c(-3, 86),
                               variance.method = "tmle",
                               abar = list(treament = int$stat1, control = int$stat0),
                               SL.library =  SL.Est_Theory
  )}, error = function(cond) return(NA))
  
  cat("base_stat estimation is over.\n")
  
  # -------- BASE ----- DYNAMIC INTERVENTION -------- #
  
  # estimation
  base_dyn <- tryCatch({ltmle(dat,
                              Lnodes = nod$L, Anodes = nod$A, Ynodes = nod$Y,
                              Qform = ltmle_prep$form_base$Q_base,
                              gform = ltmle_prep$form_base$g_base,
                              estimate.time = FALSE,
                              Yrange = c(-3, 86),
                              variance.method = "tmle",
                              abar = list(treament = int$dyn_intv, control = int$stat0),
                              SL.library =  SL.Est_Theory
  )}, error = function(cond) return(NA))
  
  cat("base_dyn estimation is over.\n")
  
  list(
    data_stat = data_stat, data_dyn = data_dyn,
    econ_stat = econ_stat, econ_dyn = econ_dyn,
    base_stat = base_stat, base_dyn = base_dyn
  )
}

start_time <- Sys.time()
res <- parLapply(cl, infl, estimation_ltmle)
end_time <- Sys.time()
stopCluster(cl)

attributes(res)$sessInfo <- sessionInfo()
attributes(res)$seed <- .Random.seed
attributes(res)$time <- end_time - start_time


