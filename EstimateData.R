################################################################################
# Philipp Baumann, Michael Schomaker, and Enzo Rossi
# Working Paper (under submission)
# Title:
# Estimating the Effect of Central Bank Independence on Inflation Using 
# Longitudinal Targeted Maximum Likelihood Estimation
################################################################################

################################################################################
# required packages
rm(list = ls())
library(parallel)
library(BaBooN) # for Rubin's Rules after Imputation
library(ltmle)
library(stringr)
library(vcd)
library(magrittr)
set.seed(1)

# insert path
setwd("/cluster/home/phibauma/CausalInflation")

# specify here how many cores are available for parallel computation (should be 5 here)
n.cluster <- 5

# load 5 imputed data.frames that are analyzed
load("causalinfl.RData")

# initiate cluster
cl <- makeCluster(n.cluster, outfile = "")
clusterSetRNGStream(cl = cl, iseed = 1)
clusterEvalQ(cl, library(ltmle))

# helper functions
source("WeightsSummary.R")

# formula to extract ATE
get_ATE <- function(out, est = "tmle") unclass(summary(out, estimator = est))$effect.measures$ATE

clusterExport(cl = cl, list(
  "learner_weights_summary_g",
  "learner_weights_summary_Q", "get_ATE"
))

estimation_ltmle <- function(dat, path = path) {
  
  source("LearnerLibrary.R")
  load("SelectedLearners.RData")

  # gbm failes frequently, takes a lot of time and was not selected often during previous analyses
  SL.Est_Data <- SL.Est_Data[-c(10,21,32,43,54)] # exclude GBM
  attr(SL.Est_Data,"return.fit") <- TRUE # to access learner weights

  # load nodes, g-/Q-formulas and interventions
  load("NodesFormInterv.RData")
  nod <- ltmle_prep$nodes
  int <- ltmle_prep$invterventions

  # -------- ScreenLearn ----- Static -------- #

  ScreenLearnSta <- tryCatch(
    {
      ltmle(dat,
        Lnodes = nod$L, Anodes = nod$A, Ynodes = nod$Y,
        estimate.time = FALSE,
        Yrange = c(-3, 86),
        variance.method = "tmle",
        abar = list(treament = int$stat1, control = int$stat0),
        SL.library = SL.Est_Data
      )
    },
    error = function(cond) {
      return(NA)
    }
  )
  
  # extract learner weights from estimation
  Q_mean <- try(learner_weights_summary_Q(ScreenLearnSta), silent = TRUE)
  g_mean <- try(learner_weights_summary_g(ScreenLearnSta), silent = TRUE)
  learn_weights <- list(Qweights = Q_mean, gweights = g_mean)
  
  # extract ATEs from estimation
  ests <- list(
    ltmle = try(get_ATE(ScreenLearnSta), silent = TRUE),
    iptw = try(get_ATE(ScreenLearnSta, est = "iptw"), silent = TRUE)
  )
  ScreenLearnSta <- list(est_out = ests, weights_out = learn_weights)

  # -------- ScreenLearn ----- Dynamic -------- #

  ScreenLearnDyn <- tryCatch(
    {
      ltmle(dat,
        Lnodes = nod$L, Anodes = nod$A, Ynodes = nod$Y,
        estimate.time = FALSE,
        Yrange = c(-3, 86),
        variance.method = "tmle",
        abar = list(treament = int$dyn_intv, control = int$stat0),
        SL.library = SL.Est_Data
      )
    },
    error = function(cond) {
      return(NA)
    }
  )
  
  # extract learner weights from estimation
  Q_mean <- try(learner_weights_summary_Q(ScreenLearnDyn), silent = TRUE)
  g_mean <- try(learner_weights_summary_g(ScreenLearnDyn), silent = TRUE)
  learn_weights <- list(Qweights = Q_mean, gweights = g_mean)
  
  # extract ATEs from estimation
  ests <- list(
    ltmle = try(get_ATE(ScreenLearnDyn), silent = TRUE),
    iptw = try(get_ATE(ScreenLearnDyn, est = "iptw"), silent = TRUE)
  )
  ScreenLearnDyn <- list(est_out = ests, weights_out = learn_weights)

  # -------- EconDAG ----- Static -------- #

  EconDAGSta <- tryCatch(
    {
      ltmle(dat,
        Lnodes = nod$L, Anodes = nod$A, Ynodes = nod$Y,
        Qform = ltmle_prep$form_econ$Q_econ,
        gform = ltmle_prep$form_econ$g_econ,
        estimate.time = FALSE,
        Yrange = c(-3, 86),
        variance.method = "tmle",
        abar = list(treament = int$stat1, control = int$stat0),
        SL.library = SL.Est_Theory
      )
    },
    error = function(cond) {
      return(NA)
    }
  )
  
  # extract learner weights from estimation
  Q_mean <- try(learner_weights_summary_Q(EconDAGSta), silent = TRUE)
  g_mean <- try(learner_weights_summary_g(EconDAGSta), silent = TRUE)
  learn_weights <- list(Qweights = Q_mean, gweights = g_mean)
  
  # extract ATEs from estimation
  ests <- list(
    ltmle = try(get_ATE(EconDAGSta), silent = TRUE),
    iptw = try(get_ATE(EconDAGSta, est = "iptw"), silent = TRUE)
  )
  EconDAGSta <- list(est_out = ests, weights_out = learn_weights)

  # -------- EconDAG ----- Dynamic -------- #

  # estimation
  EconDAGDyn <- tryCatch(
    {
      ltmle(dat,
        Lnodes = nod$L, Anodes = nod$A, Ynodes = nod$Y,
        Qform = ltmle_prep$form_econ$Q_econ,
        gform = ltmle_prep$form_econ$g_econ,
        estimate.time = FALSE,
        Yrange = c(-3, 86),
        variance.method = "tmle",
        abar = list(treament = int$dyn_intv, control = int$stat0),
        SL.library = SL.Est_Theory
      )
    },
    error = function(cond) {
      return(NA)
    }
  )
  
  # extract learner weights from estimation
  Q_mean <- try(learner_weights_summary_Q(EconDAGDyn), silent = TRUE)
  g_mean <- try(learner_weights_summary_g(EconDAGDyn), silent = TRUE)
  learn_weights <- list(Qweights = Q_mean, gweights = g_mean)
  
  # extract ATEs from estimation
  ests <- list(
    ltmle = try(get_ATE(EconDAGDyn), silent = TRUE),
    iptw = try(get_ATE(EconDAGDyn, est = "iptw"), silent = TRUE)
  )
  EconDAGDyn <- list(est_out = ests, weights_out = learn_weights)

  # -------- PlainDAG ----- Static -------- #

  # estimation
  PlainDAGSta <- tryCatch(
    {
      ltmle(dat,
        Lnodes = nod$L, Anodes = nod$A, Ynodes = nod$Y,
        Qform = ltmle_prep$form_base$Q_base,
        gform = ltmle_prep$form_base$g_base,
        estimate.time = FALSE,
        Yrange = c(-3, 86),
        variance.method = "tmle",
        abar = list(treament = int$stat1, control = int$stat0),
        SL.library = SL.Est_Theory
      )
    },
    error = function(cond) {
      return(NA)
    }
  )
  
  # extract learner weights from estimation
  Q_mean <- try(learner_weights_summary_Q(PlainDAGSta), silent = TRUE)
  g_mean <- try(learner_weights_summary_g(PlainDAGSta), silent = TRUE)
  learn_weights <- list(Qweights = Q_mean, gweights = g_mean)
  
  # extract ATEs from estimation
  ests <- list(
    ltmle = try(get_ATE(PlainDAGSta), silent = TRUE),
    iptw = try(get_ATE(PlainDAGSta, est = "iptw"), silent = TRUE)
  )
  PlainDAGSta <- list(est_out = ests, weights_out = learn_weights)

  # -------- PlainDAG ----- Dynamic -------- #

  # estimation
  PlainDAGDyn <- tryCatch(
    {
      ltmle(dat,
        Lnodes = nod$L, Anodes = nod$A, Ynodes = nod$Y,
        Qform = ltmle_prep$form_base$Q_base,
        gform = ltmle_prep$form_base$g_base,
        estimate.time = FALSE,
        Yrange = c(-3, 86),
        variance.method = "tmle",
        abar = list(treament = int$dyn_intv, control = int$stat0),
        SL.library = SL.Est_Theory
      )
    },
    error = function(cond) {
      return(NA)
    }
  )
  
  # extract learner weights from estimation
  Q_mean <- try(learner_weights_summary_Q(PlainDAGDyn), silent = TRUE)
  g_mean <- try(learner_weights_summary_g(PlainDAGDyn), silent = TRUE)
  learn_weights <- list(Qweights = Q_mean, gweights = g_mean)
  
  # extract ATEs from estimation
  ests <- list(
    ltmle = try(get_ATE(PlainDAGDyn), silent = TRUE),
    iptw = try(get_ATE(PlainDAGDyn, est = "iptw"), silent = TRUE)
  )
  PlainDAGDyn <- list(est_out = ests, weights_out = learn_weights)

  list(
    ScreenLearnSta = ScreenLearnSta, ScreenLearnDyn = ScreenLearnDyn,
    EconDAGSta = EconDAGSta, EconDAGDyn = EconDAGDyn,
    PlainDAGSta = PlainDAGSta, PlainDAGDyn = PlainDAGDyn
  )
}

# run estimation
t_ime <- system.time({
res <- parLapply(cl, infl, estimation_ltmle)
})
stopCluster(cl)

# save results
(attributes(res)$time <- t_ime)
(attributes(res)$sessinfo <- sessionInfo())
(attributes(res)$seed <- .Random.seed)
saveRDS(res, file = "EstResults.RDS")

# retrieve results
prep_res <- function(res) {

  # extract results from all ltmle outputs and put in result matrix
  est <- unlist(res["estimate", ])
  std <- unlist(res["std.dev", ])
  CI <- do.call("rbind", res["CI", ])
  pvalues <- unlist(res["pvalue", ])
  cbind(est, std, CI, pvalues)
}

# LTMLE: combine/average the results from above with Rubins Rule and p2s
ltmle_res <- lapply(res, function(analyzed_set) {
  sapply(analyzed_set, function(est_strategy) est_strategy$est_out$ltmle)
})
ltmle_res <- lapply(ltmle_res, prep_res)
ltmle_res <- do.call("cbind", ltmle_res)
(ltmle_res <- sapply(rownames(ltmle_res), function(est_strategy) {
  MI.inference(
    thetahat = ltmle_res[est_strategy, colnames(ltmle_res) == "est"],
    varhat.thetahat = ltmle_res[est_strategy, colnames(ltmle_res) == "std"]^2
  )
}))
