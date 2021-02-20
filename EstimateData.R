# Please install the following packages prior to running the file: 
# ltmle, SuperLearner, vcd, arm, rpart, nnet, glmnet, stringr, magrittr, randomForest,
# earth, gbm, gam, mgcv, reshape2, dplyr, data.table, polspline
# plyr, tibble, scales, BaBooN, simcausal (currently only on CRAN archive)

# load in memory and attach to serachpath
rm(list = ls())
library(parallel)
library(BaBooN) # for Rubin's Rules after Imputation
library(ltmle)
library(vcd)
set.seed(1)
path <- "/cluster/home/phibauma/CausalInflation"

# insert working directory
setwd(path)

# specify here how many cores are available for parallel computation (should be 5 here)
n_cluster <- 5

# load 5 imputed data.frames that are analyzed
load("causalinfl_revised.RData")

# extract learner weights
source("WeightsSummary.R")

# formula to extract ATE
get_ATE <- function(out, est = "tmle") unclass(summary(out, estimator = est))$effect.measures$ATE

# Learner Sets: SL.Est_Theory, SL.Est_Data
load("SelectedLearner.RData")

# initiate cluster
cl <- makeCluster(n_cluster, outfile = "")
clusterSetRNGStream(cl = cl, iseed = 1)
clusterEvalQ(cl, library(ltmle))

#SL.Est_Theory <- SL.Est_Data <- c("SL.glm","SL.mean")

clusterExport(cl = cl, list(
  "learner_weights_summary_g",
  "learner_weights_summary_Q", "get_ATE","SL.Est_Theory","SL.Est_Data"
))

est_ScreenLearnSta <- function(dat) {
  
  source("LearnerLibrary.R")
  
  # load nodes, g-/Q-formulas and interventions
  load("NodesFormInterv.RData")
  nod <- ltmle_prep$nodes
  int <- ltmle_prep$invterventions
  
  ids <- rownames(dat)
  int$dyn_intv <- int$dyn_intv[rownames(int$dyn_intv) %in% ids,]
  int$stat0 <- int$stat0[1:nrow(dat),]
  int$stat1 <- int$stat1[1:nrow(dat),]

  cat("# -------- ScreenLearn ----- Static -------- # \n")

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
  
  list(est_out = ests, weights_out = learn_weights, cc = try(cc_trunc(ScreenLearnSta), silent = TRUE))
  
}

res <- parLapply(cl, infl$all, est_ScreenLearnSta)
attr(res,"info") <- sessionInfo()
attr(res,"seed") <- .Random.seed
attr(res,"time") <- Sys.time()
saveRDS(res, file = paste0(path,"/results/ScreenLearnSta_all"))
rm(res)

res <- parLapply(cl, infl$low, est_ScreenLearnSta)
saveRDS(res, file = paste0(path,"/results/ScreenLearnSta_low"))
rm(res)
res <- parLapply(cl, infl$high, est_ScreenLearnSta)
saveRDS(res, file = paste0(path,"/results/ScreenLearnSta_high"))
rm(res)

est_ScreenLearnDyn <- function(dat) {
  
  source("LearnerLibrary.R")
  
  # load nodes, g-/Q-formulas and interventions
  load("NodesFormInterv.RData")
  nod <- ltmle_prep$nodes
  int <- ltmle_prep$invterventions
  
  ids <- rownames(dat)
  int$dyn_intv <- int$dyn_intv[rownames(int$dyn_intv) %in% ids,]
  int$stat0 <- int$stat0[1:nrow(dat),]
  int$stat1 <- int$stat1[1:nrow(dat),]

  cat("# -------- ScreenLearn ----- Dynamic -------- # \n")

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
  list(est_out = ests, weights_out = learn_weights, cc = try(cc_trunc(ScreenLearnDyn), silent = TRUE))
}

res <- parLapply(cl, infl$all, est_ScreenLearnDyn)
saveRDS(res, file = paste0(path,"/results/ScreenLearnDyn_all"))
rm(res)
res <- parLapply(cl, infl$low, est_ScreenLearnDyn)
saveRDS(res, file = paste0(path,"/results/ScreenLearnDyn_low"))
rm(res)
res <- parLapply(cl, infl$high, est_ScreenLearnDyn)
saveRDS(res, file = paste0(path,"/results/ScreenLearnDyn_high"))
rm(res)

est_EconDAGSta <- function(dat) {
  
  source("LearnerLibrary.R")

  # load nodes, g-/Q-formulas and interventions
  load("NodesFormInterv.RData")
  nod <- ltmle_prep$nodes
  int <- ltmle_prep$invterventions
  
  ids <- rownames(dat)
  int$dyn_intv <- int$dyn_intv[rownames(int$dyn_intv) %in% ids,]
  int$stat0 <- int$stat0[1:nrow(dat),]
  int$stat1 <- int$stat1[1:nrow(dat),]

  cat("# -------- EconDAG ----- Static -------- # \n")

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
  
  list(est_out = ests, weights_out = learn_weights, cc = try(cc_trunc(EconDAGSta), silent = TRUE))
  
}

res <- parLapply(cl, infl$all, est_EconDAGSta)
saveRDS(res, file = paste0(path,"/results/EconDAGSta_all"))
rm(res)
res <- parLapply(cl, infl$low, est_EconDAGSta)
saveRDS(res, file = paste0(path,"/results/EconDAGSta_low"))
rm(res)
res <- parLapply(cl, infl$high, est_EconDAGSta)
saveRDS(res, file = paste0(path,"/results/EconDAGSta_high"))
rm(res)

est_EconDAGDyn <- function(dat) {
  
  source("LearnerLibrary.R")
  
  # load nodes, g-/Q-formulas and interventions
  load("NodesFormInterv.RData")
  nod <- ltmle_prep$nodes
  int <- ltmle_prep$invterventions
  
  ids <- rownames(dat)
  int$dyn_intv <- int$dyn_intv[rownames(int$dyn_intv) %in% ids,]
  int$stat0 <- int$stat0[1:nrow(dat),]
  int$stat1 <- int$stat1[1:nrow(dat),]
  
  cat("# -------- EconDAG ----- Dynamic -------- # \n")

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
  
  list(est_out = ests, weights_out = learn_weights, cc = try(cc_trunc(EconDAGDyn), silent = TRUE))
}

res <- parLapply(cl, infl$all, est_EconDAGDyn)
saveRDS(res, file = paste0(path,"/results/EconDAGDyn_all"))
rm(res)
res <- parLapply(cl, infl$low, est_EconDAGDyn)
saveRDS(res, file = paste0(path,"/results/EconDAGDyn_low"))
rm(res)
res <- parLapply(cl, infl$high, est_EconDAGDyn)
saveRDS(res, file = paste0(path,"/results/EconDAGDyn_high"))
rm(res)

est_PlainDAGSta <- function(dat) {
  
  source("LearnerLibrary.R")
  
  # load nodes, g-/Q-formulas and interventions
  load("NodesFormInterv.RData")
  nod <- ltmle_prep$nodes
  int <- ltmle_prep$invterventions
  
  ids <- rownames(dat)
  int$dyn_intv <- int$dyn_intv[rownames(int$dyn_intv) %in% ids,]
  int$stat0 <- int$stat0[1:nrow(dat),]
  int$stat1 <- int$stat1[1:nrow(dat),]

  cat("# -------- PlainDAG ----- Static -------- # \n")

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
  learn_weights <- list(Qweights = Q_mean, gweights = g_mean, cc = try(cc_trunc(PlainDAGSta), silent = TRUE))
  
  # extract ATEs from estimation
  ests <- list(
    ltmle = try(get_ATE(PlainDAGSta), silent = TRUE),
    iptw = try(get_ATE(PlainDAGSta, est = "iptw"), silent = TRUE)
  )
  
  list(est_out = ests, weights_out = learn_weights, cc = try(cc_trunc(PlainDAGSta), silent = TRUE))
}

res <- parLapply(cl, infl$all, est_PlainDAGSta)
saveRDS(res, file = paste0(path,"/results/PlainDAGSta_all"))
rm(res)
res <- parLapply(cl, infl$low, est_PlainDAGSta)
saveRDS(res, file = paste0(path,"/results/PlainDAGSta_low"))
rm(res)
res <- parLapply(cl, infl$high, est_PlainDAGSta)
saveRDS(res, file = paste0(path,"/results/PlainDAGSta_high"))
rm(res)

est_PlainDAGDyn <- function(dat) {
  
  source("LearnerLibrary.R")
  
  # load nodes, g-/Q-formulas and interventions
  load("NodesFormInterv.RData")
  nod <- ltmle_prep$nodes
  int <- ltmle_prep$invterventions
  
  ids <- rownames(dat)
  int$dyn_intv <- int$dyn_intv[rownames(int$dyn_intv) %in% ids,]
  int$stat0 <- int$stat0[1:nrow(dat),]
  int$stat1 <- int$stat1[1:nrow(dat),]

  cat("# -------- PlainDAG ----- Dynamic -------- # \n")

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
  list(est_out = ests, weights_out = learn_weights, cc = try(cc_trunc(PlainDAGDyn), silent = TRUE))
  
}

res <- parLapply(cl, infl$all, est_PlainDAGDyn)
saveRDS(res, file = paste0(path,"/results/PlainDAGDyn_all"))
rm(res)
res <- parLapply(cl, infl$low, est_PlainDAGDyn)
saveRDS(res, file = paste0(path,"/results/PlainDAGDyn_low"))
rm(res)
res <- parLapply(cl, infl$high, est_PlainDAGDyn)
saveRDS(res, file = paste0(path,"/results/PlainDAGDyn_high"))
rm(res)

stopCluster(cl)

# # retrieve results
# prep_res <- function(res) {
# 
#   # extract results from all ltmle outputs and put in result matrix
#   est <- unlist(res["estimate", ])
#   std <- unlist(res["std.dev", ])
#   CI <- do.call("rbind", res["CI", ])
#   pvalues <- unlist(res["pvalue", ])
#   cbind(est, std, CI, pvalues)
# }
# 
# # LTMLE: combine/average the results from above with Rubins Rule and p2s
# ltmle_res <- lapply(res_all, function(analyzed_set) {
#   sapply(analyzed_set, function(est_strategy) est_strategy$est_out$ltmle)
# })
# ltmle_res <- lapply(ltmle_res, prep_res)
# ltmle_res <- do.call("cbind", ltmle_res)
# (ltmle_res <- sapply(rownames(ltmle_res), function(est_strategy) {
#   MI.inference(
#     thetahat = ltmle_res[est_strategy, colnames(ltmle_res) == "est"],
#     varhat.thetahat = ltmle_res[est_strategy, colnames(ltmle_res) == "std"]^2
#   )
# }))
