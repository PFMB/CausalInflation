# Please install the following packages prior to running the file: 
# ltmle, SuperLearner, vcd, arm, rpart, nnet, glmnet, stringr, magrittr, randomForest,
# earth, gbm, gam, mgcv, reshape2, dplyr, data.table, polspline
# plyr, tibble, scales, BaBooN, simcausal (currently only on CRAN archive)

# load in memory and attach to serachpath
rm(list = ls())
library(parallel)
set.seed(1)
path <- "/cluster/home/phibauma/CausalInflation"
#path <- "/Users/flipst3r/RStHomeDir/GitHub/CausalInflation/"

# specify here how many cores are available for parallel computation (should be 5 here)
n_cluster <- 5

# load 5 imputed data.frames that are analyzed
load(paste0(path,"causalinfl_revised.RData"))

# extract learner weights
source(paste0(path,"WeightsSummary.R"))

# formula to extract ATE
get_ATE <- function(out, est = "tmle") unclass(summary(out, estimator = est))$effect.measures$ATE

# Learner Sets: SL.Est_Theory, SL.Est_Data
load(paste0(path,"SelectedLearner.RData"))

# initiate cluster
cl <- makeCluster(n_cluster, outfile = "")
clusterSetRNGStream(cl = cl, iseed = 1)
clusterEvalQ(cl, library(ltmle))

#SL.Est_Theory <- SL.Est_Data <- c("SL.glm","SL.mean") # for fast test runs
SL.Est_Theory <- SL.Est_Theory[-10] # gbm consumes too much memory and does not help much either

clusterExport(cl = cl, list(
  "learner_weights_summary_g",
  "learner_weights_summary_Q", "get_ATE","SL.Est_Theory","SL.Est_Data","path"
))

est_ScreenLearnSta <- function(dat) {

  source(paste0(path,"LearnerLibrary.R"))

  # load nodes, g-/Q-formulas and interventions
  load(paste0(path,"NodesFormInterv.RData"))
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
  cat("Object size:", format(object.size(ScreenLearnSta), units = "auto"),"\n")

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
attr(res,"MI") <- get_results(res)
saveRDS(res, file = paste0(path,"/results/ScreenLearnSta_all.RDS"))
rm(res)

res <- parLapply(cl, infl$low, est_ScreenLearnSta)
attr(res,"MI") <- get_results(res)
saveRDS(res, file = paste0(path,"/results/ScreenLearnSta_low.RDS"))
rm(res)

res <- parLapply(cl, infl$high, est_ScreenLearnSta)
attr(res,"MI") <- get_results(res)
saveRDS(res, file = paste0(path,"/results/ScreenLearnSta_high.RDS"))
rm(res)

est_ScreenLearnDyn <- function(dat) {

  source(paste0(path,"LearnerLibrary.R"))
  
  # load nodes, g-/Q-formulas and interventions
  load(paste0(path,"NodesFormInterv.RData"))
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
attr(res,"MI") <- get_results(res)
saveRDS(res, file = paste0(path,"/results/ScreenLearnDyn_all.RDS"))
rm(res)

res <- parLapply(cl, infl$low, est_ScreenLearnDyn)
attr(res,"MI") <- get_results(res)
saveRDS(res, file = paste0(path,"/results/ScreenLearnDyn_low.RDS"))
rm(res)

res <- parLapply(cl, infl$high, est_ScreenLearnDyn)
attr(res,"MI") <- get_results(res)
saveRDS(res, file = paste0(path,"/results/ScreenLearnDyn_high.RDS"))
rm(res)

est_EconDAGSta <- function(dat) {
  
  source(paste0(path,"LearnerLibrary.R"))
  
  # load nodes, g-/Q-formulas and interventions
  load(paste0(path,"NodesFormInterv.RData"))
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
attr(res,"MI") <- get_results(res)
saveRDS(res, file = paste0(path,"/results/EconDAGSta_all.RDS"))
rm(res)

res <- parLapply(cl, infl$low, est_EconDAGSta)
attr(res,"MI") <- get_results(res)
saveRDS(res, file = paste0(path,"/results/EconDAGSta_low.RDS"))
rm(res)

res <- parLapply(cl, infl$high, est_EconDAGSta)
attr(res,"MI") <- get_results(res)
saveRDS(res, file = paste0(path,"/results/EconDAGSta_high.RDS"))
rm(res)

est_EconDAGDyn <- function(dat) {
  
  source(paste0(path,"LearnerLibrary.R"))
  
  # load nodes, g-/Q-formulas and interventions
  load(paste0(path,"NodesFormInterv.RData"))
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
attr(res,"MI") <- get_results(res)
saveRDS(res, file = paste0(path,"/results/EconDAGDyn_all.RDS"))
rm(res)

res <- parLapply(cl, infl$low, est_EconDAGDyn)
attr(res,"MI") <- get_results(res)
saveRDS(res, file = paste0(path,"/results/EconDAGDyn_low.RDS"))
rm(res)

res <- parLapply(cl, infl$high, est_EconDAGDyn)
attr(res,"MI") <- get_results(res)
saveRDS(res, file = paste0(path,"/results/EconDAGDyn_high.RDS"))
rm(res)

est_PlainDAGSta <- function(dat) {
  
  source(paste0(path,"LearnerLibrary.R"))
  
  # load nodes, g-/Q-formulas and interventions
  load(paste0(path,"NodesFormInterv.RData"))
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
attr(res,"MI") <- get_results(res)
saveRDS(res, file = paste0(path,"/results/PlainDAGSta_all.RDS"))
rm(res)

res <- parLapply(cl, infl$low, est_PlainDAGSta)
attr(res,"MI") <- get_results(res)
saveRDS(res, file = paste0(path,"/results/PlainDAGSta_low.RDS"))
rm(res)

res <- parLapply(cl, infl$high, est_PlainDAGSta)
attr(res,"MI") <- get_results(res)
saveRDS(res, file = paste0(path,"/results/PlainDAGSta_high.RDS"))
rm(res)

est_PlainDAGDyn <- function(dat) {
  
  source(paste0(path,"LearnerLibrary.R"))
  
  # load nodes, g-/Q-formulas and interventions
  load(paste0(path,"NodesFormInterv.RData"))
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
attr(res,"MI") <- get_results(res)
saveRDS(res, file = paste0(path,"/results/PlainDAGDyn_all.RDS"))
rm(res)

res <- parLapply(cl, infl$low, est_PlainDAGDyn)
attr(res,"MI") <- get_results(res)
saveRDS(res, file = paste0(path,"/results/PlainDAGDyn_low.RDS"))
rm(res)

res <- parLapply(cl, infl$high, est_PlainDAGDyn)
attr(res,"MI") <- get_results(res)
saveRDS(res, file = paste0(path,"/results/PlainDAGDyn_high.RDS"))
rm(res)

stopCluster(cl)
