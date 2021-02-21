# Please install the following packages prior to running the file: 
# ltmle, SuperLearner, vcd, arm, rpart, nnet, glmnet, stringr, magrittr, randomForest,
# earth, gbm, gam, mgcv, reshape2, dplyr, data.table, polspline
# plyr, tibble, scales, BaBooN, simcausal (currently only on CRAN archive)

# load in memory and attach to serachpath
rm(list = ls())
library(parallel)
set.seed(1)
path <- "/cluster/home/phibauma/CausalInflation/"
#path <- "/Users/flipst3r/RStHomeDir/GitHub/CausalInflation/"

# specify here how many cores are available for parallel computation (should be 5 here)
n_cluster <- 5

# load 5 imputed data.frames that are analyzed
load(paste0(path,"causalinfl_revised.RData"))
all <- infl$all # needed so that get() works with characters/strings in the estimation environment
high <- infl$high
low <- infl$low

# extract learner weights
source(paste0(path,"WeightsSummary.R"))

# formula to extract ATE
get_ATE <- function(out, est = "tmle") unclass(summary(out, estimator = est))$effect.measures$ATE

# Learner Sets: SL.Est_Theory, SL.Est_Data
load(paste0(path,"SelectedLearner.RData"))

# Nodes, g-/Q-formulas and interventions
load(paste0(path,"NodesFormInterv.RData")) 
int <- ltmle_prep$invterventions
stat_intv_0 <- int$stat0
stat_intv_1 <- int$stat1
dyn_intv <- int$dyn_intv
nod <- ltmle_prep$nodes
g_econ <- ltmle_prep$form_econ$g_econ
Q_econ <- ltmle_prep$form_econ$Q_econ
g_base <- ltmle_prep$form_base$g_base
Q_base <- ltmle_prep$form_base$Q_base

# initiate cluster
cl <- makeCluster(n_cluster, outfile = "")
clusterSetRNGStream(cl = cl, iseed = 1)
clusterEvalQ(cl, library(ltmle))

source(paste0(path,"LearnerLibrary.R"))
#SL.Est_Theory <- SL.Est_Data <- c("SL.glm","SL.mean") # for fast test runs
SL.Est_Theory <- SL.Est_Theory[-10] # gbm consumes too much memory and does not help much either

clusterExport(cl = cl, list(
  "learner_weights_summary_g", "learner_weights_summary_Q","get_ATE","cc_trunc","nod","path"
))

est <- function(d_s, Q_fm, g_fm, treat, cntrl, SL_lib) {
  
  # workhorse for all ltmle estimations
  
  seed <- .Random.seed
  source(paste0(path,"LearnerLibrary.R")) # for all manual learners
  
  res <- tryCatch(
    {
      ltmle(d_s,
            Qform = Q_fm, gform = g_fm,
            Lnodes = nod$L, Anodes = nod$A, Ynodes = nod$Y,
            estimate.time = FALSE,
            Yrange = c(-3, 86),
            variance.method = "tmle",
            abar = list(treament = treat, control = cntrl),
            SL.library = SL_lib
      )
    },
    error = function(cond) {
      return(NA)
    }
  )
  si_ze <- format(object.size(res), units = "auto")
  
  # extract learner weights from estimation
  Q_mean <- try(learner_weights_summary_Q(res), silent = TRUE)
  g_mean <- try(learner_weights_summary_g(res), silent = TRUE)
  learn_weights <- list(Qweights = Q_mean, gweights = g_mean)
  
  # extract ATEs from estimation
  ests <- list(
    ltmle = try(get_ATE(res), silent = TRUE),
    iptw = try(get_ATE(res, est = "iptw"), silent = TRUE)
  )
  
  res <- list(est_out = ests, weights_out = learn_weights, cc = try(cc_trunc(res), silent = TRUE))
  attr(res,"info") <- sessionInfo()
  attr(res,"seed") <- seed
  attr(res,"time") <- Sys.time()
  attr(res,"size") <- si_ze
  res
}

res <- lapply(est_spec, function(sp) {
  
  # iterate over every estimation strategy to reduce the copy/paste amount regarding
  # est(). est_spec is a list where character entries can be used in the estimation
  # environment to load the corresponding estimation specfication (e.g. Q-formula)

  SL_lib <- get(sp$Sl_lib)
  treat <- get(sp$treat)
  cntrl <- get(sp$cntrl)
  g_fm <- if (is.null(sp$g_form)) NULL else get(sp$g_form)
  Q_fm <- if (is.null(sp$Q_form)) NULL else get(sp$Q_form)
  d <- get(sp$d_set)
  
  # handle dynamic treatment in subgroups (i.e. high-income) by rownames of matrix
  # which are the row indices (1-124) of the imputed data.frame
  ids <- rownames(d[[1]])
  is_static <- is.null(rownames(treat))
  if (!is_static) treat <- treat[rownames(treat) %in% ids,] 
  if (is_static & length(ids) < 59) treat <- treat[1:length(ids),]
  cntrl <- cntrl[1:nrow(d[[1]]),]
  r <- parLapply(cl, d, est, Q_fm = Q_fm, g_fm = g_fm, treat = treat, cntrl = cntrl, SL_lib = SL_lib)
  attr(r,"MI") <- get_results(r)
  r
  
})
stopCluster(cl)
saveRDS(res, file = paste0(path,"/results/Estimations.RDS"))
res <- sapply(res, function(x) attributes(x)$MI)
saveRDS(res, file = paste0(path,"/results/ATEs.RDS"))
