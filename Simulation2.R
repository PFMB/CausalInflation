# Please install the following packages prior to running the file: 
# ltmle, SuperLearner, vcd, arm, rpart, nnet, glmnet, randomForest,
# earth, gbm, gam, mgcv, reshape2, dplyr, data.table, polspline
# plyr, tibble, scales, BaBooN, simcausal (currently only on CRAN archive)

# load in memory and attach to serachpath
rm(list = ls())
library(simcausal)  # needs to be installed from CRAN (Archive):
                    # https://cran.r-project.org/web/packages/simcausal/index.html
library(parallel)
set.seed(1)

# insert working directory here
setwd("/cluster/home/scstepha/CausalInflation")

# setup
runs      <- 1500  # number of simulation runs (take into account failing estimations)
n_cluster <- 24    # specify here how many cores are available for parallel computation
n_obs <- 1e3       # number of iid observations in simulated (observed) data sets

# ------- DEFINE DGP ------- #

## Simulation 2
t.start <- 2000L
t.end <- 2005L

Sim2_DAG <- DAG.empty() +
  node("L.1", t = t.start, distr = "rnorm", mean = 0, sd = 0.5) +
  node("L.1", t = (t.start + 1):t.end, distr = "rnorm", mean = L.7[t - 1], sd = 0.5) +
  node("A", t = t.start, distr = "rbern", prob = plogis(L.1[t])) +
  node("A", t = (t.start + 1):t.end, distr = "rbern", prob = plogis(0.25 * L.1[t] + 0.25 * L.6[t - 1])) +
  node("Y", t = t.start, distr = "rnorm", mean = A[t] + L.1[t], sd = 3) +
  node("Y", t = (t.start + 1):t.end, distr = "rnorm", mean = A[t] + L.1[t] + L.9[t - 1] + L.10[t - 1] * 0.05, sd = 0.5) +
  node("L.2", t = t.start:(t.end - 1), distr = "rnorm", mean = A[t] + L.1[t], sd = 0.5) +
  node("L.3", t = t.start:(t.end - 1), distr = "rnorm", mean = Y[t] + L.2[t], sd = 1) +
  node("L.4", t = t.start:(t.end - 1), distr = "rnorm", mean = A[t], sd = 0.5) +
  node("L.5", t = t.start, distr = "rnorm", mean = Y[t], sd = 1.5) +
  node("L.5", t = (t.start + 1):(t.end - 1), distr = "rnorm", mean = Y[t] + L.10[t - 1], sd = 1.5) +
  node("L.6", t = t.start:(t.end - 1), distr = "rnorm", mean = L.4[t], sd = 0.5) +
  node("L.7", t = t.start:(t.end - 1), distr = "rnorm", mean = L.2[t], sd = 0.5) +
  node("L.8", t = t.start:(t.end - 1), distr = "rnorm", mean = L.5[t], sd = 0.5) +
  node("L.9", t = t.start:(t.end - 1), distr = "rnorm", mean = L.3[t], sd = 1) +
  node("L.10", t = t.start:(t.end - 1), distr = "rnorm", mean = L.8[t] + L.9[t], sd = 0.5)

D <- set.DAG(Sim2_DAG)

# Set Interventions

# Add A1 to DAG: A1=(1 1 1 1 1 1)
action_A1 <- c(node("A",
  t = t.start:t.end, distr = "rbern",
  prob = 1
))
D <- D + action("A1", nodes = action_A1)

# Add A0 to DAG: A0=(0 0 0 0 0 0)
action_A0 <- c(node("A",
  t = t.start:t.end, distr = "rbern",
  prob = 0
))
D <- D + action("A0", nodes = action_A0)

# ------- SIMULATE OBSERVED DATA ------- #

# preallocate
# Simulation 2: vector with number of intended runs
Obs_dat <- vector("list", runs)

# simulate observed data
for (ind in seq(Obs_dat)) {
  Obs_dat[[ind]] <- sim(D, n = n_obs, verbose = FALSE) # observed data
}

# ltmle does not need ID later on
Obs_dat <- lapply(Obs_dat, function(x) x[, -1, drop = FALSE])

# ------- SIMULATE COUNTERFACTUAL DATA ------- #

# counterfactual data set with 1e6 draws for the true ATE
counter_dat <- sim(D, n = 1e6, actions = c("A0", "A1"), verbose = FALSE, rndseed = 123)

# Define parameter of interest: ATE
# For 3 time points - A1: 1 1 1 1 1 1 vs A0: 0 0 0 0 0 0
D <- set.targetE(D, outcome = "Y", t = t.end, param = "A1-A0")
true_ATE <- eval.target(D, data = counter_dat)$res

# ----- correct Q-FORMULA ----- #

new_Q_form <- vector(length = 6)
names(new_Q_form) <- c("Y_2000", "Y_2001", "Y_2002", "Y_2003", "Y_2004", "Y_2005")
new_Q_form[1] <- c("Q.kplus1 ~ A_2000 + L.1_2000")
new_Q_form[2] <- c("Q.kplus1 ~ A_2001 + L.1_2001 + L.9_2000 + L.10_2000")
new_Q_form[3] <- c("Q.kplus1 ~ A_2002 + L.1_2002 + L.9_2001 + L.10_2001")
new_Q_form[4] <- c("Q.kplus1 ~ A_2003 + L.1_2003 + L.9_2002 + L.10_2002")
new_Q_form[5] <- c("Q.kplus1 ~ A_2004 + L.1_2004 + L.9_2003 + L.10_2003")
new_Q_form[6] <- c("Q.kplus1 ~ A_2005 + L.1_2005 + L.9_2004 + L.10_2004")

# ----- correct g-FORMULA ----- #

new_g_form <- vector(length = 6)
names(new_g_form) <- c("A_2000", "A_2001", "A_2002", "A_2003", "A_2004", "A_2005")
new_g_form[1] <- c("A_2000 ~ L.1_2000")
new_g_form[2] <- c("A_2001 ~ L.1_2001 + L.6_2000")
new_g_form[3] <- c("A_2002 ~ L.1_2002 + L.6_2001")
new_g_form[4] <- c("A_2003 ~ L.1_2003 + L.6_2002")
new_g_form[5] <- c("A_2004 ~ L.1_2004 + L.6_2003")
new_g_form[6] <- c("A_2005 ~ L.1_2005 + L.6_2004")

correct_forms <- list(gforms = new_g_form, Qforms = new_Q_form)

# ----- incorrect Q-formula ----- #

incor_Q_form <- new_Q_form
incor_Q_form[1] <- c("Q.kplus1 ~ A_2000 ")
incor_Q_form[2] <- c("Q.kplus1 ~ A_2001 + L.9_2000 + L.10_2000")
incor_Q_form[3] <- c("Q.kplus1 ~ A_2002 + L.9_2001 + L.10_2001")
incor_Q_form[4] <- c("Q.kplus1 ~ A_2003 + L.9_2002 + L.10_2002")
incor_Q_form[5] <- c("Q.kplus1 ~ A_2004 + L.9_2003 + L.10_2003")
incor_Q_form[6] <- c("Q.kplus1 ~ A_2005 + L.9_2004 + L.10_2004")

incor_forms <- list(gforms = new_g_form, Qforms = incor_Q_form)

# formula to extract ATE
get_ATE <- function(out, est = "tmle") unclass(summary(out, estimator = est))$effect.measures$ATE

# Initiate cluster
cl <- makeCluster(n_cluster, outfile = "")
clusterSetRNGStream(cl = cl, iseed = 1)

# define static interventions
treatment_1 <- matrix(1, nrow = n_obs, ncol = 6)
control_0 <- matrix(0, nrow = n_obs, ncol = 6)

load("SelectedLearner.RData") # predefined Learner Sets

clusterEvalQ(cl, library(ltmle))
clusterExport(cl = cl, list(
  "correct_forms", "incor_forms", "treatment_1",
  "control_0", "get_ATE"
))

# ----- Run Simulation ----- #

exe <- function(y) {
  cat("Estimation with correct Q-formulas starts.\n")

  source("LearnerLibrary.R") # individual learner

  L_set <- try(y$learner, silent = TRUE)
  x <- try(y$data, silent = TRUE)

  cor <- try(ltmle(x,
    Anodes = grep("A", names(x)),
    Ynodes = grep("Y", names(x)),
    Lnodes = grep("L", names(x)),
    Qform = correct_forms$Qforms,
    gform = correct_forms$gforms,
    Yrange = c(-102, 148),
    gbounds = c(0.01, 1),
    abar = list(treament = treatment_1, control = control_0),
    SL.library = L_set
  ), silent = TRUE)

  # extract ATEs from estimation
  cor_ests <- list(
    ltmle = try(get_ATE(cor), silent = TRUE),
    iptw = try(get_ATE(cor, est = "iptw"), silent = TRUE)
  )

  cat("Estimation with INCORRECT Q-formulas starts.\n")

  incor <- try(ltmle(x,
    Anodes = grep("A", names(x)),
    Ynodes = grep("Y", names(x)),
    Lnodes = grep("L", names(x)),
    Qform = incor_forms$Qforms,
    gform = incor_forms$gforms,
    Yrange = c(-102, 148),
    abar = list(treament = treatment_1, control = control_0),
    SL.library = L_set
  ), silent = TRUE)

  # extract ATEs from estimation
  incor_ests <- list(
    ltmle = try(get_ATE(incor), silent = TRUE),
    iptw = try(get_ATE(incor, est = "iptw"), silent = TRUE)
  )

  list(correct = cor_ests, incorrect = incor_ests)
}

# assemble in list of lists since we can only handle one alternating argument in
# parlapply()
list1 <- c(Obs_dat, Obs_dat, Obs_dat)
reps <- rep(1, runs)
list2 <- c(list(SL.Set1)[reps], list(SL.Set2)[reps], list(SL.Set3)[reps])
data_n_learner <- lapply(seq(list1), function(idx) {
  list(data = list1[[idx]], learner = list2[[idx]])
})

# run estimation
t_ime <- system.time({
  Sim2 <- parLapply(cl, data_n_learner, exe)
})
stopCluster(cl)

# save results
(attributes(Sim2)$time <- t_ime)
(attributes(Sim2)$sessinfo <- sessionInfo())
(attributes(Sim2)$seed <- .Random.seed)
saveRDS(Sim2, file = "Sim2.RDS")

# ------- PRINT RESULTS ------- #

source("CalcBiasCP.R")

# seperate results by learner sets
Sim2_L1 <- Sim2[1:runs]
Sim2_L2 <- Sim2[(runs + 1):(2 * runs)]
Sim2_L3 <- Sim2[(2 * runs + 1):(3 * runs)]

cat("L1: correct and incorrect Q-formula: \n")
res_all_L1 <- do.call("c", Sim2_L1)
res_corr_L1 <- do.call("c", unname(res_all_L1[names(res_all_L1) == "correct"]))
res_incorr_L1 <- do.call("c", unname(res_all_L1[names(res_all_L1) == "incorrect"]))
(get_bias_cp(res_corr_L1, true_ATE))
(get_bias_cp(res_incorr_L1, true_ATE))

cat("L2: correct and incorrect Q-formula: \n")
res_all_L2 <- do.call("c", Sim2_L2)
res_corr_L2 <- do.call("c", unname(res_all_L2[names(res_all_L2) == "correct"]))
res_incorr_L2 <- do.call("c", unname(res_all_L2[names(res_all_L2) == "incorrect"]))
(get_bias_cp(res_corr_L2, true_ATE))
(get_bias_cp(res_incorr_L2, true_ATE))

cat("L3: correct and incorrect Q-formula: \n")
res_all_L3 <- do.call("c", Sim2_L3)
res_corr_L3 <- do.call("c", unname(res_all_L3[names(res_all_L3) == "correct"]))
res_incorr_L3 <- do.call("c", unname(res_all_L3[names(res_all_L3) == "incorrect"]))
(get_bias_cp(res_corr_L3, true_ATE))
(get_bias_cp(res_incorr_L3, true_ATE))

