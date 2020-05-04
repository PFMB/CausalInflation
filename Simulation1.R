# load in memory and attach to serachpath
rm(list = ls())
library(simcausal)  # needs to be installed from CRAN (Archive):
                    # https://cran.r-project.org/web/packages/simcausal/index.html
library(parallel)
library(ltmle)
set.seed(1)

# insert working directory here
setwd("")

# setup
runs      <- 2000  # number of simulation runs
n_cluster <-  48   # specify here how many cores are available for parallel computation
n_obs <- 1e3       # number of iid observations in simulated (observed) data sets

# ------- DEFINE DGP ------- #

## Simulation 1
t.start <- 2000L
t.end <- 2002L

Sim1_DAG <- DAG.empty() +
  node("L", t = t.start, distr = "rnorm", mean = 0, sd = 0.5) +
  node("L", t = (t.start + 1):t.end, distr = "rnorm", mean = L[t - 1] + A[t - 1], sd = 0.5) +
  node("A", t = t.start, distr = "rbern", prob = plogis(L[t])) +
  node("A", t = (t.start + 1):t.end, distr = "rbern", prob = plogis(L[t] + A[t - 1] * 2 - L[t - 1])) +
  node("Y", t = t.start, distr = "rnorm", mean = A[t] * 50 + L[t], sd = 0.5) +
  node("Y", t = (t.start + 1):t.end, distr = "rnorm", mean = A[t] * 50 + L[t - 1] + L[t] + Y[t - 1], sd = 0.25)

D <- set.DAG(Sim1_DAG)

## Set Interventions
# Add A1 to DAG: A1=(1 1 1)
action_A1 <- c(node("A",
  t = t.start:t.end, distr = "rbern",
  prob = 1
))
D <- D + action("A1", nodes = action_A1)

# Add A0 to DAG: A0=(0 0 0)
action_A0 <- c(node("A",
  t = t.start:t.end, distr = "rbern",
  prob = 0
))
D <- D + action("A0", nodes = action_A0)

# ------- SIMULATE OBSERVED DATA ------- #

# preallocate
# Simulation 1: vector with number of intended runs
Obs_dat <- vector("list", runs)

# simulate observed data
for (ind in seq(Obs_dat)) {
  Obs_dat[[ind]] <- sim(D, n = n_obs, verbose = FALSE) # observed data
}

# ltmle does not need ID later on
Obs_dat <- lapply(Obs_dat, function(x) x[, -1, drop = FALSE])

# ------- SIMULATE COUNTERFACTUAL DATA ------- #

# counterfactual data set with 1 million draws for the true ATE
counter_dat <- sim(D, n = 1e6, actions = c("A0", "A1"), verbose = FALSE, rndseed = 1)

# Define parameter of interest: ATE
# For 3/6 time points - A1: 1 1 1 vs A0: 0 0 0
D <- set.targetE(D, outcome = "Y", t = t.end, param = "A1-A0")
true_ATE <- eval.target(D, data = counter_dat)$res

# ------- RUN ANALYSIS ------- #

## Q- and g- Model are correctly specified

# Initiate cluster
cl <- makeCluster(n_cluster, outfile = "")
clusterSetRNGStream(cl = cl, iseed = 1)
clusterEvalQ(cl, library(ltmle))

exe_Sim1 <- function(x) {
  SL.Set <- c("SL.glm")
  attr(SL.Set, "return.fit") <- TRUE

  # define static interventions
  treatment_1 <- matrix(1, nrow = nrow(x), ncol = length(grep("A", names(x))))
  control_0 <- matrix(0, nrow = nrow(x), ncol = length(grep("A", names(x))))

  output <- try(ltmle(x,
    Anodes = grep("A", names(x)),
    Ynodes = grep("Y", names(x)),
    Lnodes = grep("L", names(x)),
    Cnodes = NULL,
    Yrange = c(-15, 169),
    gbounds = c(0.01, 1),
    abar = list(treament = treatment_1, control = control_0),
    SL.library = SL.Set
  ), silent = FALSE)

  ltmle_ATE <- unclass(summary(output)$effect.measures$ATE)
  iptw_ATE <- unclass(summary(output, estimator = "iptw")$effect.measures$ATE)
  est_output_temp <- list(ltmle = ltmle_ATE, iptw = iptw_ATE)

  return(est_output_temp)
}

# run estimation
t_ime <- system.time({
  Sim1 <- parLapply(cl, Obs_dat, exe_Sim1)
})

# save results
(attributes(Sim1)$time <- t_ime)
(attributes(Sim1)$sessinfo <- sessionInfo())
(attributes(Sim1)$seed <- .Random.seed)
saveRDS(Sim1, file = "Sim1.RDS")

## g-Model is correctly and Q-Model is misspecified

mis_Q_form <- vector(length = 3)
names(mis_Q_form) <- c("Y_2000", "Y_2001", "Y_2002")
mis_Q_form[1] <- c("Q.kplus1 ~ A_2000")
mis_Q_form[2] <- c("Q.kplus1 ~ A_2000 + Y_2000 + A_2001")
mis_Q_form[3] <- c("Q.kplus1 ~ A_2000 + Y_2000 + A_2001 + Y_2001 + A_2002")

clusterExport(cl, "mis_Q_form")

exe_Sim1_mis <- function(x) {
  SL.Set <- c("SL.glm")

  # define static interventions
  treatment_1 <- matrix(1, nrow = nrow(x), ncol = length(grep("A", names(x))))
  control_0 <- matrix(0, nrow = nrow(x), ncol = length(grep("A", names(x))))

  output <- try(ltmle(x,
    Anodes = grep("A", names(x)),
    Ynodes = grep("Y", names(x)),
    Lnodes = grep("L", names(x)),
    Cnodes = NULL,
    Yrange = c(-15, 169),
    gbounds = c(0.01, 1),
    Qform = mis_Q_form,
    abar = list(treament = treatment_1, control = control_0),
    SL.library = SL.Set
  ), silent = FALSE)

  ltmle_ATE <- unclass(summary(output)$effect.measures$ATE)
  iptw_ATE <- unclass(summary(output, estimator = "iptw")$effect.measures$ATE)
  list(ltmle = ltmle_ATE, iptw = iptw_ATE)
}

# run estimation
t_ime <- system.time({
  Sim1_mis <- parLapply(cl, Obs_dat, exe_Sim1_mis)
})

# save results
(attributes(Sim1_mis)$time <- t_ime)
(attributes(Sim1_mis)$sessinfo <- sessionInfo())
(attributes(Sim1_mis)$seed <- .Random.seed)
saveRDS(Sim1_mis, file = "Sim1Mis.RDS")

stopCluster(cl)

# ------- PRINT RESULTS ------- #

source("CalcBiasCP.R")

cat("GLM: correct Q-formula: \n")
(get_bias_cp(do.call("c", Sim1), true_ATE))
cat("GLM: incorrect Q-formula: \n")
(get_bias_cp(do.call("c", Sim1_mis), true_ATE))
