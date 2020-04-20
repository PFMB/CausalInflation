########################################################################################################
# Philipp Baumann, Michael Schomaker, and Enzo Rossi
# Working Paper (under submission)
# Title:
# Estimating the Effect of Central Bank Independence on Inflation Using 
# Longitudinal Targeted Maximum Likelihood Estimation
########################################################################################################

# Copyright (c) 2020  <Philipp Baumann, Michael Schomaker, and Enzo Rossi>
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

##############################################################################################
# required packages
rm(list = ls())
library(simcausal)  # needs to be installed from CRAN (Archive):
                    # https://cran.r-project.org/web/packages/simcausal/index.html
library(parallel)
set.seed(123)

# insert your working directory here
setwd("C:/temp/")

# setup
runs      <- 2000                      # number of simulation runs
n.cluster <- parallel::detectCores()   # specify here how many cores are available for parallel computation

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
# Simulation 1: 2000 df with n = 1000
Obs_dat <- vector("list", runs)

# simulate observed data
for (ind in seq(Obs_dat)) {
  Obs_dat[[ind]] <- sim(D, n = 1000, verbose = FALSE) # observed data
}

# ltmle does not need ID later on
Obs_dat <- lapply(Obs_dat, function(x) x[, -1, drop = FALSE])

# ------- SIMULATE COUNTERFACTUAL DATA ------- #

# counterfactual data set with 1 million draws for the true ATE
counter_dat <- sim(D, n = 1e6, actions = c("A0", "A1"), verbose = FALSE, rndseed = 123)

# Define parameter of interest: ATE
# For 3/6 time points - A1: 1 1 1 vs A0: 0 0 0
D <- set.targetE(D, outcome = "Y", t = t.end, param = "A1-A0")
true_ATE <- eval.target(D, data = counter_dat)$res

# ------- RUN ANALYSIS ------- #

## Q- and g- Model are correctly specified

# Initiate cluster
cl <- makeCluster(n.cluster, outfile = "")
clusterSetRNGStream(cl = cl, 123)
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

start_time <- Sys.time()
Sim1 <- parLapply(cl, Obs_dat, exe_Sim1)
end_time <- Sys.time()

cat("Elapsed time:", end_time - start_time,"sec. \n")

## g-Model is correctly and Q-Model is misspecified

mis_Q_form <- vector(length = 3)
names(mis_Q_form) <- c("Y_2000", "Y_2001", "Y_2002")
mis_Q_form[1] <- c("Q.kplus1 ~ A_2000")
mis_Q_form[2] <- c("Q.kplus1 ~ A_2000 + Y_2000 + A_2001")
mis_Q_form[3] <- c("Q.kplus1 ~ A_2000 + Y_2000 + A_2001 + Y_2001 + A_2002")

clusterExport(cl, "mis_Q_form")

exe_Sim1_mis <- function(x) {
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
    Qform = mis_Q_form,
    abar = list(treament = treatment_1, control = control_0),
    SL.library = SL.Set
  ), silent = FALSE)

  ltmle_ATE <- unclass(summary(output)$effect.measures$ATE)
  iptw_ATE <- unclass(summary(output, estimator = "iptw")$effect.measures$ATE)
  list(ltmle = ltmle_ATE, iptw = iptw_ATE)
}

start_time <- Sys.time()
Sim1_mis <- parLapply(cl, Obs_dat, exe_Sim1_mis)
end_time <- Sys.time()

stopCluster(cl)

sessionInfo()
.Random.seed

cat("Elapsed time:", end_time - start_time,"sec. \n")

# ------- GET RESULTS ------- #

source("CalcBiasCP.R")
(get_bias_cp(do.call("c", Sim1), true_ATE))
(get_bias_cp(do.call("c", Sim1_mis), true_ATE))
