get_bias_cp <- function(est_output, true_ATE, take_first = 10L) {

  # est_output: list with estimation output summaries of LTMLE and IPTW (estimator="iptw").
  #             est_output should contain at least 1000 successful model fits but can
  #             also handle failed estimations (i.e. try-error classes).
  #
  # true_ATE: scalar containing the result of simcausal on a counterfactual dataset
  #           with 1e6 draws (imitating the population)
  #
  # take_first: scalar indicating how many successful estimations should be taken
  #             to calculate the abs. bias and coverage probability
  #
  # ouput: a named list indicating the bias and coverage probability
  #        for the LTMLE and the IPTW estimator.

  estimates <- sapply(est_output, get_estimates)
  
  # we want only estimates where both estimations (LTMLE and IPTW) were successful
  grps <- cut(seq(ncol(estimates)), seq(0, ncol(estimates), 2), include.lowest = TRUE)
  grps <- do.call("cbind", split(seq(ncol(estimates)), grps))
  succ_grps <- apply(grps, 2, function(est_pair) {
    # FALSE if either LTMLE or IPTW failes to guarantee a fair comparison
    if (any(is.na(estimates[, est_pair]))) {
      return(FALSE)
    }
    return(TRUE)
  })
  cat("Share of either failed LMTLE or IPTW estimations:", 1 - mean(succ_grps),"\n")
  succ_grps <- grps[, succ_grps]
  succ_grps <- succ_grps[, 1:take_first]
  succ_estim <- estimates[, as.vector(succ_grps)]

  list(
    "ltmle" = calc_res(succ_estim[, colnames(succ_estim) == "ltmle"]),
    "iptw" = calc_res(succ_estim[, colnames(succ_estim) == "iptw"])
  )
}

get_estimates <- function(est) {

  # function to extract ATE and CI estimates from the summary fits

  if (class(est) == "try-error") {
    return(setNames(rep(NA, 3), c("ATE", "2.5%", "97.5%")))
  }
  ate <- est$estimate
  ci <- as.vector(est$CI)
  setNames(c(ate, ci), c("ATE", "2.5%", "97.5%"))
}

calc_res <- function(est_vec) {
  
  # Absolute Bias and coverage probability
  # est_vec matrix of successfull fits of either LTMLE or IPTW
  
  ate <- est_vec["ATE", ]
  lower <- est_vec["2.5%", ]
  upper <- est_vec["97.5%", ]
  c("abs_bias" = abs(mean(ate - true_ATE)), "cp" = mean(true_ATE >= lower & true_ATE <= upper))
}