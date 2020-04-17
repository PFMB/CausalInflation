get_bias_cp <- function(est_output, true_ATE) {

  #        x: list with ltmle output summaries from the LTMLE est and estimator="iptw"
  #           estimate. ltmle estimations were run on simulated data sets.
  #           The list should contain at least 1000 successful model fits.
  #
  # true_ATE: Result of SimCausal on a counterfactual dataset with 10000000 draws (imitating the population)
  #
  # ouput:    a matrix indicating the bias, coverage probability and the share of
  #           NA's (failed estimations) and further more for the ltml and the iptw estimator

  # seperate ltmle results from iptw results since both are alternately contained
  raw_ltmle_results <- est_output[which(names(est_output) == "ltmle")]
  raw_iptw_results <- est_output[which(names(est_output) == "iptw")]

  get_estimates <- function(x) {

    # function to extract ATE and CI estimates from the summary fits

    ate <- x$estimate
    ci <- as.vector(x$CI)
    res <- c(ate, ci)
    names(res) <- c("ATE", "2.5%", "97.5%")
    return(res)
  }

  # LTMLE
  all_ltmle_results <- lapply(raw_ltmle_results, get_estimates)
  all_ltmle_ate_results <- unlist(lapply(all_ltmle_results, `[`, "ATE"))

  # IPTW
  all_iptw_results <- lapply(raw_iptw_results, get_estimates)
  all_iptw_ate_results <- unlist(lapply(all_iptw_results, `[`, "ATE"))

  # Are there failed entries (NA's) among all estimations?
  failed_share_ltmle <- (sum(is.na(all_ltmle_ate_results)) / length(all_ltmle_ate_results)) * 100
  failed_share_iptw <- (sum(is.na(all_iptw_ate_results)) / length(all_iptw_ate_results)) * 100

  ## find and extract the first 1000 successful estimations where lmtle AND iptw are not NA's

  good_est_idx <- which((!is.na(all_ltmle_ate_results)) & (!is.na(all_iptw_ate_results)))
  if (length(good_est_idx) < 1000) warning("Less than 1000 successful estimations.")

  # take the first 1000 or take whats available (might be far less than 1000)
  if (length(good_est_idx) > 1000) good_est_idx <- good_est_idx[seq(1000)]

  # LTMLE
  ltmle_results <- lapply(raw_ltmle_results[good_est_idx], get_estimates)
  ltmle_ate_results <- unlist(lapply(ltmle_results, `[`, "ATE"))
  ltmle_lower_ci <- unlist(lapply(ltmle_results, `[`, "2.5%"))
  ltmle_upper_ci <- unlist(lapply(ltmle_results, `[`, "97.5%"))

  # IPTW
  iptw_results <- lapply(raw_iptw_results[good_est_idx], get_estimates)
  iptw_ate_results <- unlist(lapply(iptw_results, `[`, "ATE"))
  iptw_lower_ci <- unlist(lapply(iptw_results, `[`, "2.5%"))
  iptw_upper_ci <- unlist(lapply(iptw_results, `[`, "97.5%"))

  ## LTMLE ##

  # Abs. bias of all ltmle's ATE
  mean_ate_ltmle <- mean(ltmle_ate_results)
  bias_ate_ltmle <- abs(mean(ltmle_ate_results - true_ATE))
  names(bias_ate_ltmle) <- "Bias LTMLE"

  # Coverage Probability of Ltmle's CI
  bin_ind_ltmle <- rep(0, length(ltmle_lower_ci))
  bin_ind_ltmle[(ltmle_lower_ci <= true_ATE) & (true_ATE <= ltmle_upper_ci)] <- 1
  cov_prob_ltmle <- mean(bin_ind_ltmle) * 100
  names(cov_prob_ltmle) <- "CP (%) LTMLE"

  ## IPTW ##

  # Abs. bias of all IPTW's ATE
  mean_ate_iptw <- mean(iptw_ate_results)
  bias_ate_iptw <- abs(mean(iptw_ate_results - true_ATE))
  names(bias_ate_iptw) <- "Bias IPTW"

  # Coverage Probability of IPTW's CI
  bin_ind_iptw <- rep(0, length(iptw_lower_ci))
  bin_ind_iptw[(iptw_lower_ci <= true_ATE) & (true_ATE <= iptw_upper_ci)] <- 1
  cov_prob_iptw <- mean(bin_ind_iptw) * 100
  names(cov_prob_iptw) <- "CP (%) IPTW"
  
  c(bias_ate_ltmle,cov_prob_ltmle,bias_ate_iptw,cov_prob_iptw)
}
