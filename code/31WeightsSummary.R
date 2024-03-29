prep_res <- function(res) {

  # retrieve results from ltmle outputs and return matrix
  
  est <- res$estimate
  std <- res$std.dev
  CI <- res$CI
  pvalues <- res$pvalue
  cbind(est, std, CI, pvalues)
}

get_results <- function(res){
  
  # average results with Rubins Rule
  
  res_coll <- lapply(res, function(analyzed_set) analyzed_set$est_out$ltmle)
  res_coll <- lapply(res_coll, prep_res)
  res_coll <- do.call("rbind", res_coll)
  res_coll <- BaBooN::MI.inference(thetahat = res_coll[,"est"],
                                   varhat.thetahat = res_coll[,"std"]^2)
  cat(sprintf("ATE = %.2f | CI = [%.2f,%.2f]", res_coll$MI.Est, res_coll$CI.low, res_coll$CI.up),"\n")
  unlist(res_coll)
}

learner_weights_summary_Q <- function(ltmle_est, mean_tf = TRUE){
  
  # takes ltmle output and extracts the coefficient/weight of each learner in every 
  # Q-Model. It returns a vector, where each entry corresponds to the mean of this learner's weights across
  # all treatment and all control models (here: 2*11 points in time). Controls and Treatment weights for the same point
  # in time might be equal and merged anyway but this has no further consequence since the relative frequency is of interest.
  
  ltmle_est$fit$Q <- c(ltmle_est$fit$Q[[1]],ltmle_est$fit$Q[[2]]) # treatment and control
  bad_preds <- sapply(ltmle_est$fit$Q, function(x) is.character(x[[1]]))
  if(any(bad_preds)) {
    ltmle_est$fit$Q[bad_preds] <- NULL # no estimation occured
    message("Some Q-Learner had constant Ys and hence no estimation occured.")
  } 
  learner_weights <- sapply(ltmle_est$fit$Q, function(x) {
    try(return(x[,"Coef"]), silent = TRUE)
    x$coef
  })
  if(mean_tf == FALSE) return(learner_weights)
  mean_weights <- rowMeans(learner_weights)
  return(mean_weights)
}

learner_weights_summary_g <- function(ltmle_est, mean_tf = TRUE){
  
  # takes ltmle output and extracts the coefficient/weight of each learner in every 
  # g-Model. It returns a vector, where each entry corresponds to the mean of this learner's weights across
  # all treatment and all control models (here: 2*11 points in time). Controls and Treatment weights for the same point
  # in time might be equal and merged anyway but this has no further consequence since the relative frequency is of interest.
  
  ltmle_est$fit$g <- c(ltmle_est$fit$g[[1]],ltmle_est$fit$g[[2]]) # treatment and control
  bad_preds <- sapply(ltmle_est$fit$g, function(x) is.character(x[[1]]))
  if(any(bad_preds)) {
    ltmle_est$fit$g[bad_preds] <- NULL # no estimation occured
    message("Some g-Learner had constant Ys and hence no estimation occured.")
  } 
  learner_weights <- sapply(ltmle_est$fit$g, function(x) {
    try(return(x[,"Coef"]), silent = TRUE)
    x$coef
  })
  if(mean_tf == FALSE) return(learner_weights)
  mean_weights <- rowMeans(learner_weights)
  return(mean_weights)
}

# Different estimation strategies/specifications which are iterated in EstimateData.R

no_fm <- NULL
est_spec <- list("ScreenLearnSta_all" = list("d_set" = "all", "treat" = "stat_intv_1", "cntrl" = "stat_intv_0", "Sl_lib" = "SL.Est_Data", "g_form" = no_fm, "Q_form" = no_fm),
                 "ScreenLearnSta_high" = list("d_set" = "high", "treat" = "stat_intv_1", "cntrl" = "stat_intv_0", "Sl_lib" = "SL.Est_Data", "g_form" = no_fm, "Q_form" = no_fm),
                 "ScreenLearnSta_low" = list("d_set" = "low", "treat" = "stat_intv_1", "cntrl" = "stat_intv_0", "Sl_lib" = "SL.Est_Data", "g_form" = no_fm, "Q_form" = no_fm),
                 "ScreenLearnDyn_all" = list("d_set" = "all", "treat" = "dyn_intv", "cntrl" = "stat_intv_0", "Sl_lib" = "SL.Est_Data", "g_form" = no_fm, "Q_form" = no_fm),
                 "ScreenLearnDyn_high" = list("d_set" = "high", "treat" = "dyn_intv", "cntrl" = "stat_intv_0", "Sl_lib" = "SL.Est_Data", "g_form" = no_fm, "Q_form" = no_fm),
                 "ScreenLearnDyn_low" = list("d_set" = "low", "treat" = "dyn_intv", "cntrl" = "stat_intv_0", "Sl_lib" = "SL.Est_Data", "g_form" = no_fm, "Q_form" = no_fm),
                 "EconDAGSta_all" = list("d_set" = "all", "treat" = "stat_intv_1", "cntrl" = "stat_intv_0", "Sl_lib" = "SL.Est_Theory", "g_form" = "g_econ", "Q_form" = "Q_econ"),
                 "EconDAGSta_high" = list("d_set" = "high", "treat" = "stat_intv_1", "cntrl" = "stat_intv_0", "Sl_lib" = "SL.Est_Theory", "g_form" = "g_econ", "Q_form" = "Q_econ"),
                 "EconDAGSta_low" = list("d_set" = "low", "treat" = "stat_intv_1", "cntrl" = "stat_intv_0", "Sl_lib" = "SL.Est_Theory", "g_form" = "g_econ", "Q_form" = "Q_econ"),
                 "EconDAGDyn_all" = list("d_set" = "all", "treat" = "dyn_intv", "cntrl" = "stat_intv_0", "Sl_lib" = "SL.Est_Theory", "g_form" = "g_econ", "Q_form" = "Q_econ"),
                 "EconDAGDyn_high" = list("d_set" = "high", "treat" = "dyn_intv", "cntrl" = "stat_intv_0", "Sl_lib" = "SL.Est_Theory", "g_form" = "g_econ", "Q_form" = "Q_econ"),
                 "EconDAGDyn_low" = list("d_set" = "low", "treat" = "dyn_intv", "cntrl" = "stat_intv_0", "Sl_lib" = "SL.Est_Theory", "g_form" = "g_econ", "Q_form" = "Q_econ"),
                 "PlainDAGSta_all" = list("d_set" = "all", "treat" = "stat_intv_1", "cntrl" = "stat_intv_0", "Sl_lib" = "SL.Est_Theory", "g_form" = "g_base", "Q_form" = "Q_base"),
                 "PlainDAGSta_high" = list("d_set" = "high", "treat" = "stat_intv_1", "cntrl" = "stat_intv_0", "Sl_lib" = "SL.Est_Theory", "g_form" = "g_base", "Q_form" = "Q_base"),
                 "PlainDAGSta_low" = list("d_set" = "low", "treat" = "stat_intv_1", "cntrl" = "stat_intv_0", "Sl_lib" = "SL.Est_Theory", "g_form" = "g_base", "Q_form" = "Q_base"),
                 "PlainDAGDyn_all" = list("d_set" = "all", "treat" = "dyn_intv", "cntrl" = "stat_intv_0", "Sl_lib" = "SL.Est_Theory", "g_form" = "g_base", "Q_form" = "Q_base"),
                 "PlainDAGDyn_high" = list("d_set" = "high", "treat" = "dyn_intv", "cntrl" = "stat_intv_0", "Sl_lib" = "SL.Est_Theory", "g_form" = "g_base", "Q_form" = "Q_base"),
                 "PlainDAGDyn_low" = list("d_set" = "low", "treat" = "dyn_intv", "cntrl" = "stat_intv_0", "Sl_lib" = "SL.Est_Theory", "g_form" = "g_base", "Q_form" = "Q_base"))
