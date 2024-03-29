######
# library of wrapper functions/learner for variable selection and estimation/prediction
# for use in SuperLearner/LTMLE
######

#### Screening ####

### 1) Random Forest with 'randomForest::randomForest'

screen.randomForest_base <- function(Y, X, family = list(), nVar = 8, ntree = 200, mtry = ifelse(family$family ==
                                       "gaussian", floor(sqrt(ncol(X))), max(floor(ncol(X) / 3), 1)),
                                     nodesize = ifelse(family$family == "gaussian", 5, 1), maxnodes = NULL,
                                     ...) {
  # chose family dependent upon response variable
  Y <- as.vector(as.matrix(Y))
  if (all(Y == 0 | Y == 1)) {
    family$family <- "binomial"
  } else {
    family$family <- "gaussian"
  }
  if (ncol(X) > 8) {
    cat("- screen.randomForest (ntree =", ntree, ") was started and screens", nVar, " vars. -\n")
    t_ime <- system.time({
      SuperLearner:::.SL.require("randomForest")
      if (family$family == "gaussian") {
        rank.rf.fit <- randomForest::randomForest(Y ~ .,
          data = X,
          ntree = ntree, mtry = mtry, nodesize = nodesize,
          keep.forest = FALSE, maxnodes = maxnodes, importance = TRUE
        )
        # variables with the largest %IncMSE are the most important ones
        # negative scores mean zero or low importance
        imp_measure <- rank.rf.fit$importance[, "%IncMSE"]
        return(rank(-imp_measure, ties.method = "random") <= nVar)
      }
      if (family$family == "binomial") {
        rank.rf.fit <- randomForest::randomForest(as.factor(Y) ~ .,
          data = X,
          ntree = ntree, mtry = mtry, nodesize = nodesize,
          keep.forest = FALSE, maxnodes = maxnodes, importance = TRUE
        )
        # variables with the largest mean decrease in accuracy are the most important ones
        # negative scores mean zero or low importance
        imp_measure <- rank.rf.fit$importance[, "MeanDecreaseAccuracy"]
        return(rank(-imp_measure, ties.method = "random") <= nVar)
      }
    })
    cat("- screen.randomForest finished - \n")
  }
  rep(TRUE, ncol(X))
}

# Try different Hyperparamters
tuneGrid <- expand.grid(ntree = c(500, 1000))
for (i in seq(nrow(tuneGrid))) {
  eval(parse(text = paste0(
    "screen.randomForest_grid", tuneGrid[i, 1],
    "<- function(..., ntree = ", tuneGrid[i, 1], ")
    {screen.randomForest_base(..., ntree = ntree)}"
  )))
}

### 2) Cramers V with vcd::assocstats

screen.cramersv_base <- function(Y, X, nscreen = 4, num_cat = 10, ...) {
  cat("- screen.cramersv screens", nscreen, "variables -\n")
  require(vcd)

  # Selects nscreen variables among X that have the highest association with Y.
  # Metric vars (X or Y) are automatically turned into factors when the number
  # of unique observations falls below num_cat.

  if (ncol(X) > 8) {
    dat <- cbind(Y, X)
    contin_var <- apply(dat, 2, function(var) length(unique(var)) > num_cat)
    make_categ <- function(var) {
      cut(var, unique(quantile(var, prob = seq(0, 1, 0.2))), include.lowest = T)
    }
    if (any(contin_var)) {
      dat[, contin_var] <- apply(dat[, contin_var, drop = FALSE], 2, make_categ)
    }
    calc_cram_v <- function(x_var, y_var) assocstats(table(y_var, x_var))$cramer
    cramers_v <- apply(dat[, !colnames(dat) %in% "Y"], 2, calc_cram_v, y_var = dat[, "Y"])
    # cramers_v is between 0 and 1, the higher the more relevant
    return(unname(rank(-cramers_v, ties.method = "random") <= nscreen))
  }
  rep(TRUE, ncol(X))
}

# Try different Hyperparamters
tuneGrid <- expand.grid(nscreen = seq(4, 8, 4))
for (i in seq(nrow(tuneGrid))) {
  eval(parse(text = paste0(
    "screen.cramersv_grid", tuneGrid[i, 1],
    "<- function(..., nscreen = ", tuneGrid[i, 1], ")
                           {screen.cramersv_base(..., nscreen = nscreen)}"
  )))
}

### 3) Elastic Net with glmnet::cv.glmnet

screen.glmnet_nVar <- function(Y, X, family = list(), alpha = 0.75, nfolds = 5, nlambda = 150, nVar = 8, ...) {
  SuperLearner:::.SL.require("glmnet")

  # relevant for column names but shouldnt be a matrix anyways
  X <- as.data.frame(X)

  # chose family dependent upon response
  Y <- as.vector(as.matrix(Y))
  if (all(Y == 0 | Y == 1)) {
    family$family <- "binomial"
  } else {
    family$family <- "gaussian"
  }

  # needed for var names to select from levels of factors later on
  if (ncol(X) > 26 * 27) stop("Find further column names for X!")
  let <- c(letters, sort(do.call("paste0", expand.grid(letters, letters[1:26]))))
  names(X) <- let[1:ncol(X)]

  # factors are coded as dummies which are standardized in cv.glmnet()
  # intercept is not in model.matrix() because its already in cv.glmnet()
  is_fact_var <- sapply(X, is.factor)
  X <- try(model.matrix(~ -1 + ., data = X), silent = FALSE)

  # cv.glmnet() calls glmnet(), thus arguments are given to glmnet()
  if (ncol(X) > 8) {
    fitCV <- try(glmnet::cv.glmnet(
      x = X, y = Y, lambda = NULL, type.measure = "deviance",
      nfolds = nfolds, family = family$family, alpha = alpha,
      nlambda = nlambda, keep = T
    ), silent = TRUE)
    # if no variable was selected, penalization might have been too strong, try log(lambda)
    if (all(fitCV$nzero == 0) | all(is.na(fitCV$nzero))) {
      fitCV <- try(glmnet::cv.glmnet(
        x = X, y = Y, lambda = log(fitCV$glmnet.fit$lambda + 1), type.measure = "deviance",
        nfolds = nfolds, family = family$family, alpha = alpha, keep = T
      ), silent = TRUE)
    }
    # if nVar is not available, take the closest to nVar available
    if (all(fitCV$nzero != nVar)) {
      lambda_index_with_nVar <- min(which(abs(fitCV$nzero - nVar) == min(abs(fitCV$nzero - nVar))))
      # nVar is available
    } else if (any(fitCV$nzero == nVar)) {
      lambda_index_with_nVar <- min(which(fitCV$nzero == nVar))
    }
    coefs <- coef(fitCV$glmnet.fit, s = fitCV$glmnet.fit$lambda[lambda_index_with_nVar])
    var_nms <- coefs@Dimnames[[1]]

    # Instead of Group Lasso:
    # If any level of a dummy coded factor is selected, the whole factor is selected
    if (any(is_fact_var)) {
      nms_fac <- names(which(is_fact_var))
      is_selected <- coefs[-1] != 0 # drop intercept
      # model.matrix adds numbers to dummy coded factors which we need to get rid of
      var_nms_sel <- gsub("[^::a-z::]", "", var_nms[-1][is_selected])
      sel_fac <- nms_fac[nms_fac %in% var_nms_sel]
      sel_numer <- var_nms_sel[!var_nms_sel %in% sel_fac]
      all_sel_vars <- c(sel_fac, sel_numer)
      whichVariable <- names(is_fact_var) %in% all_sel_vars
    } else {
      # metric variables only
      whichVariable <- coefs[-1] != 0
    }
  } else {
    # no screening
    whichVariable <- rep(TRUE, ncol(X))
  }
  if (nVar != sum(whichVariable)) {
    cat("- Screening with glmnet:", nVar, "vars were not available, instead", sum(whichVariable), "vars were screened. - \n")
  }
  return(whichVariable)
}

### 4) Pearson Correlation Coef. with cor.test()$p.value

screen.corPearson <- function(Y, X, family, obsWeights, id, method = "pearson", minPvalue = 0.01,
                              minscreen = 3, maxscreen = 8, ...) {
  if (ncol(X) > 8) {
    p_val <- apply(X, 2, function(var) {
      # factors are only selected when no driving metrics seem to be there
      if (length(unique(var)) < 5 | is.factor(var)) {
        return(minPvalue + 1e-4)
      }
      # predictors with zero variance are not meaningful
      if (var(var) == 0) {
        return(1)
      }
      cor.test(var, y = Y, method = "pearson")$p.value
    })
    no_sel <- sum(p_val <= minPvalue)
    if (no_sel > maxscreen) {
      return(rank(p_val, ties.method = "random") <= maxscreen)
    }
    if (no_sel < minscreen) {
      return(rank(p_val, ties.method = "random") <= minscreen)
    }
    return(p_val <= minPvalue)
  }
  rep(TRUE, ncol(X))
}

#### Prediction ####

# 1) multivariate adaptive regression splines with 'earth'

SL.earth2 <- function(Y, X, newX = NULL, family = list(), obsWeights = NULL, id = NULL, degree = 2, penalty = 3,
                      nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", nfold = 0,
                      ncross = 1, minspan = 0, endspan = 0, ...) {
  cat(" - earth2 was started and is making predictions - \n")

  # chose family dependent upon response
  Y <- as.vector(as.matrix(Y))
  if (all(Y == 0 | Y == 1)) {
    family$family <- "binomial"
  } else {
    family$family <- "gaussian"
  }

  earth_time <- system.time({
    SuperLearner:::.SL.require("earth")
    if (family$family == "gaussian") {
      fit.earth <- earth::earth(
        x = X, y = Y, degree = degree,
        nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold,
        ncross = ncross, minspan = minspan, endspan = endspan
      )
    }
    if (family$family == "binomial" & all(Y %in% c(0, 1))) {
      fit.earth <- earth::earth(
        x = X, y = Y, degree = degree,
        nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold,
        ncross = ncross, minspan = minspan, endspan = endspan,
        glm = list(family = binomial)
      )
    }
    if (family$family == "binomial" & all(Y %in% c(0, 1)) == FALSE) {
      fit.earth <- earth::earth(
        x = X, y = Y, degree = degree,
        nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold,
        ncross = ncross, minspan = minspan, endspan = endspan
      )
    }
    pred <- predict(fit.earth, newdata = newX, type = "response")
    fit <- list(object = fit.earth)
    out <- list(pred = pred, fit = fit)
    class(out$fit) <- c("SL.earth")
  })
  cat("- earth2 was finished lasting:", round(unclass(earth_time)["elapsed"], 2), " - \n")
  return(out)
}

# 2) Random Forest with different numbers of Trees

SL.randomForest_base <- function(Y, X, newX = NULL, family = list(), mtry = ifelse(family$family ==
                                   "gaussian", max(floor(ncol(X) / 3), 1), floor(sqrt(ncol(X)))),
                                 ntree = 100, nodesize = ifelse(family$family == "gaussian",
                                   5, 1
                                 ), maxnodes = NULL, importance = FALSE, ...) {
  cat(" - randomForest was started (ntree =", ntree, ") and is making predictions - \n")
  randomForest_time <- system.time({
    SuperLearner:::.SL.require("randomForest")

    # avoid infinite search for split points in trees
    if (all(apply(X,2,var) == 0)) {
      fit.rf <- "Empty"
      attr(fit.rf, "class") <- "try-error"
      pred <- rep(mean(Y), nrow(Xnew))
      fit <- list(object = fit.rf)
      cat("- Failed random forest - \n")
    }
    
    if (family$family == "gaussian" & !exists("fit.rf")) {
      fit.rf <- randomForest::randomForest(Y ~ .,
        data = X,
        ntree = ntree, xtest = newX, keep.forest = TRUE,
        mtry = mtry, nodesize = nodesize, maxnodes = maxnodes,
        importance = importance
      )
      try(pred <- fit.rf$test$predicted, silent = TRUE)
      if (any(class(fit.rf) == "try-error")) {
        pred <- rep(mean(Y), nrow(Xnew))
        cat("- Failed random forest - \n")
      }
      fit <- list(object = fit.rf)
    }
    if (family$family == "binomial" & !exists("fit.rf")) {
      fit.rf <- randomForest::randomForest(
        y = as.factor(Y),
        x = X, ntree = ntree, xtest = newX, keep.forest = TRUE,
        mtry = mtry, nodesize = nodesize, maxnodes = maxnodes,
        importance = importance
      )
      try(pred <- fit.rf$test$votes[, 2], silent = TRUE)
      if (any(class(fit.rf) == "try-error")) {
        pred <- rep(mean(Y), nrow(Xnew))
        cat("- Failed random forest - \n")
      }
      fit <- list(object = fit.rf)
    }
    out <- list(pred = pred, fit = fit)
    class(out$fit) <- c("SL.randomForest")
  })
  cat("- randomForest was finished lasting: ", round(unclass(randomForest_time)["elapsed"], 2), " - \n")
  return(out)
}

# Try different Hyperparamters
tuneGrid <- expand.grid(ntree = c(500, 1000))
for (i in seq(nrow(tuneGrid))) {
  eval(parse(text = paste0(
    "SL.randomForest_grid", tuneGrid[i, 1],
    "<- function(..., ntree = ", tuneGrid[i, 1], ")
                           {SL.randomForest_base(..., ntree = ntree)}"
  )))
}

# 3) Generalized Boosted Regression with tuning grid
# failes frequently and takes a lot of time

SL.gbm_base <- function(Y, X, newX = NULL, family = list(), obsWeights = NULL, gbm.trees = 10,
                        interaction.depth = 1, shrinkage = 0.001, ...) {
  cat(" - GBM started and is making predictions, interaction depth = ", interaction.depth, " and gbm.trees = ", gbm.trees, " - \n")
  SuperLearner:::.SL.require("gbm")
  gbm.model <- as.formula(paste("Y~", paste(colnames(X), collapse = "+")))

  # chose family dependent upon response
  Y <- as.vector(as.matrix(Y))
  if (all(Y == 0 | Y == 1)) {
    family$family <- "binomial"
  } else {
    family$family <- "gaussian"
  }

  if (family$family == "gaussian") {
    fit.gbm <- gbm::gbm(
      formula = gbm.model, data = X, distribution = "gaussian",
      n.trees = gbm.trees, interaction.depth = interaction.depth,
      shrinkage = shrinkage, cv.folds = 5, keep.data = TRUE, bag.fraction = 0.75,
      weights = obsWeights, verbose = FALSE
    )
  }
  if (family$family == "binomial" & all(Y %in% c(0, 1))) {
    fit.gbm <- gbm::gbm(
      formula = gbm.model, data = X, distribution = "bernoulli",
      n.trees = gbm.trees, interaction.depth = interaction.depth,
      shrinkage = shrinkage, cv.folds = 5, keep.data = TRUE, bag.fraction = 0.75,
      weights = obsWeights, verbose = FALSE
    )
  }
  if (family$family == "binomial" & all(Y %in% c(0, 1)) == FALSE) {
    fit.gbm <- gbm::gbm(
      formula = gbm.model, data = X, distribution = "gaussian",
      n.trees = gbm.trees, interaction.depth = interaction.depth,
      shrinkage = shrinkage, cv.folds = 5, keep.data = TRUE, bag.fraction = 0.75,
      weights = obsWeights, verbose = FALSE
    )
  }
  best.iter <- gbm::gbm.perf(fit.gbm, method = "cv", plot.it = FALSE)
  pred <- try(predict(fit.gbm, newdata = newX, best.iter, type = "response"), silent = TRUE)
  if (any(class(pred) == "try-error")) {
    cat("- GBM Unsuccessful! - \n")
    pred <- rep(median(Y), length(Y))
  }
  fit <- list(object = fit.gbm, n.trees = best.iter)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.gbm")
  cat(" - GBM finished - \n")
  return(out)
}

# Try different Hyperparamters
tuneGrid <- expand.grid(gbm.trees = c(10), interaction.depth = c(1, 2))
for (i in seq(nrow(tuneGrid))) {
  eval(parse(text = paste0(
    "SL.gbm_grid", tuneGrid[i, 1], "_interaction.depth", tuneGrid[i, 2],
    "<- function(..., gbm.trees = ", tuneGrid[i, 1], ", interaction.depth = ", tuneGrid[i, 2], ")
                           {SL.gbm_base(..., gbm.trees = gbm.trees, interaction.depth = interaction.depth)}"
  )))
}

# 4) GAM from 'gam' package just slighlty modified for fluent estimation procedure
# no substantial changes

SL.gam2 <- function(Y, X, newX = NULL, family = list(), obsWeights = NULL, deg.gam = 2, cts.num = 5,
                    ...) {

  # deg.gam == 2 is needed due to frequently occurring convergence failure
  # chose family dependent upon response variable
  Y <- as.vector(as.matrix(Y))
  if (all(Y == 0 | Y == 1)) {
    family$family <- binomial()
  } else {
    family$family <- gaussian() # scat() also reasonable but computationally intense
  }

  cat(" - gam::gam was started and is making predictions - \n")
  gam_time <- system.time({
    SuperLearner:::.SL.require("gam")
    s <- gam:::s # s() is also used by 'mgcv' package - avoid clash

    # adjust model formula for metric and categorical predictors
    metric_var <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
    if (sum(metric_var) != 0 & sum(metric_var) != length(metric_var)) {
      # metric and categorical variables
      gam.model <- as.formula(paste("Y~", paste(paste("s(",
        colnames(X[, metric_var, drop = FALSE]), ",", deg.gam,
        ")",
        sep = ""
      ), collapse = "+"), "+", paste(colnames(X[,
        !metric_var,
        drop = FALSE
      ]), collapse = "+")))
    }
    if (all(metric_var)) {
      # metric variables only
      gam.model <- as.formula(paste("Y~", paste(paste("s(",
        colnames(X[, metric_var, drop = FALSE]), ",", deg.gam,
        ")",
        sep = ""
      ), collapse = "+")))
    } else {
      # all categorical
      gam.model <- as.formula(paste("Y~", paste(colnames(X),
        collapse = "+"
      ), sep = ""))
    }
    fit.gam <- gam::gam(gam.model,
      data = X, family = family$family,
      control = gam::gam.control(maxit = 50, bf.maxit = 50),
      weights = obsWeights
    )
    # or predict.gam depending on version
    pred <- gam::predict.Gam(fit.gam, newdata = newX, type = "response")
    fit <- list(object = fit.gam)
    out <- list(pred = pred, fit = fit)
    class(out$fit) <- c("SL.gam")
  })
  cat("- gam::gam was finished lasting: ", round(unclass(gam_time)["elapsed"], 2), " - \n")
  return(out)
}

# 6) glm.interaction with informative output

SL.glm.interaction_info <- function(Y, X, newX = NULL, family = list(), obsWeights = NULL, ...) {
  cat("- GLM Interaction was started and is making predictions - \n")
  glm_time <- system.time({

    # chose family dependent upon response variable
    Y <- as.vector(as.matrix(Y))
    if (all(Y == 0 | Y == 1)) {
      family$family <- binomial()
    } else {
      family$family <- gaussian() # scat() also reasonable but computationally intense
    }

    if (is.matrix(X)) {
      X <- as.data.frame(X)
    }
    fit.glm <- glm(Y ~ .^2, data = X, family = family$family, weights = obsWeights)
    if (is.matrix(newX)) {
      newX <- as.data.frame(newX)
    }
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm"
    out <- list(pred = pred, fit = fit)
  })

  cat("- GLM Interaction was finished lasting: ", round(unclass(glm_time)["elapsed"], 2), " - \n")
  return(out)
}

# Diagnostic Function for the Clever Covariate

cc_trunc <- function(ltmle_est){
  
  # input: ltmle output (list), two (A1 and A0) binary matrices (same dimensions as 'abar' from ltmle)
  # output: matrix with share of truncated prob. and a default summary() statistics for the cc
  
  # Function that calculates a summary statistics of the "clever" covariate (cc)
  # and the trunacted share of intervention probabilities/propensity scores 
  # from an ltmle output. The summary statistics of the cc is for the last point in time.
  
  # for treatment: A1 - Dimension of cum.g/cum.g.unbounded: [observations,time-points,treatment]
  cum_prob_A1_used <- ltmle_est$cum.g[,ncol(ltmle_est$cum.g),1] # last column = last point in time
  cum_prob_A1_est  <- ltmle_est$cum.g.unbounded[,ncol(ltmle_est$cum.g),1] # last column = last point in time
  followed_A1 <- ltmle_est$cum.g.used[,ncol(ltmle_est$cum.g),1] # which subjects acutally followed A1, hence are used for the updating step (since all are uncensored here)
  cc_A1 <- 1/cum_prob_A1_used # inverse probabilites 
  cc_A1[!followed_A1] <- 0 # those who did NOT follow A1 in every point in time til the very end are set to 0
  cc_A1_stat <- unclass(summary(cc_A1)) # summary statistics for the clever covariate
  
  # how many probabilities were truncated (among those who followed A1 and are uncensored)
  trunc_share_A1 <- mean(cum_prob_A1_est[followed_A1]!=cum_prob_A1_used[followed_A1])
  
  # for control: A0 - Dimension of cum.g/cum.g.unbounded: [observations,time-points,control]
  cum_prob_A0_used <- ltmle_est$cum.g[,ncol(ltmle_est$cum.g),2] # last column = last point in time
  cum_prob_A0_est  <- ltmle_est$cum.g.unbounded[,ncol(ltmle_est$cum.g),2] # last column = last point in time
  followed_A0 <- ltmle_est$cum.g.used[,ncol(ltmle_est$cum.g),2] # which subjects acutally followed A0, hence are used for the updating step (since all are uncensored here)
  cc_A0 <- 1/cum_prob_A0_used # inverse probabilites 
  cc_A0[!followed_A0] <- 0 # those who did NOT follow A0 in every point in time til the very end are set to 0
  cc_A0_stat <- unclass(summary(cc_A0)) # summary statistics for the clever covariate
  
  # how many probabilites were truncated (among those who followed A0 and are uncensored)
  trunc_share_A0 <- mean(cum_prob_A0_est[followed_A0]!=cum_prob_A0_used[followed_A0])
  
  # collect results for output
  cc_trunc_matrix <- matrix(c(cc_A0_stat,trunc_share_A0,cc_A1_stat,trunc_share_A1),7,2)
  rownames(cc_trunc_matrix) <- c(names(cc_A0_stat),"TruncShare")
  colnames(cc_trunc_matrix) <- c("Control","Treatment")
  
  cc_info <- list(cc_trunc_matrix=cc_trunc_matrix,dist_cc_A1=cc_A1,dist_cc_A0=cc_A0)
  
  return(cc_info)
}