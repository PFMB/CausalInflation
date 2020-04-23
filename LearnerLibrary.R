######
# library of wrapper functions/learner for variable selection and estimation/prediction
# for use in SuperLearner/LTMLE
######

#### Screening ####

# 1) Random Forest with 'randomForest::randomForest'
# -- Hyperparamter Tuning via grid search and 'caret' package

screen.randomForest_base <- function(Y, X, family = list(), nVar = 5, ntree = ntree, mtry = ifelse(family$family ==
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
        # X <- matrix(rnorm(12 * 1000), ncol = 1000)
        # Y <- t(c(5, 6, -4, -3, 2, -0.001, 0.001, 0.002, 4, 7, -2, 5) %*% X + rnorm(1000))
        
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
        # X <- matrix(rnorm(12 * 1000), ncol = 1000)
        # Y <- rbinom(1000, 1, plogis(c(5, 6, -4, -3, 2, -0.001, 0.001, 0.002, 4, 7, -2, 5) %*% X + rnorm(1000)))
        
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
  rep(TRUE, ncol(X)) # select all if less than 11 vars are contained in X
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

# 2) Cramers V with vcd::assocstats
# -- Hyperparamter Tuning via grid search and 'caret' package

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

# 3) Elastic Net with glmnet::cv.glmnet

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

# 4) Pearson Correlation Coef. with cor.test()$p.value

screen.corPearson <- function(Y, X, family, obsWeights, id, method = "pearson", minPvalue = 0.01,
                              minscreen = 2, maxscreen = 8, ...) {
  if (ncol(X) > 8) {
    p_val <- sapply(X, function(var) {
      # factors are only selected when no driving metrics seem to be there
      if (length(unique(var)) < 5 | is.factor(var)) return(minPvalue + 1e-4)
      # predictors with zero variance are not meaningful
      if (var(var) == 0) return(1)
      cor.test(var, y = Y, method = "pearson")$p.value
    })
    no_sel <- sum(p_val <= minPvalue)
    if (no_sel > maxscreen) return(rank(p_val, ties.method = "random") <= maxscreen)
    if (no_sel < minscreen) return(rank(p_val, ties.method = "random") <= minscreen)
    return(p_val <= minPvalue)
  }
  rep(TRUE, ncol(X))
}

#### Prediction ####

# 1) Multivariate Adaptive Regression Splines with 'earth'

SL.earth2 <- function(Y, X, newX, family, obsWeights, id, degree = 2, penalty = 3,
                      nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", nfold = 0,
                      ncross = 1, minspan = 0, endspan = 0, ...) {
  cat(" - earth2 was started and is making predictions - \n")
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

SL.randomForest_base <- function(Y, X, newX, family, mtry = ifelse(family$family ==
                                   "gaussian", max(floor(ncol(X) / 3), 1), floor(sqrt(ncol(X)))),
                                 ntree = ntree, nodesize = ifelse(family$family == "gaussian",
                                   5, 1
                                 ), maxnodes = NULL, importance = FALSE, ...) {
  cat(" - randomForest was started (ntree =", ntree, ") and is making predictions - ")
  randomForest_time <- system.time({
    SuperLearner:::.SL.require("randomForest")

    if (family$family == "gaussian") {
      fit.rf <- randomForest::randomForest(Y ~ .,
        data = X,
        ntree = ntree, xtest = newX, keep.forest = TRUE,
        mtry = mtry, nodesize = nodesize, maxnodes = maxnodes,
        importance = importance
      )
      try(pred <- fit.rf$test$predicted, silent = TRUE)
      if (any(class(fit.rf) == "try-error")) {
        pred <- rep(mean(Y), length(Y))
        cat(paste0("FAIL!"))
      }
      fit <- list(object = fit.rf)
    }
    if (family$family == "binomial") {
      fit.rf <- randomForest::randomForest(
        y = as.factor(Y),
        x = X, ntree = ntree, xtest = newX, keep.forest = TRUE,
        mtry = mtry, nodesize = nodesize, maxnodes = maxnodes,
        importance = importance
      )
      try(pred <- fit.rf$test$votes[, 2], silent = TRUE)
      if (any(class(fit.rf) == "try-error")) {
        pred <- rep(mean(Y), length(Y))
        cat(paste0("FAIL!"))
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

SL.gbm_base <- function(Y, X, newX, family, obsWeights, gbm.trees = gbm.trees,
                        interaction.depth = interaction.depth, shrinkage = 0.001, ...) {
  cat(" - GBM started and is making predictions, interaction depth = ", interaction.depth, " and gbm.trees = ", gbm.trees, " - \n")
  SuperLearner:::.SL.require("gbm")
  gbm.model <- as.formula(paste("Y~", paste(colnames(X), collapse = "+")))

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
    cat("GBM Unsuccessful!")
    pred <- rep(median(Y), length(Y))
  }
  fit <- list(object = fit.gbm, n.trees = best.iter)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.gbm")
  cat(" - GBM finished - ")
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

# 5) GAM from 'gam' package just slighlty modified for fluent estimation procedure
# no substantial changes

SL.gam2 <- function(Y, X, newX, family, obsWeights, deg.gam = 2, cts.num = 4,
                    ...) {

  # chose family dependent upon response variable
  Y <- as.vector(as.matrix(Y))
  if (all(Y == 0 | Y == 1)) {
    family$family <- binomial()
  } else {
    family$family <- gaussian() # scat() also reasonable but computationally intense
  }

  cat(" - GAM2 was started and is making predictions - \n")
  gam_time <- system.time({
    SuperLearner:::.SL.require("gam")
    s <- gam:::s # s() is also used by 'mgcv' package - avoid clash

    cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
    if (sum(!cts.x) > 0) {
      gam.model <- as.formula(paste("Y~", paste(paste("s(",
        colnames(X[, cts.x, drop = FALSE]), ",", deg.gam,
        ")",
        sep = ""
      ), collapse = "+"), "+", paste(colnames(X[,
        !cts.x,
        drop = FALSE
      ]), collapse = "+")))
    }
    else {
      gam.model <- as.formula(paste("Y~", paste(paste("s(",
        colnames(X[, cts.x, drop = FALSE]), ",", deg.gam,
        ")",
        sep = ""
      ), collapse = "+")))
    }
    if (sum(!cts.x) == length(cts.x)) {
      gam.model <- as.formula(paste("Y~", paste(colnames(X),
        collapse = "+"
      ), sep = ""))
    }
    fit.gam <- gam::gam(gam.model,
      data = X, family = family,
      control = gam::gam.control(maxit = 50, bf.maxit = 50),
      weights = obsWeights
    )
    pred <- gam::predict.Gam(fit.gam, newdata = newX, type = "response") # or predict.gam depending on version
    fit <- list(object = fit.gam)
    out <- list(pred = pred, fit = fit)
    class(out$fit) <- c("SL.gam")
  })
  cat("- GAM2 was finished lasting: ", round(unclass(gam_time)["elapsed"], 2), " - \n")
  return(out)
}

# 6) GLM.Interaction with informative output
SL.glm.interaction_info <- function(Y, X, newX, family, obsWeights, ...) {
  cat("- GLM Interaction was started and is making predictions - \n")
  glm_time <- system.time({
    if (is.matrix(X)) {
      X <- as.data.frame(X)
    }
    fit.glm <- glm(Y ~ .^2, data = X, family = family, weights = obsWeights)
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