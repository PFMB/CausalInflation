library(stringr) # remove quotationmarks with str_replace_all()
library(vcd) #
library(magrittr) # from factors to numeric
# libarary of wrapper functions/learner for variable selection and estimation/prediction
# for use in SuperLearner/LTMLE

#### Screening ####

# 1) Random Forest with 'randomForest'
# -- Hyperparamter Tuning via grid search and 'caret' package

screen.randomForest_base <- function(Y, X, family, nVar = 10, ntree = ntree, mtry = ifelse(family$family ==
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
  cat(paste0(" - screen.randomForest (ntree =", ntree, ") was started and screens - "))
  if (ncol(X) > 10) {
    screen.randomForest_time <- system.time({
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
        (whichVariable <- (rank(-(imp_measure)) <= nVar))
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
        (whichVariable <- (rank(-(imp_measure)) <= nVar))
      }
    })
    cat("- screen.randomForest finished lasting: ", round(unclass(screen.randomForest_time)["elapsed"], 2), " - ")
  } else {
    whichVariable <- rep(TRUE, ncol(X)) # select all
  }
  return(whichVariable)
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
# search for algorithms in the enviroment
all_objects <- ls()
my_randomForest_screen_objects <- grep("screen.randomForest_grid", all_objects)
all_randomForests_screen <- all_objects[my_randomForest_screen_objects]
rm(tuneGrid)
rm(all_objects)


# 2) Different correlation coeeficients with 'cor.test'
# -- "pearson", "kendall", or "spearman"

screen.corRank_base <- function(Y, X, method = method, rank = 2, ...) {
  cat("-screen.corRank was started and screens with method = ", method, "- \n")
  screen.corRank_time <- system.time({
    # from factor to numeric
    factor_to_numeric <- function(x) {
      if (is.factor(x)) {
        x <- unclass(x) %>% as.numeric()
      } else {
        x <- x
      }
    }
    X <- sapply(X, factor_to_numeric)

    listp <- apply(X, 2, function(x, Y, method) {
      ifelse(var(x) <= 0, 1, cor.test(x, y = Y, method = method)$p.value)
    }, Y = Y, method = method)
    whichVariable <- (rank(listp) <= rank)
  })
  cat("-screen.corRank finished lasting:", round(unclass(screen.corRank_time)["elapsed"], 2), "-\n")
  return(whichVariable)
}

# Try different Hyperparamters
tuneGrid <- expand.grid(methods = paste0('"', c("pearson", "kendall", "spearman"), '"'))
for (i in seq(nrow(tuneGrid))) {
  eval(parse(text = paste0(
    "screen.corRank_grid", str_replace_all(tuneGrid[i, 1], "[[:punct:]]", ""),
    "<- function(..., method = ", tuneGrid[i, 1], ")
                           {screen.corRank_base(..., method = method)}"
  )))
}
# search for algorithms in the enviroment
all_objects <- ls()
my_corRank_objects <- grep("screen.corRank_grid", all_objects)
all_corRank <- all_objects[my_corRank_objects]
rm(tuneGrid)
rm(all_objects)


# 3) Cramers V
screen.cramersv_base <- function(Y, X, nscreen = 4, num_cat = 10, ...) {
  cat("- screen.cramersv screens", nscreen, "variables -\n")
  require(vcd)

  # Selects nscreen variables among X that have the highest association with Y.
  # Metric vars (X or Y) are automatically turned into factors when the number
  # of unique observations falls below num_cat.

  if (ncol(X) > nscreen) {
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
    return(unname(rank(-cramers_v) <= nscreen))
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
# search for algorithms in the enviroment
all_objects <- ls()
my_cramersv_objects <- grep("screen.cramersv_grid", all_objects)
all_cramersv <- all_objects[my_cramersv_objects]
rm(all_objects)


# 4) Screen variables with the elastic net from 'glmnet' package
screen.glmnet_base <- function(Y, X, family, alpha = alpha, minscreen = 2, pw = T,
                               maxtries = 4, nfolds = 10, nlambda = nlambda, ...) {
  cat(paste0(" - screen.glmnet was started and screens with alpha = ", alpha, " and nlambda = ", nlambda, " - "))
  screen.glmnet_time <- system.time({
    if (family$family == "binomial" & all(Y %in% c(0, 1)) == FALSE) {
      myfamily <- "gaussian"
    } else {
      myfamily <- family$family
    }
    if (!is.matrix(X)) {
      X <- try(model.matrix(~ -1 + ., data = X), silent = T)
    }
    successfulfit <- FALSE
    fitCV <- try(glmnet::cv.glmnet(
      x = X, y = Y, lambda = NULL, type.measure = "deviance",
      nfolds = nfolds, family = myfamily, alpha = alpha,
      nlambda = nlambda, keep = T
    ), silent = T)
    if (class(fitCV) == "try-error") {
      i <- 2
      while (successfulfit == FALSE & i <= maxtries) {
        if (pw == T) {
          cat(paste("glmnet failed, new try #", i, "\n"))
        }
        fitCV <- try(glmnet::cv.glmnet(
          x = X, y = Y, lambda = NULL, type.measure = "deviance",
          nfolds = nfolds, family = myfamily, alpha = alpha, nfolds = 4,
          nlambda = nlambda * (i + 3), keep = T, foldid = sample(fitCV$foldid)
        ), silent = T)
        i <- i + 1
        if (class(fitCV) == "try-error") {
          successfulfit <- FALSE
        } else {
          successfulfit <- TRUE
        }
      }
    } else {
      successfulfit <- TRUE
    }
    whichVariable <- NULL
    if (successfulfit == TRUE) {
      si <- try(abs((max(fitCV$nzero) - (dim(X)[2] / 2)) - fitCV$nzero), silent = T)
      whichVariable2 <- try((as.numeric(coef(fitCV$glmnet.fit, s = fitCV$lambda.min))[-1] != 0), silent = T)
      whichVariable3 <- try(as.numeric(glmnet::coef.glmnet(fitCV, s = fitCV$lambda[si == min(si)][1]))[-1] != 0, silent = T)
      if (sum(whichVariable2) > 1) {
        whichVariable <- whichVariable2
      } else {
        if (sum(whichVariable3) > 1 & sum(whichVariable3) < (dim(X)[2] / 2)) {
          whichVariable <- whichVariable3
        }
      }
    }
    if (is.null(whichVariable)) {
      whichVariable <- c(T, rep(F, dim(X)[2] - 1))
      if (pw == T) {
        cat("Lasso failed and only first variable was screened \n")
      }
    }
    if (pw == T) {
      cat("Number of included variables:", sum(whichVariable), "\n")
    }
  })
  cat("- screen.glmnet was finished lasting: ", round(unclass(screen.glmnet_time)["elapsed"], 2), " - ")
  return(whichVariable)
}

# Try different Hyperparamters
tuneGrid <- expand.grid(alpha = seq(0.5, 1, 0.25), nlambda = c(300, 400))
for (i in seq(nrow(tuneGrid))) {
  eval(parse(text = paste0(
    "screen.glmnet_grid", tuneGrid[i, 1], "_nlambda", tuneGrid[i, 2],
    "<- function(..., alpha = ", tuneGrid[i, 1], ", nlambda = ", tuneGrid[i, 2], ")
                           {screen.glmnet_base(..., alpha = alpha, nlambda = nlambda)}"
  )))
}
# search for algorithms in the enviroment
all_objects <- ls()
my_glmnet_objects <- grep("screen.glmnet_grid", all_objects)
all_glmnets <- all_objects[my_glmnet_objects]
rm(tuneGrid)
rm(all_objects)

# GLMNET with nVar

screen.glmnet_nVar <- function(Y, X, family, alpha = 0.75, nfolds = 5, nlambda = 200, nVar = 8, ...) {
  SuperLearner:::.SL.require("glmnet")

  # chose family dependent upon response variable
  Y <- as.vector(as.matrix(Y))
  if (all(Y == 0 | Y == 1)) {
    family$family <- "binomial"
  } else {
    family$family <- "gaussian"
  }
  if (!is.matrix(X)) {
    X <- try(model.matrix(~ -1 + ., data = X), silent = T)
  }
  if (ncol(X) > 1) { # if only one variable, no screening is conducted
    whichVariable <- rep(FALSE, ncol(X)) # no intercept
    fitCV <- try(glmnet::cv.glmnet(
      x = X, y = Y, lambda = NULL, type.measure = "deviance",
      nfolds = nfolds, family = family$family, alpha = alpha,
      nlambda = nlambda, keep = T
    ), silent = TRUE)
    if (all(fitCV$nzero == 0) | all(is.na(fitCV$nzero))) { # failed -> transform lamdas from prior fit and try again
      fitCV <- try(glmnet::cv.glmnet(
        x = X, y = Y, lambda = log(fitCV$glmnet.fit$lambda + 1), type.measure = "deviance",
        nfolds = nfolds, family = family$family, alpha = alpha, keep = T
      ), silent = TRUE)
    }
    if (all(fitCV$nzero != nVar)) { # nVar is not available, take the closest to nVar available
      lambda_index_with_nVar <- min(which(abs(fitCV$nzero - nVar) == min(abs(fitCV$nzero - nVar))))
    } else if (any(fitCV$nzero == nVar)) { # lambda value for the no. of coef. as defined in nVar
      lambda_index_with_nVar <- min(which(fitCV$nzero == nVar))
    }
    coefList <- coef(fitCV$glmnet.fit, s = fitCV$glmnet.fit$lambda[lambda_index_with_nVar])
    non_zero_index <- which(coefList != 0)[-1] - 1 # drop intercept and adjust index reducing one
    whichVariable[non_zero_index] <- TRUE
  } else {
    whichVariable <- rep(TRUE, ncol(X)) # no screening
  }
  if (nVar != sum(whichVariable)) {
    cat("nVar (", nVar, ") is not available, instead ", sum(whichVariable), " vars were screened.")
  }
  return(whichVariable)
}

######

screen.corPearson <- function(Y, X, family, obsWeights, id, method = "pearson", minPvalue = 0.01,
                              minscreen = 2, ...) {
  if (ncol(X) > 5) {
    listp <- apply(X, 2, function(x, Y, method) {
      ifelse(var(x) <= 0, 1, cor.test(x, y = Y, method = method)$p.value)
    }, Y = Y, method = method)
    whichVariable <- (listp <= minPvalue)
    if (sum(whichVariable) < minscreen) {
      cat("number of variables with p value less than minPvalue is less than minscreen")
      whichVariable[rank(listp) <= minscreen] <- TRUE
    }
  } else {
    whichVariable <- rep(TRUE, ncol(X))
  }
  cat("CorPearson screens ", sum(whichVariable), " variables.")
  return(whichVariable)
}

#### ----------- Prediction ----------- ####

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
  cat("- earth2 was finished lasting: ", round(unclass(earth_time)["elapsed"], 2), " - \n")
  return(out)
}


# 2) Random Forest with different numbers of Trees

SL.randomForest_base <- function(Y, X, newX, family, mtry = ifelse(family$family ==
                                   "gaussian", max(floor(ncol(X) / 3), 1), floor(sqrt(ncol(X)))),
                                 ntree = ntree, nodesize = ifelse(family$family == "gaussian",
                                   5, 1
                                 ), maxnodes = NULL, importance = FALSE, ...) {
  cat(paste0(" - randomForest was started (ntree =", ntree, ") and is making predictions - "))
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
  cat(paste0("- randomForest was finished lasting: ", round(unclass(randomForest_time)["elapsed"], 2), " - "))
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
# search for algorithms in the enviroment
all_objects <- ls()
my_randomForest_objects_Pred <- grep("SL.randomForest_grid", all_objects)
all_randomForests_Pred <- all_objects[my_randomForest_objects_Pred]
rm(tuneGrid)
rm(all_objects)


# 3) Non-Linear effects via 'mgcv'

SL.mgcv <- function(Y, X, newX, family, obsWeights, ...) {

  # X <- X[,which(apply(X, 2, var)>0)] # columns need variation
  print(X)
  print(Y)
  cat(paste0("Correlation X:"))
  print(cor(X))
  cat(paste0("NEW X:"))
  print(newX)

  cat(paste0(" - MGCV started and is making predictions - "))
  mgcv_time <- system.time({
    SuperLearner:::.SL.require("mgcv")
    s <- mgcv:::s # s() is also used by 'gam' package - avoid clash

    # chose family dependent upon response variable
    Y <- as.vector(as.matrix(Y))
    if (all(Y == 0 | Y == 1)) {
      family$family <- binomial()
    } else {
      family$family <- gaussian() # scat() also reasonable but computationally intense
    }

    ## build model formula
    X <- data.frame(X)
    Y <- data.frame(Y)
    cat(paste0("X:"))
    print(X)
    cat(paste0("Y:"))
    print(Y)

    # considered categorical if less than 11 levels obseverd
    categorical_columns <- which(lengths(apply(X, 2, unique), use.names = TRUE) <= 10)
    # categorical variables are not used for smoothing (fished out by index)
    non_categ_index <- which((!colnames(X) %in% names(categorical_columns)))
    categ_index <- which((colnames(X) %in% names(categorical_columns)))

    # build model formula with thin plate smooth main terms (default)
    if (length(categ_index) > 0) {
      main_term <- paste0(paste0("s(", colnames(X)[non_categ_index], ",bs = \"tp\")", collapse = " + "),
        sep = " + ", paste0(colnames(X)[categ_index], collapse = " + ")
      )
    } else {
      main_term <- paste0("s(", colnames(X)[non_categ_index], ",bs = \"tp\")", collapse = " + ")
    }
    gam_formula <- as.formula(paste(names(Y), " ~ ", main_term))

    data <- cbind(Y, X)
    fit.gam <- try(mgcv::bam(gam_formula, family = family$family, data = data), silent = TRUE)

    # if fit was unsuccessful:
    # try and error with different number of knots and cubic regression splines
    # -> sparse matrix results
    num_knots <- c(9, 6, 3) # knots under consideration
    successful_gam <- FALSE # exit loop if TRUE
    if (any(class(fit.gam) == "try-error")) {
      index_i <- 1
      while ((successful_gam == FALSE) & (!is.na(num_knots[index_i]))) {
        cat(paste("Thin plate failed, try cubic regression splines with", num_knots[index_i], "knots."))
        if (length(categ_index) > 0) {
          main_term <- paste0(paste0("s(", colnames(X)[non_categ_index], ",bs = \"cs\", k = ", as.character(num_knots[index_i]), ")", collapse = " + "),
            sep = " + ", paste0(colnames(X)[categ_index], collapse = " + ")
          )
        } else {
          main_term <- paste0("s(", colnames(X)[non_categ_index], ",bs = \"cs\", k = ", num_knots[index_i], ")", collapse = " + ")
        }
        gam.formula <- as.formula(paste(names(Y), " ~ ", main_term))
        fit.gam <- try(mgcv::bam(gam.formula, data = data, family = family$family, drop.intercept = FALSE), silent = TRUE)

        if (any(class(fit.gam) == "try-error")) {
          successful_gam <- FALSE
        } else {
          successful_gam <- TRUE
        }

        # if everything fails -> Try GCV.Cp instead of fast.REML
        if ((successful_gam == FALSE) & min(num_knots) == num_knots[index_i]) {
          if (length(categ_index) > 0) {
            main_term <- paste0(paste0("s(", colnames(X)[non_categ_index], ",bs = \"cs\", k = 3)", collapse = " + "),
              sep = " + ", paste0(colnames(X)[categ_index], collapse = " + ")
            )
          } else {
            main_term <- paste0("s(", colnames(X)[non_categ_index], ",bs = \"cs\", k = 3)", collapse = " + ")
          }
          gam.formula <- as.formula(paste(names(Y), " ~ ", main_term))
          fit.gam <- try(mgcv::bam(gam.formula, data = data, family = family$family, method = "GCV.Cp", drop.intercept = FALSE), silent = TRUE)
          successful_gam <- TRUE
          cat(paste("Last try with estimation via generalized CV and 3 knots instead of fast REML."))
        }

        index_i <- index_i + 1
      }
    }

    cat(paste0("FORMULA:"))
    print(gam.formula)

    # predict unseen Y with seen X
    pred <- try(mgcv::predict.bam(fit.gam, newX, type = "response"), silent = TRUE)
    cat(paste0("PRED:"))
    print(pred)
    # cat(paste0("PRED1"))
    # print(pred)
    if (any(class(pred) == "try-error")) {
      cat(paste("MGCV Unsuccessful!"))
      pred <- rep(mean(Y), length(Y))
    }
    cat(paste0("PRED2"))
    print(pred)
    fit <- list(object = fit.gam)
    cat(paste0("FIT"))
    print(fit)
    out <- list(pred = pred, fit = fit)
    class(out$fit) <- c("SL.glm")
  })
  cat(paste0("- MGCV finished lasting: ", round(unclass(mgcv_time)["elapsed"], 2), " - "))
  return(out)
}

# 4) Boosting with tuning grid

SL.gbm_base <- function(Y, X, newX, family, obsWeights, gbm.trees = gbm.trees,
                        interaction.depth = interaction.depth, shrinkage = 0.001, ...) {
  cat(paste0(" - GBM started and is making predictions, interaction depth = ", interaction.depth, " and gbm.trees = ", gbm.trees, " - "))
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
    cat(paste("GBM Unsuccessful!"))
    pred <- rep(median(Y), length(Y))
  }
  fit <- list(object = fit.gbm, n.trees = best.iter)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.gbm")
  cat(paste0(" - GBM finished - "))
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
# search for algorithms in the enviroment
all_objects <- ls()
my_gbm_objects_Pred <- grep("SL.gbm_grid", all_objects)
all_gbm_Pred <- all_objects[my_gbm_objects_Pred]
rm(tuneGrid)
rm(all_objects)

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

  cat(paste0(" - GAM2 was started and is making predictions - "))
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
  cat(paste0("- GAM2 was finished lasting: ", round(unclass(gam_time)["elapsed"], 2), " - "))
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

SL.glm2 <- function(Y, X, newX, family, obsWeights, model = TRUE, ...) {
  Y <- as.vector(as.matrix(Y))
  if (all(Y == 0 | Y == 1)) {
    family <- binomial(link = "logit")
  } else {
    family <- gaussian(link = "identity")
  }
  print(family)
  if (is.matrix(X)) {
    X <- as.data.frame(X)
  }
  fit.glm <- glm(Y ~ .,
    data = X, family = family, weights = obsWeights,
    model = model
  )
  if (is.matrix(newX)) {
    newX <- as.data.frame(newX)
  }
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm"
  out <- list(pred = pred, fit = fit)
  return(out)
}