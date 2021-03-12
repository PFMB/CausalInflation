rm(list = ls())

## Tests the functionality of the screeners constructed in 32LearnerLibrary.R
## with regard to different input types (e.g. factor vs numeric) and ensures
## a standardized input and output of each screening algorithm. Just paste every
## screener after line #52 and run the script.

gen_data <- function(n = 1e2, p = 129, p_fac = 20, p_num = p - p_fac, p_rel = 5, mk_mat = FALSE) {
  
  ## generate some test data

  # p: total no. of preds
  # n: no. of observations
  # p_fac: no. of factor variables
  # p_num: no. of metric variables
  # p_rel: no. of X vars that drive Y
  # mk_mat: make matrix, check if screen/pred can handle matrix and df equally

  if (p_fac > 0 & mk_mat) stop("Not meaningful: Factors are lost in matrix!")

  # generate some cont and categ data where only a part of X (rel_col) drives Y
  X_fac <- matrix(sample(c(3, -5, 10), size = p_fac * n, replace = TRUE), nrow = n)
  X_num <- matrix(rnorm(n * p_num), ncol = p_num)
  X <- cbind(X_num, X_fac)
  X <- X[, sample(1:ncol(X), ncol(X))] # shuffle columns

  # randomly pick columns that drive Y
  rel_col <- sample(seq(ncol(X)), p_rel, replace = FALSE)
  if (p_rel == 0) {
    lin_pred <- rep(0, n)
  } else {
    # building the lin_pred like that might not be ideal for factor variables
    lin_pred <- X[, rel_col, drop = FALSE] %*% rep(1,p_rel)
  }

  # turn categor. into factors and make df
  X <- as.data.frame(X)
  X <- lapply(X, function(var) {
    if (length(unique(var)) < 5) {
      return(as.factor(var))
    }
    var
  })

  X <- as.data.frame(X)

  # add disturbance: variable that has variance 0
  X[, sample(1:ncol(X), 1)] <- rep(1, n)

  if (mk_mat) X <- as.matrix(X)
  list(preds = X, Y_num = lin_pred + rnorm(n), Y_bin = rbinom(n, 1, plogis(lin_pred + rnorm(n))), rel_col = rel_col)
}

#############

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


#############

res <- sapply(1:100, function(b) {
  # if your algo is not aimed at trying to select a predefined number, insert
  # p_rel = rpois(1,5) instead, or a constant respectively
  dat <- gen_data(p = 100, p_rel = rpois(1,5), p_fac = rpois(1,10))
  # relative selection frequencies
  Y <- dat$Y_num
  X <- dat$preds
  SL.glm.interaction_info(Y,X, newX = X)$pred
})
