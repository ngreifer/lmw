get_w_from_X <- function(X, base.weights = NULL, s.weights = NULL) {

  #Remove linearly dependent columns
  qr_X <- qr(X)
  X <- X[, qr_X$pivot[seq_len(qr_X$rank)], drop = FALSE]

  #Treatment always in second column
  if (is.null(base.weights) && is.null(s.weights)) {
    weights <- X %*% chol2inv(chol(crossprod(X)))[,2]
  }
  else {
    if (is.null(s.weights)) s.weights <- rep(1, nrow(X))
    if (is.null(base.weights)) base.weights <- rep(1, nrow(X))
    wX <- s.weights * base.weights * X
    weights <- wX %*% chol2inv(chol(crossprod(X, wX)))[,2]
  }

  weights[X[,2] == 0] <- -weights[X[,2] == 0]

  return(weights)
}

get_X_from_formula_URI <- function(formula, data, treat, estimand, target, s.weights) {
  formula <- terms(formula, data = data)

  tt.factors <- attr(formula, "factors")

  if (is.null(treat)) {
    if (!any(colSums(tt.factors != 0) == 1)) {
      stop("Please supply an argument to 'treat' to identify the treatment variable.", call. = FALSE)
    }
    treat <- colnames(tt.factors)[which(colSums(tt.factors != 0) == 1)[1]]
    message(paste0("Using \"", treat, "\" as the treatment variable. If this is incorrect, please supply an argument to 'treat' to identify the treatment variable."))
  }

  #Remove treat main effect from tt.factors
  treat_main_effect <- any(colnames(tt.factors) == treat)
  if (!treat_main_effect) {
    stop("The model formula must include a main effect of treatment when type = \"URI\".", call. = FALSE)
  }
  tt.factors <- tt.factors[, colnames(tt.factors) != treat, drop = FALSE]

  #Extract terms that interact w/ treat
  new.f.terms <- colnames(tt.factors)
  interacts_with_treat <- tt.factors[treat,] > 0

  #Remove treat from interactions
  for (i in seq_along(new.f.terms)[interacts_with_treat]) {
    new.f.terms[i] <- paste(rownames(tt.factors)[rownames(tt.factors) != treat & tt.factors[, i] != 0], collapse = ":")
  }

  #Reconstruct formula and dataset without treat
  formula_without_treat <- terms(reformulate(new.f.terms, intercept = TRUE))

  mf <- model.frame(formula_without_treat, data = data, na.action = "na.pass")

  mf <- process_mf(mf)

  covs <- model.matrix(formula_without_treat, data = mf)

  assign <- attr(formula_without_treat, "term.labels")[attr(covs, "assign")[-1]]
  covs <- covs[,-1, drop = FALSE]

  #Extract treatment variable
  t <- model.response(model.frame(reformulate("0", treat), data = data, na.action = "na.pass"))

  t <- process_treat(t)

  #Center covariates at mean based on estimand; only affects weights if there
  #are interactions w/ treatment
  covs <- scale_covs(covs, t, estimand, target, s.weights)

  #Reconstruct X from centered covs by multiplying covs that interact with treat
  #by t
  X <- do.call("cbind", lapply(seq_along(new.f.terms), function(i) {
    if (interacts_with_treat[i]) {
      out <- t * covs[,assign == new.f.terms[i], drop = FALSE]
      colnames(out) <- paste(treat, colnames(covs)[assign == new.f.terms[i]], sep = ":")
    }
    else {
      out <- covs[,assign == new.f.terms[i], drop = FALSE]
    }
    return(out)
  }))

  #Add treatment and intercept to X
  X <- cbind(1, t, X)
  colnames(X)[1:2] <- c("(Intercept)", treat)

  attr(mf, "terms") <- NULL
  attr(t, "treat_name") <- treat

  return(list(X = X, t = t, mf = mf))
}

get_X_from_formula_MRI <- function(formula, data, treat, estimand, target, s.weights) {
  formula <- terms(formula, data = data)

  tt.factors <- attr(formula, "factors")

  if (is.null(treat)) {
    if (!any(colSums(tt.factors != 0) == 1)) {
      stop("Please supply an argument to 'treat' to identify the treatment variable.", call. = FALSE)
    }
    treat <- colnames(tt.factors)[which(colSums(tt.factors != 0) == 1)[1]]
    message(paste0("Using \"", treat, "\" as the treatment variable. If this is incorrect, please supply an argument to 'treat' to identify the treatment variable."))
  }
  tt.factors <- tt.factors[, colnames(tt.factors) != treat, drop = FALSE]

  #Extract terms that interact w/ treat
  new.f.terms <- colnames(tt.factors)
  interacts_with_treat <- tt.factors[treat,] > 0

  #Remove treat from interactions
  for (i in seq_along(new.f.terms)[interacts_with_treat]) {
    new.f.terms[i] <- paste(rownames(tt.factors)[rownames(tt.factors) != treat & tt.factors[, i] != 0], collapse = ":")
  }

  #Reconstruct formula and dataset without treat
  formula_without_treat <- terms(reformulate(new.f.terms, intercept = TRUE))

  mf <- model.frame(formula_without_treat, data = data, na.action = "na.pass")

  mf <- process_mf(mf)

  covs <- model.matrix(formula_without_treat, data = mf)

  assign <- attr(formula_without_treat, "term.labels")[attr(covs, "assign")[-1]]
  covs <- covs[,-1, drop = FALSE]

  #Extract treatment variable
  t <- model.response(model.frame(reformulate("0", treat), data = data, na.action = "na.pass"))

  t <- process_treat(t)

  #Center covariates at mean based on estimand
  covs <- scale_covs(covs, t, estimand, target, s.weights)

  #Reconstruct X from centered covs by multiplying all covs by t and adding
  #treat and intercept
  X <- cbind(1, t, covs, t * covs)
  colnames(X) <- c("(Intercept)", treat,
                   colnames(covs),
                   paste(treat, colnames(covs), sep = ":"))

  attr(mf, "terms") <- NULL
  attr(t, "treat_name") <- treat

  return(list(X = X, t = t, mf = mf))
}

scale_covs <- function(covs, treat, estimand, target, s.weights) {
  if (estimand == "ATE") {
    scaled_covs <- sweep(covs, 2L, colMeans_w(covs, s.weights), check.margin = FALSE)
  }
  else if (estimand == "ATT") {
    scaled_covs <- sweep(covs, 2L, colMeans_w(covs, s.weights, subset = treat == 1), check.margin = FALSE)
  }
  else if (estimand == "ATC") {
    scaled_covs <- sweep(covs, 2L, colMeans_w(covs, s.weights, subset = treat == 0), check.margin = FALSE)
  }
  else if (estimand == "CATE") {
    # scaled_covs <- sweep(covs, 2L, target, check.margin = FALSE)
  }

  scaled_covs
}
