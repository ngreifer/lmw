lmw_est.lmw_iv <- function(x, outcome, data = NULL, robust = TRUE, cluster = NULL, ...) {

  call <- match.call()

  data <- get_data(data, x)

  obj1 <- get_1st_stage_X_from_formula_iv(x$formula, data = data, treat = x$treat,
                                          iv = x$iv, method = x$method, estimand = x$estimand, target = x$target,
                                          s.weights = x$s.weights, target.weights = attr(x$target, "target.weights"),
                                          focal = x$focal)

  X <- obj1$X
  t <- as.numeric(x$treat == levels(x$treat)[2])

  outcome <- do.call("get_outcome", list(substitute(outcome), data, x$formula,
                                         obj1$X))
  outcome_name <- attr(outcome, "outcome_name")

  attributes(outcome) <- NULL

  ## model dimensions
  p <- ncol(obj1$X)

  weighted <- !is.null(x$s.weights) || !is.null(x$base.weights)

  if (weighted) {
    if (is.null(x$s.weights)) x$s.weights <- rep(1, nrow(obj1$X))
    if (is.null(x$base.weights)) x$base.weights <- rep(1, nrow(obj1$X))
    w <- x$s.weights * x$base.weights
    pos_w <- which(w > 0)
  }
  else {
    w <- rep(1, nrow(obj1$X))
    pos_w <- seq_along(w)
  }

  if (!is.null(x$fixef)) {
    for (i in seq_len(ncol(obj1$X))[-1]) {
      obj1$X[,i] <- demean(obj1$X[,i], x$fixef, w)
    }
    for (i in seq_len(ncol(obj1$A))) {
      obj1$A[,i] <- demean(obj1$A[,i], x$fixef, w)
    }
  }

  auxreg <- lm.wfit(x = obj1$X,
                    y = obj1$A,
                    w = w)

  t_fitted <- array(0, dim = dim(obj1$A), dimnames = dimnames(obj1$A))
  t_fitted[] <- auxreg$fitted.values

  obj2 <- get_2nd_stage_X_from_formula_iv(x$formula, data = data, treat = x$treat, treat_fitted = t_fitted,
                             method = x$method, estimand = x$estimand, target = x$target,
                             s.weights = x$s.weights, target.weights = attr(x$target, "target.weights"),
                             focal = x$focal)
  XZ <- obj2$X
  XZ[,colnames(t_fitted)] <- t_fitted

  if (!is.null(x$fixef)) {
    for (i in seq_len(ncol(XZ))[-1]) {
      XZ[,i] <- demean(XZ[,i], x$fixef, w)
    }
    outcome <- demean(outcome, x$fixef, w)
  }

  ## main regression
  fit <- lm.wfit(x = XZ, y = outcome, w = w)

  class(fit) <- c("lmw_est_iv", "lmw_est")

  fit$model.matrix <- XZ

  ok <- which(!is.na(fit$coefficients))

  #Replace t_fitted with A to compute residuals
  XZ[,colnames(obj1$A)] <- obj1$A
  fit$fitted.values <- drop(XZ[, ok, drop = FALSE] %*% fit$coefficients[ok])
  res <- outcome - fit$fitted.values
  fit$residuals <- res

  if (!is.null(x$fixef)) {
    fit$df.residual <- fit$df.residual - length(unique(x$fixef[pos_w])) + 1
    fit$fixef <- x$fixef
  }

  rss <- sum(fit$weights * res^2)
  fit$sigma <- sqrt(rss/fit$df.residual)

  #Subset model outputs to those with positive weights
  #for compatibility with vcovHC
  fit_sub <- subset_fit(fit)

  if (isTRUE(robust)) {
    if (is.null(cluster)) robust <- "HC3"
    else robust <- "HC1"
  }
  else if (isFALSE(robust)) {
    if (is.null(cluster)) robust <- "const"
    else {
      robust <- "HC1"
      warning("Setting robust = \"HC1\" because 'cluster' is non-NULL.", call. = FALSE)
    }
  }
  else if (!is.character(robust) || length(robust) != 1 ||
           !robust %in% eval(formals(sandwich::vcovHC.default)$type)) {
    stop("'robust' must be TRUE, FALSE, or one of the allowable inputs to the 'type' argument of sandwich::vcovHC().", call. = FALSE)
  }

  if (is.null(cluster)) {
    fit$vcov <- sandwich::vcovHC(fit, type = robust, ...)
  }
  else {
    if (inherits(cluster, "formula")) {
      cluster <- model.frame(cluster,
                             data = data[pos_w,, drop = FALSE],
                             na.action = na.pass)
    }
    else {
      cluster <- as.data.frame(cluster)
    }

    if (nrow(cluster) == nrow(data)) cluster <- cluster[pos_w,, drop = FALSE]
    else if (nrow(cluster) != length(pos_w)) {
      stop("'cluster' must have the same number of rows as the original data set.", call. = FALSE)
    }

    withCallingHandlers({
      fit$vcov <- sandwich::vcovCL(fit, type = robust, cluster = cluster, ...)
    },
    warning = function(w) {
      if (conditionMessage(w) != "clustered HC2/HC3 are only applicable to (generalized) linear regression models") warning(w)
      invokeRestart("muffleWarning")
    })
  }

  if (!is.null(x$fixef) && robust %in% c("const", "HC1")) {
    n <- length(pos_w)
    fit$vcov <- fit$vcov * (n - ncol(fit$model.matrix))/fit$df.residual
  }

  fit$lmw.weights <- x$weights
  fit$call <- call
  fit$estimand <- x$estimand
  fit$focal <- x$focal
  fit$method <- x$method
  fit$robust <- robust
  fit$outcome <- outcome_name
  fit$treat_levels <- levels(x$treat)

  fit
}

bread.lmw_est_iv <- function(x) {
  p <- x$rank
  p1 <- seq_len(p)
  Qr <- x$qr
  cov.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])

  n <- {
    if (is.null(x$weights)) nrow(x$model.matrix)
    else sum(x$weights > 0)
  }

  b <- cov.unscaled * n
  dimnames(b) <- list(names(x$coefficients)[p1], names(x$coefficients)[p1])
  return(b)
}
