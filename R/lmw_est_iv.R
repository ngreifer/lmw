lmw_est.lmw_iv <- function(x, outcome, data = NULL, robust = TRUE, cluster = NULL, ...) {

  call <- match.call()

  if (!inherits(x, "lmw_iv")) {
    stop("'x' must be an lmw_iv object.", call. = FALSE)
  }

  data <- get_data(data, x)

  outcome <- do.call("get_outcome", list(substitute(outcome), data, x$formula))
  outcome_name <- attr(outcome, "outcome_name")

  attributes(outcome) <- NULL

  obj1 <- get_1st_stage_X_from_formula_iv(x$formula, data = data, treat = x$treat,
                                          iv = x$iv, type = x$type, estimand = x$estimand, target = x$target,
                                          s.weights = x$s.weights, focal = x$focal)

  X <- obj1$X
  t <- as.numeric(x$treat == levels(x$treat)[2])

  ## model dimensions
  n <- NROW(outcome)
  p <- ncol(X)

  w <- NULL
  rn <- rownames(x$X)

  if (is.null(x$s.weights) && is.null(x$base.weights)) {
    pos_w <- seq_len(n)
  }
  else {
    if (is.null(x$s.weights)) x$s.weights <- rep(1, length(outcome))
    if (is.null(x$base.weights)) x$base.weights <- rep(1, length(outcome))
    w <- x$s.weights * x$base.weights
    pos_w <- which(w > 0)
  }

  offset <- rep(0, n)

  if (is.null(w)) {
    auxreg <- lm.fit(X, t)
  }
  else {
    auxreg <- lm.wfit(X[pos_w,, drop = FALSE], t[pos_w], w[pos_w])
  }

  t_fitted <- auxreg$fitted.values

  obj2 <- get_2nd_stage_X_from_formula_iv(x$formula, data = data, treat = x$treat, treat_fitted = t_fitted,
                                          type = x$type, estimand = x$estimand, target = x$target,
                                          s.weights = x$s.weights, focal = x$focal)

  XZ <- obj2$X

  ## main regression
  if (is.null(w)) {
    fit <- lm.fit(x = XZ, y = outcome)
  }
  else {
    fit <- lm.wfit(x = XZ[pos_w,, drop = FALSE], y = outcome[pos_w], w = w[pos_w])

    non_pos_w <- which(w <= 0)
    fit$na.action <- setNames(non_pos_w, rn[non_pos_w])
    class(fit$na.action) <- "omit"
  }

  class(fit) <- c("lmw_est_iv", "lmw_est")

  fit$model.matrix <- XZ

  ok <- which(!is.na(fit$coefficients))
  res <- numeric(n)

  #Replace t_fitted with t to compute residuals
  XZ[,2] <- t
  res[pos_w] <- outcome[pos_w] - drop(XZ[, ok, drop = FALSE] %*% fit$coefficients[ok])
  fit$residuals <- res

  rss <- if (is.null(w)) sum(res^2) else sum(w * res^2)
  fit$sigma <- sqrt(rss/fit$df.residual)

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
    fit$vcov <- sandwich::vcovCL(fit, type = robust, cluster = cluster, ...)
  }

  fit$call <- call
  fit$estimand <- x$estimand
  fit$focal <- x$focal
  fit$type <- x$type
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
  dimnames(b) <- list(names(x$coefficients), names(x$coefficients))
  return(b)
}
