lmw_est <- function(x, outcome, data = NULL) {

  call <- match.call()

  if (!inherits(x, "lmw")) {
    stop("'x' must be an lmw object.", call. = FALSE)
  }

  f_env <-environment(x$formula)

  if (is.null(data)) {
    data <- try(eval(x$call$data, envir = f_env), silent = TRUE)
    if (length(data) == 0 || inherits(data, "try-error") || length(dim(data)) != 2 || nrow(data) != length(x[["treat"]])) {
      data <- try(eval(x$call$data, envir = parent.frame()), silent = TRUE)
      if (length(data) == 0 || inherits(data, "try-error") || length(dim(data)) != 2 || nrow(data) != length(x[["treat"]])) {
        data <- x[["model"]][["data"]]
        if (length(data) == 0 || nrow(data) != length(x[["treat"]])) {
          stop("A valid dataset could not be found. Please supply an argument to 'data' containing the original dataset used to estimate the weights.", call. = FALSE)
        }
      }
    }
  }

  if (!is.data.frame(data)) {
    if (is.matrix(data)) data <- as.data.frame.matrix(data)
    else stop("'data' must be a data frame.", call. = FALSE)
  }
  if (nrow(data) != length(x$treat)) {
    stop("'data' must have as many rows as there were units in the original call to lmw().", call. = FALSE)
  }

  tt <- terms(x$formula, data = data)

  outcome_char <- deparse1(substitute(outcome))
  if (is.null(outcome)) {
    if (attr(tt, "response") == 0) {
      stop("'outcome' must be supplied.", call. = FALSE)
    }
    mf <- model.frame(tt, data = data)
    outcome <- model.response(mf)
  }
  else {
    if ((is.numeric(outcome) || is.logical(outcome)) && length(outcome) == nrow(data)) {

    }
    else if (is.character(outcome) && length(outcome) == 1) {
      outcome_char <- outcome
      outcome <- try(eval(str2expression(outcome), data))
      if (length(outcome) == 0 || inherits(outcome, "try-error")) {
        stop("The outcome variable must be present in the dataset.", call. = FALSE)
      }
    }
    else {
      stop("'outcome' must be the name of the outcome variable or a vector containing the outcome values for all units.", call. = FALSE)
    }
  }

  if (!is.numeric(outcome) && !is.logical(outcome)) {
    stop("The outcome variable must be numeric.", call. = FALSE)
  }

  obj <- switch(x$type,
                "URI" = get_X_from_formula_URI(x$formula, data = data, treat = attr(x$treat, "treat_name"),
                                               estimand = x$estimand, target = x$target, s.weights = x$s.weights),
                "MRI" = get_X_from_formula_MRI(x$formula, data = data, treat = attr(x$treat, "treat_name"),
                                               estimand = x$estimand, target = x$target, s.weights = x$s.weights)
  )

  if (is.null(x$s.weights)) {
    fit <- lm.fit(x = obj$X, y = outcome)
  }
  else {
    fit <- lm.fit(x = obj$X, y = outcome, w = x$s.weights)
  }

  fit$call <- call
  class(fit) <- "lmw_est"
  fit
}

summary.lmw_est <- function(object, robust = TRUE, cluster = NULL, ...) {
  coefs <- object$coefficients

  est <- coefs[2]
  EY0 <- coefs[1]
  EY1 <- EY0 + est

  if (isTRUE(robust)) {
    if (is.null(cluster)) robust <- "HC3"
    else robust <- "HC1"
  }
  else if (isFALSE(robust)) {
    cluster <- NULL
    robust <- "const"
  }
  else if (!is.character(robust) || length(robust) != 1 ||
           !robust %in% eval(formals(sandwich::vcovHC.default)$type)) {
    stop("'robust' must be TRUE, FALSE, or one of the allowable inputs to the 'type' argument of sandwich::vcovHC().", call. = FALSE)
  }
}

#Extract bread from lm.fit() output
bread.lm.fit <- function(x) {
  p <- x$rank
  p1 <- seq_len(p)
  Qr <- x$qr
  cov.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  df <- c(p, x$df.residual)

  return(cov.unscaled * as.vector(sum(df)))
}
meat.lm.fit <- function(x, mm, type = "HC3") {
  X <- mm
  if (any(alias <- is.na(x$coefficients)))
    X <- X[, !alias, drop = FALSE]
  attr(X, "assign") <- NULL
  n <- NROW(X)
  diaghat <- hat_fast(X, x$weights)
  df <- n - NCOL(X)
  ef <- estfun(x)
  res <- rowMeans(ef/X, na.rm = TRUE)
  all0 <- apply(abs(ef) < .Machine$double.eps, 1L, all)
  res[all0] <- 0
  if (any(all0) && substr(type, 1L, 1L) == "c") {
    if (inherits(x, "glm")) {
      res <- as.vector(residuals(x, "working")) * weights(x,
                                                          "working")
      if (!(substr(x$family$family, 1L, 17L) %in% c("poisson",
                                                    "binomial", "Negative Binomial"))) {
        res <- res * sum(weights(x, "working"), na.rm = TRUE)/sum(res^2,
                                                                  na.rm = TRUE)
      }
    }
    else if (inherits(x, "lm")) {
      res <- as.vector(residuals(x))
      if (!is.null(weights(x)))
        res <- res * weights(x)
    }
  }

  type <- match.arg(type)
  if (type == "HC")
    type <- "HC0"
  switch(type, const = {
    omega <- function(residuals, diaghat, df) rep(1,
                                                  length(residuals)) * sum(residuals^2)/df
  }, HC0 = {
    omega <- function(residuals, diaghat, df) residuals^2
  }, HC1 = {
    omega <- function(residuals, diaghat, df) residuals^2 *
      length(residuals)/df
  }, HC2 = {
    omega <- function(residuals, diaghat, df) residuals^2/(1 -
                                                             diaghat)
  }, HC3 = {
    omega <- function(residuals, diaghat, df) residuals^2/(1 -
                                                             diaghat)^2
  }, HC4 = {
    omega <- function(residuals, diaghat, df) {
      n <- length(residuals)
      p <- as.integer(round(sum(diaghat), digits = 0))
      delta <- pmin(4, n * diaghat/p)
      residuals^2/(1 - diaghat)^delta
    }
  }, HC4m = {
    omega <- function(residuals, diaghat, df) {
      gamma <- c(1, 1.5)
      n <- length(residuals)
      p <- as.integer(round(sum(diaghat), digits = 0))
      delta <- pmin(gamma[1], n * diaghat/p) + pmin(gamma[2],
                                                    n * diaghat/p)
      residuals^2/(1 - diaghat)^delta
    }
  }, HC5 = {
    omega <- function(residuals, diaghat, df) {
      k <- 0.7
      n <- length(residuals)
      p <- as.integer(round(sum(diaghat), digits = 0))
      delta <- pmin(n * diaghat/p, pmax(4, n * k *
                                          max(diaghat)/p))
      residuals^2/sqrt((1 - diaghat)^delta)
    }
  })

  omega <- omega(res, diaghat, df)
  rval <- sqrt(omega) * X
  rval <- crossprod(rval)/n
  return(rval)
}
