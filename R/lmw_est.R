lmw_est <- function(x, ...) {
  UseMethod("lmw_est")
}

lmw_est.lmw <- function(x, outcome, data = NULL, robust = TRUE, cluster = NULL, ...) {

  call <- match.call()

  if (!inherits(x, "lmw")) {
    stop("'x' must be an lmw object.", call. = FALSE)
  }

  data <- get_data(data, x)

  outcome <- do.call("get_outcome", list(substitute(outcome), data, x$formula))
  outcome_name <- attr(outcome, "outcome_name")

  attributes(outcome) <- NULL

  #Get model matrix
  obj <- get_X_from_formula(x$formula, data = data, treat = x$treat,
                            method = x$method, estimand = x$estimand, target = x$target,
                            s.weights = x$s.weights, target.weights = attr(x$target, "target.weights"),
                            focal = x$focal)

  #Fit regression model; note use lm.[w]fit() instead of lm() because
  #we already have the model matrix (obj$X)
  w <- NULL
  rn <- rownames(obj$X)
  if (is.null(x$s.weights) && is.null(x$base.weights)) {
    if (!is.null(x$fixef)) {
      for (i in seq_len(ncol(obj$X))[-1]) {
        obj$X[,i] <- demean(obj$X[,i], x$fixef)
      }
      outcome <- demean(outcome, x$fixef)
    }
    fit <- lm.fit(x = obj$X, y = outcome)
    pos_w <- seq_along(outcome)
  }
  else {
    if (is.null(x$s.weights)) x$s.weights <- rep(1, length(outcome))
    if (is.null(x$base.weights)) x$base.weights <- rep(1, length(outcome))
    w <- x$s.weights * x$base.weights
    pos_w <- which(w > 0)
    obj$X <- obj$X[pos_w,, drop = FALSE]
    outcome <- outcome[pos_w]

    if (!is.null(x$fixef)) {
      for (i in seq_len(ncol(obj$X))[-1]) {
        obj$X[,i] <- demean(obj$X[,i], x$fixef[pos_w], w[pos_w])
      }
      outcome <- demean(outcome, x$fixef[pos_w], w[pos_w])
    }

    fit <- lm.wfit(x = obj$X, y = outcome, w = w[pos_w])

    non_pos_w <- which(w <= 0)
    fit$na.action <- setNames(non_pos_w, rn[non_pos_w])
    class(fit$na.action) <- "omit"
  }

  class(fit) <- c("lmw_est")

  fit$model.matrix <- obj$X

  n <- length(residuals(fit))

  if (!is.null(x$fixef)) {
    fit$df.residual <- fit$df.residual - nlevels(x$fixef) + 1
    fit$fixef <- x$fixef
  }

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

  if (robust == "const") { #Regular OLS vcov; faster than sandwich
    fit$vcov <- {
      if (is.null(w)) bread(fit)/length(residuals(fit)) * sum(residuals(fit)^2)/fit$df.residual
      else bread(fit)/length(residuals(fit)) * sum(w[pos_w] * residuals(fit)^2)/fit$df.residual
    }
  }
  else if (is.null(cluster)) {
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

  #Need to correct df for HC1 when fixed effects are present
  if (!is.null(x$fixef) && robust %in% c("const", "HC1")) {
    fit$vcov <- fit$vcov * (n - ncol(fit$model.matrix))/fit$df.residual
  }

  fit$call <- call
  fit$estimand <- x$estimand
  fit$focal <- x$focal
  fit$method <- x$method
  fit$robust <- robust
  fit$outcome <- outcome_name
  fit$treat_levels <- levels(x$treat)


  fit
}

print.lmw_est <- function(x, ...) {
  cat(sprintf("An %s object\n", class(x)[1]))
  cat(" - outcome:", x$outcome, "\n")
  cat(" - standard errors:", if (hasName(x$call, "cluster") && !is.null(x$call[["cluster"]])) "cluster",
      if (x$robust == "const") "usual" else sprintf("robust (%s)", x$robust), "\n")
  if (!inherits(x, "lmw_est_iv")) cat(" - estimand:", x$estimand, "\n")
  cat(" - method:", x$method, "\n")
  if (!is.null(x$fixef)) cat(" - fixed effects:", attr(x$fixef, "fixef_name"), "\n")
  cat("\n")
  cat("Use summary() to examine estimates, standard errors, p-values, and confidence intervals.", "\n")
  invisible(x)
}

lmw_est.formula <- function(x, data = NULL, estimand = "ATE", method = "URI", treat = NULL, base.weights = NULL,
                            s.weights = NULL, dr.method = "WLS", obj = NULL, target = NULL, target.weights = NULL,
                            contrast = NULL, focal = NULL,
                            outcome, robust = TRUE, cluster = NULL, ...) {
  call <- match.call()

  formula <- x

  method <- match_arg(method, c("URI", "MRI"))

  estimand <- process_estimand(estimand, target, obj)

  data <- process_data(data, obj)

  base.weights <- process_base.weights(base.weights, obj)

  s.weights <- process_s.weights(s.weights, obj)

  dr.method <- process_dr.method(dr.method, base.weights, method)

  treat_name <- process_treat_name(treat, formula, data, method, obj)

  #treat changes meaning from treatment name to treatment vector
  treat <- process_treat(treat_name, data)

  contrast <- process_contrast(contrast, treat, method)

  #treat_contrast has levels re-ordered so contrast is first
  treat_contrast <- apply_contrast_to_treat(treat, contrast)

  focal <- process_focal(focal, treat_contrast, estimand)

  # data <- get_data(data, x)

  outcome <- do.call("get_outcome", list(substitute(outcome), data, formula))
  outcome_name <- attr(outcome, "outcome_name")

  attributes(outcome) <- NULL

  #Get model matrix; note: use original treat, unlike
  obj <- get_X_from_formula(formula, data, treat, method, estimand,
                            target, s.weights, target.weights, focal)

  #Fit regression model; note use lm.[w]fit() instead of lm() because
  #we already have the model matrix (obj$X)
  w <- NULL
  rn <- rownames(obj$X)
  if (is.null(s.weights) && is.null(base.weights)) {
    fit <- lm.fit(x = obj$X, y = outcome)
    pos_w <- seq_along(outcome)
  }
  else {
    if (is.null(s.weights)) s.weights <- rep(1, length(outcome))
    if (is.null(base.weights)) base.weights <- rep(1, length(outcome))
    w <- s.weights * base.weights
    pos_w <- which(w > 0)
    obj$X <- obj$X[pos_w,, drop = FALSE]

    fit <- lm.wfit(x = obj$X, y = outcome[pos_w], w = w[pos_w])

    non_pos_w <- which(w <= 0)
    fit$na.action <- setNames(non_pos_w, rn[non_pos_w])
    class(fit$na.action) <- "omit"
  }

  class(fit) <- "lmw_est"

  fit$model.matrix <- obj$X

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
  fit$estimand <- estimand
  fit$focal <- focal
  fit$method <- method
  fit$robust <- robust
  fit$outcome <- outcome_name
  fit$treat_levels <- levels(treat)

  fit
}
