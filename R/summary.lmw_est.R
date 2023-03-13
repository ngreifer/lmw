#' Extract effect estimates and standard errors from `lmw_est` fits
#'
#' @description
#' `summary()` computes the treatment effect and potential outcome mean
#' estimates from the supplied `lmw_est` object. It functions similarly to
#' [summary.lm()] in producing estimate tables with the estimates,
#' standard errors, t-statistics, and p-values. Other model statistics can be
#' additionally requested.
#'
#' @details
#' `summary.lmw_est()` produces a table of treatment effect estimates
#' corresponding to all possible pairwise contrasts between the treatment
#' levels. These treatment effects generalize to the population implied by the
#' regression weights, which depends on the supplied estimand, whether sampling
#' weights were provided, and which of the MRI or URI models was requested. The
#' treatment effects are computed using linear contrasts of the outcome model
#' coefficients.
#'
#' When `method = "MRI"`, the potential outcome mean estimates are also
#' reported. These correspond to the potential outcome means in the population
#' implied by the regression weights. When `method = "URI"`, only the
#' treatment effects are estimated; the model-implied outcome means do not
#' correspond to the potential outcome means for the population implied by the
#' regression weights. That is, while the treatment effect generalizes to the
#' population defined by the regression weights, the estimated potential
#' outcome means do not and so are not reported.
#'
#' When `model = TRUE`, the model coefficients and their tests statistics
#' are additionally produced. It is inappropriate to interpret or report these
#' values as they have no causal interpretation. This is especially true when
#' using AIPW, as the model coefficients do not incorporate the augmentation
#' terms.
#'
#' @param object an `lmw_est` object; the output of a call to
#' `lmw_est`.
#' @param model `logical`; whether to produce a coefficient table for the
#' outcome model coefficients. Note that these values should not be interpreted
#' or reported so they are not produced by default.
#' @param ci `logical`; whether to include confidence intervals in the
#' output.
#' @param alpha when `ci = TRUE`, the alpha value used to compute the
#' critical test statistic for the confidence interval; equivalently, 1 minus
#' the confidence level (e.g., for a 99\% confidence interval, `alpha =
#' .01` should be specified). Default is .05 for a 95\% confidence interval.
#' @param \dots ignored.
#'
#' @return
#' A `summary.lmw_est` object with the following components:
#' \item{call}{the original call to `lmw_est()`}
#' \item{means}{a matrix
#' containting the estimated potential outcome means, their standard errors,
#' confidence interval limits (if requested with `ci = TRUE`),
#' t-statistics, and p-values. Omitted when `method = "URI"` or
#' `fixef` is not `NULL` and for `lmw_iv` objects.}
#' \item{coefficients}{a matrix containing the treatment effect estimates and
#' their standard errors, t-statistics, and p-values.When `ci = TRUE`, the
#' confidence limits `"95\%" CI L` (lower) and `"95\%" CI U` (upper)
#' will be included between the standard error and t-statistic columns. When
#' AIPW is used, z-statistics and z-tests are reported instead.}
#' \item{model.coefficients}{when `model = TRUE`, the coefficient table of
#' the model coefficients, which has the same columns as `coefficients.`}
#' \item{aliased}{when `model = TRUE`, a named logical vector showing if
#' the original coefficients are aliased (i.e., `NA`).}
#' \item{sigma, df,
#' r.squared, adj.r.squared}{the residual standard deviation, degrees of
#' freedom components, R-squared, and adjusted R-squared. See
#' [summary.lm()]. When AIPW is used, `sigma` and `df` are
#' omitted.}
#'
#' Other components containing information for printing are also included.
#'
#' @seealso
#' [lmw_est()] for fitting the outcome regression model,
#' [summary.lm()] for more information on the output components
#'
#' @examples
#' # See examples at `help("lmw_est")`

#' @exportS3Method summary lmw_est
summary.lmw_est <- function(object, model = FALSE, ci = TRUE, alpha = .05, ...) {
  treat_coef_inds <- seq_along(object$treat_levels)

  object$vcov <- .vcov.aliased(is.na(object$coefficients), object$vcov)

  coefs <- object$coefficients[treat_coef_inds]
  vcov <- object$vcov[treat_coef_inds, treat_coef_inds]
  rdf <- object$df.residual

  treat_name <- treat_name_from_coefs(names(coefs)[-1], object$treat_levels)

  coef_levels <- treat_levels_from_coefs(names(coefs)[-1], object$treat_levels,
                                         treat_name)

  treat_level_inds <- setNames(if (can_str2num(object$treat_levels)) str2num(object$treat_levels)
                               else seq_along(object$treat_levels),
                               object$treat_levels)

  model_coefs <- object$coefficients

  #Get means vector and vcov
  a <- diag(length(coef_levels))
  a[,1] <- 1

  means_order <- match(object$treat_levels, coef_levels)
  means <- drop(a %*% coefs)[means_order]
  means_vcov <- (a %*% vcov %*% t(a))[means_order, means_order]

  if (!is.null(object$focal)) {
    treat_levels <- c(object$focal, setdiff(object$treat_levels, object$focal))
    contrasts <- utils::combn(treat_levels, 2, simplify = FALSE)
  }
  else {
    contrasts <- utils::combn(object$treat_levels, 2, simplify = FALSE)
  }

  a0 <- setNames(rep(0, length(object$treat_levels)), object$treat_levels)
  a <- do.call("rbind", lapply(contrasts, function(i) {
    a0[i] <- switch(object$estimand, "ATT" = c(1, -1), c(-1, 1))
    a0
  }))

  effects <- drop(a %*% means)
  effects_vcov <- a %*% means_vcov %*% t(a)
  effects_se <- sqrt(diag(effects_vcov))
  effects_tval <- effects/effects_se
  effects_pval <- 2 * pt(abs(effects_tval), rdf, lower.tail = FALSE)

  effects_mat <- cbind(Estimate = effects,
                       `Std. Error` = effects_se,
                       `t value` = effects_tval,
                       `Pr(>|t|)` = effects_pval)

  if (ci) {
    effects_ci <- matrix(nrow = nrow(effects_mat), ncol = 2)
    colnames(effects_ci) <- paste0(round(100*(1-alpha), 1), "% CI " , c("L", "U"))
    t.crit <- abs(qt(alpha/2, rdf))
    effects_ci[,1] <- effects - effects_se*t.crit
    effects_ci[,2] <- effects + effects_se*t.crit

    effects_mat <- cbind(effects_mat[,1:2,drop=FALSE], effects_ci, effects_mat[,-(1:2),drop=FALSE])
  }

  rownames(effects_mat) <- lapply(contrasts, function(i) {
    if (object$estimand == "ATT") {
      sprintf("E[Y%s-Y%s]", treat_level_inds[i[1]], treat_level_inds[i[2]])
    }
    else {
      sprintf("E[Y%s-Y%s]", treat_level_inds[i[2]], treat_level_inds[i[1]])
    }
  })

  means_mat <- NULL

  if (object$method == "MRI" && is.null(object$fixef) && !inherits(object, "lmw_est_iv")) {

    means_se <- sqrt(diag(means_vcov))
    means_tval <- means/means_se
    means_pval <- 2 * pt(abs(means_tval), rdf, lower.tail = FALSE)

    means_mat <- cbind(Estimate = means,
                       `Std. Error` = means_se,
                       `t value` = means_tval,
                       `Pr(>|t|)` = means_pval)

    if (ci) {
      means_ci <- matrix(nrow = nrow(means_mat), ncol = 2)
      colnames(means_ci) <- paste0(round(100*(1-alpha), 1), "% CI " , c("L", "U"))
      t.crit <- abs(qt(alpha/2, rdf))
      means_ci[,1] <- means - means_se*t.crit
      means_ci[,2] <- means + means_se*t.crit

      means_mat <- cbind(means_mat[,1:2,drop=FALSE], means_ci, means_mat[,-(1:2),drop=FALSE])
    }

    rownames(means_mat) <- paste0("E[Y", treat_level_inds, "]")
  }

  f <- object$fitted.values
  r <- object$residuals
  y <- f + r
  w <- object$weights
  n <- length(f)
  m <- mean_w(f, w)
  if (is.null(w)) w <- rep(1, n)

  rss <- sum(w * r^2)
  tss <- sum(w * (y - mean_w(y, w))^2)

  r.squared <- 1 - rss/tss
  adj.r.squared <- 1 - (1 - r.squared) * ((n - 1)/rdf)

  sigma <- sqrt(rss/rdf)

  model_mat <- aliased <- NULL
  if (model) {
    model_se <- sqrt(diag(object$vcov))
    model_tval <- model_coefs/model_se

    model_mat <- cbind(Estimate = model_coefs,
                       `Std. Error` = model_se,
                       `t value` = model_tval,
                       `Pr(>|t|)` = 2 * pt(abs(model_tval), rdf, lower.tail = FALSE))
    aliased <- is.na(model_coefs)

    if (ci) {
      conf <- matrix(nrow = nrow(model_mat), ncol = 2)
      colnames(conf) <- paste0(round(100*(1-alpha), 1), "% CI " , c("L", "U"))

      conf[,1] <- model_coefs - model_se*t.crit
      conf[,2] <- model_coefs + model_se*t.crit

      model_mat <- cbind(model_mat[,1:2,drop=FALSE], conf, model_mat[,-(1:2),drop=FALSE])
    }
  }

  ans <- list(call = object$call,
              means = means_mat,
              coefficients = effects_mat,
              model.coefficients = model_mat,
              aliased = aliased,
              sigma = sigma,
              df = c(object$rank, rdf, NCOL(object$qr[["qr"]])),
              r.squared = r.squared,
              adj.r.squared = adj.r.squared,
              estimand = object$estimand,
              focal = object$focal,
              treat_levels = object$treat_levels,
              fixef_name = attr(object$fixef, "fixef_name"))

  class(ans) <- "summary.lmw_est"
  return(ans)
}

#' @exportS3Method print summary.lmw_est
print.summary.lmw_est <- function(x, digits = max(3, getOption("digits") - 3),
                                  signif.stars = getOption("show.signif.stars"), ...) {
  # cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
  #     "\n", sep = "")

  df <- x$df
  rdf <- df[2L]

  treat_level_inds <- setNames(if (can_str2num(x$treat_levels)) str2num(x$treat_levels)
                               else seq_along(x$treat_levels),
                               x$treat_levels)

  if (!is.null(x$coefficients)) {
    cat("\nEffect estimates:\n")
    coefs <- na.omit(x$coefficients)
    rownames(coefs) <- num2sub(rownames(coefs))
    rownames(coefs) <- formatC(if (!is.null(x$focal)) gsub("]", paste0("|A=", treat_level_inds[x$focal],"]"), rownames(coefs), fixed = TRUE)
                               else if (x$estimand == "CATE") gsub("]", "|X=x*]", rownames(coefs), fixed = TRUE)
                               else rownames(coefs))

    ci <- ncol(coefs) == 6
    printCoefmat.args <- list(digits = digits, signif.stars = signif.stars,
                              na.print = ".", cs.ind = if (ci) 1:4 else 1:2,
                              tst.ind = if (ci) 5 else 3, ...)
    do.call("printCoefmat", c(list(coefs), printCoefmat.args))

    if (!is.null(df)) {
      cat("\nResidual standard error:", format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom")
    }
    if (!is.null(x$fixef_name)) {
      cat("\nEstimated with fixed effects for", x$fixef_name)
    }
    cat("\n")
  }

  if (!is.null(x$means)) {
    cat("\nPotential outcome means:\n")
    means <- na.omit(x$means)
    rownames(means) <- num2sub(rownames(means))
    rownames(means) <- formatC(if (!is.null(x$focal)) gsub("]", paste0("|A=", treat_level_inds[x$focal],"]"), rownames(means), fixed = TRUE)
                               else if (x$estimand == "CATE") gsub("]", "|X=x*]", rownames(means), fixed = TRUE)
                               else rownames(means))

    ci <- ncol(means) == 6
    printCoefmat.args <- list(digits = digits, signif.stars = signif.stars,
                              na.print = ".", cs.ind = if (ci) 1:4 else 1:2,
                              tst.ind = if (ci) 5 else 3, ...)
    do.call("printCoefmat", c(list(means), printCoefmat.args))
  }

  if ((!is.null(x$coefficients) || !is.null(x$means)) && !can_str2num(x$treat_levels)) {
    cat("\nKey: ")
    cat(paste(paste( treat_level_inds, "=", names(treat_level_inds)), collapse = "; "))
    cat("\n")
  }

  if (!is.null(x$model.coefficients)) {
    if (length(x$aliased) == 0L) {
      cat("\nNo Coefficients\n")
    }
    else {
      if (nsingular <- sum(is.na(x$model.coefficients)))
        cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", sep = "")
      else cat("\nModel coefficients:\n")
      do.call("printCoefmat", c(list(x$model.coefficients), printCoefmat.args))
    }
    cat("\nAll covariates centered at",
        switch(x$estimand,
               "ATE" = "their sample means.",
               "ATT" = "their means in the treated group (A=1).",
               "ATC" = "their means in the control group (A=0).",
               "CATE" = "their target values (x*)."))
    cat("\nMultiple R-squared: ", formatC(x$r.squared, digits = digits),
        ",  Adjusted R-squared: ", formatC(x$adj.r.squared,
                                           digits = digits), sep = "")

  }
  invisible(x)
}

#' @exportS3Method model.matrix lmw_est
model.matrix.lmw_est <- function(object, ...) {
  mm <- object$model.matrix
  if (called_from("meatHC", "meatCL")) {
    mm <- mm[object$weights > 0, , drop = FALSE]
  }
  mm
}

#' @exportS3Method hatvalues lmw_est
hatvalues.lmw_est <- function(model, ...) {
  h <- hat_fast(model$model.matrix, model$weights, model$fixef)
  if (called_from("meatHC", "meatCL")) {
    h <- h[model$weights > 0]
  }
  h
}

#' @exportS3Method estfun lmw_est
estfun.lmw_est <- function (x, ...) {

  x <- subset_fit(x)

  xmat <- model.matrix(x)
  xmat <- naresid(x$na.action, xmat)
  if (any(alias <- is.na(coef(x))))
    xmat <- xmat[, !alias, drop = FALSE]
  wts <- weights(x)
  if (is.null(wts)) wts <- 1
  res <- residuals(x)
  rval <- as.vector(res) * wts * xmat
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL

  return(rval)
}

#' @exportS3Method vcov lmw_est
vcov.lmw_est <- function(object, complete = TRUE, ...) {
  object$vcov
}

#' @exportS3Method bread lmw_est
bread.lmw_est <- function(x) {
  p <- x$rank
  p1 <- seq_len(p)
  Qr <- x$qr
  cov.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])

  #Note: in presence of fixed effects, n is incorrect
  # df <- c(p, x$df.residual)
  # b <- cov.unscaled * as.vector(sum(df))

  b <- cov.unscaled * sum(x$weights > 0)
  dimnames(b) <- list(names(x$coefficients[p1]), names(x$coefficients[p1]))
  return(b)
}

#' @exportS3Method weights lmw_est
weights.lmw_est <- function(object, ...) {
  wts <- object$weights
  if (called_from("meatHC", "meatCL")) {
    wts <- wts[wts > 0]
  }
  wts
}

#' @exportS3Method predict lmw_est
predict.lmw_est <- function(object, newdata, ...) {
 if (!missing(newdata)) {
   warning("`predict()` cannot be used with `lmw_est` objects to generate predictions for new observations. Returning fitted values from the original dataset.")
 }
  object$fitted.values
}
