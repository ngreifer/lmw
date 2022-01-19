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

  contrasts <- combn(object$treat_levels, 2, simplify = FALSE)

  a0 <- setNames(rep(0, length(object$treat_levels)), object$treat_levels)
  a <- do.call("rbind", lapply(contrasts, function(i) {
    a0[i] <- c(-1, 1)
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
    sprintf("E[Y%s-Y%s]", treat_level_inds[i[2]], treat_level_inds[i[1]])
  })

  means_mat <- NULL

  if (object$type == "MRI") {

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
  w <- object$weights
  n <- length(r)
  if (is.null(w)) {
    mss <- sum((f - mean(f))^2)
    rss <- sum(r^2)
  }
  else {
    m <- sum(w * f/sum(w))
    mss <- sum(w * (f - m)^2)
    rss <- sum(w * r^2)
  }

  r.squared <- mss/(mss + rss)
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
              treat_levels = object$treat_levels)

  class(ans) <- "summary.lmw_est"
  return(ans)
}

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

    cat("\nResidual standard error:", format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom")
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
      if (nsingular <- df[3L] - df[1L])
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

model.matrix.lmw_est <- function(object, ...) {
  object$model.matrix
}
hatvalues.lmw_est <- function(model, ...) {
  hat_fast(model$model.matrix, model$weights)
}
estfun.lmw_est <- function (x, ...) {
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
vcov.lmw_est <- function(object, complete = TRUE, ...) {
  object$vcov

  # p <- object$rank
  # p1 <- seq_len(p)
  # Qr <- object$qr
  # cov.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  #
  # if (is.null(object$weights)) {
  #   rss <- sum(object$residuals^2)
  # }
  # else {
  #   rss <- sum(object$weights * object$residuals^2)
  # }
  # resvar <- rss/object$df.residual
  #
  # .vcov.aliased(is.na(object$coefficients),
  #               resvar * cov.unscaled,
  #               complete)
}
bread.lmw_est <- function(x) {
  p <- x$rank
  p1 <- seq_len(p)
  Qr <- x$qr
  cov.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  df <- c(p, x$df.residual)

  b <- cov.unscaled * as.vector(sum(df))
  dimnames(b) <- list(names(x$coefficients[p1]), names(x$coefficients[p1]))
  return(b)
}