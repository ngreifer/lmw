#' @exportS3Method lmw_est lmw_aipw
#' @rdname lmw_est
lmw_est.lmw_aipw <- function(x, outcome, data = NULL, robust = TRUE, cluster = NULL, ...) {

  call <- match.call()

  if (!inherits(x, "lmw")) {
    stop("'x' must be an lmw object.", call. = FALSE)
  }

  if ("robust" %in% names(call)[-1]) {
    warning("'robust' is ignored for lmw_aipw objects. See ?lmw_est.lmw_aipw for details.", call. = FALSE)
  }
  if ("cluster" %in% names(call)[-1]) {
    warning("'cluster' is ignored for lmw_aipw objects. See ?lmw_est.lmw_aipw for details.", call. = FALSE)
  }

  data <- get_data(data, x)

  #Get model matrix
  obj <- get_X_from_formula(x$formula, data = data, treat = x$treat,
                            method = x$method, estimand = x$estimand, target = x$target,
                            s.weights = x$s.weights, target.weights = attr(x$target, "target.weights"),
                            focal = x$focal)

  outcome <- do.call("get_outcome", list(substitute(outcome), data, x$formula,
                                         obj$X))
  outcome_name <- attr(outcome, "outcome_name")

  attributes(outcome) <- NULL

  #Fit regression model; note use lm.[w]fit() instead of lm() because
  #we already have the model matrix (obj$X)

  if (is.null(x$s.weights)) {
    if (!is.null(x$fixef)) {
      for (i in seq_len(ncol(obj$X))[-1]) {
        obj$X[,i] <- demean(obj$X[,i], x$fixef)
      }
      outcome <- demean(outcome, x$fixef)
    }
    fit <- lm.fit(x = obj$X, y = outcome)
    reg_w <- rep(1, length(outcome))
    pos_w <- seq_along(outcome)
  }
  else {
    reg_w <- x$s.weights
    pos_w <- which(reg_w > 0)

    if (!is.null(x$fixef)) {
      for (i in seq_len(ncol(obj$X))[-1]) {
        obj$X[,i] <- demean(obj$X[,i], x$fixef, reg_w)
      }
      outcome <- demean(outcome, x$fixef, reg_w)
    }

    fit <- lm.wfit(x = obj$X, y = outcome, w = reg_w)

    # non_pos_w <- which(reg_w <= 0)
    # fit$na.action <- setNames(non_pos_w, rn[non_pos_w])
    # class(fit$na.action) <- "omit"
  }

  fit$model.matrix <- obj$X

  if (!is.null(x$fixef)) {
    fit$df.residual <- fit$df.residual - length(unique(x$fixef[pos_w])) + 1
    fit$fixef <- x$fixef
  }

  #Get augmentation terms
  #E[Y1] = mu1 + ipw1 - aug1
  nA <- nlevels(x$treat)
  mu <- ipw <- aug <- numeric(nA)
  yA <- vector("list", nA)

  #"estimand" weights; ensure means are taken
  #over correct estimand.
  ew <- {
    if (is.null(x$focal)) as.numeric(reg_w > 0)
    else as.numeric(reg_w > 0)*as.numeric(x$treat == x$focal)
  }

  aipw_w <- reg_w * x$base.weights

  for (i in 1:nA) {
    obj_i <- get_X_from_formula(x$formula, data = data, treat = x$treat,
                                method = x$method, estimand = x$estimand, target = x$target,
                                s.weights = x$s.weights, target.weights = attr(x$target, "target.weights"),
                                focal = x$focal, treat_fixed = levels(x$treat)[i])

    yA[[i]] <-  drop(obj_i$X[pos_w,, drop = FALSE] %*% fit$coefficients)

    mu[i] <- mean_w(yA[[i]], reg_w[pos_w] * ew[pos_w])

    ipw[i] <- mean_w(outcome[pos_w], aipw_w[pos_w], x$treat[pos_w] == levels(x$treat)[i])

    aug[i] <- mean_w(fit$fitted.values[pos_w], aipw_w[pos_w], x$treat[pos_w] == levels(x$treat)[i])
  }

  #Get covariance of coefs and augmentation terms
  p <- ncol(obj$X)
  n <- length(pos_w)

  # theta <- tcrossprod(c(fit$coefficients, mu1, mu0, ipw1 - aug1, ipw0 - aug0))
  psi <- do.call("rbind", c(
    list(t(reg_w[pos_w] * fit$residuals[pos_w] * obj$X[pos_w,, drop = FALSE])),
    lapply(1:nA, function(i) (reg_w*ew)[pos_w]*(yA[[i]] - mu[i])), #replace ew with PS to get "correct" SEs for AIPW that agree with PSweight
    lapply(1:nA, function(i) (x$treat[pos_w] == levels(x$treat)[i]) * aipw_w[pos_w] * (fit$residuals[pos_w] - (ipw[i] - aug[i])))
  ))

  #rowMeans(psi) should be all 0s

  B <- tcrossprod(psi)/n

  A <- matrix(0, nrow = p + 2*nA, ncol = p + 2*nA)
  A[1:p, 1:p] <- crossprod(obj$X[pos_w,, drop = FALSE])/n

  for (i in 1:nA) {
    A[p + i, p + i] <- sum(ew)/n
    A[p + nA + i, p + nA + i] <- sum(aipw_w[x$treat == levels(x$treat)[i]])/n
    A[p + i, c(1, i)] <- -sum(ew)/n
    A[p + nA + i, 1:p] <- colSums((x$treat == levels(x$treat)[i]) * aipw_w * obj$X)/n
  }

  # A[p + 1, 1:p] <- -colSums(obj1$X)/n
  # A[p + 2, 1:p] <- -colSums(obj0$X)/n
  #Because covariates are centered, colSums/n = 0 except for intercept and treat,
  #which = 1

  Ainv <- solve(A)

  V <- Ainv %*% B %*% t(Ainv)/n

  #VCOV of outcome model coefficients
  vcov <- V[1:p, 1:p]
  dimnames(vcov) <- list(names(fit$coefficients), names(fit$coefficients))

  #VCOV of AIPW parameters (means and augmentation)
  treat_level_inds <- setNames(if (can_str2num(levels(x$treat))) str2num(levels(x$treat))
                               else seq_along(levels(x$treat)),
                               levels(x$treat))
  coef_aipw_names <- c(paste0("mu", treat_level_inds), paste0("aug", treat_level_inds))
  vcov_aipw <- V[-(1:p), -(1:p)]
  dimnames(vcov_aipw) <- list(coef_aipw_names, coef_aipw_names)
  coef_aipw <- setNames(unname(c(mu, ipw - aug)), coef_aipw_names)

  fit$vcov <- vcov
  fit$lmw.weights <- x$weights
  fit$call <- call
  fit$estimand <- x$estimand
  fit$focal <- x$focal
  fit$method <- x$method
  fit$robust <- "HC0"
  fit$outcome <- outcome_name
  fit$treat_levels <- levels(x$treat)
  fit$coef_aipw <- coef_aipw
  fit$vcov_aipw <- vcov_aipw

  class(fit) <- c("lmw_est_aipw", "lmw_est")
  fit
}

#' @exportS3Method summary lmw_est_aipw
#' @rdname summary.lmw_est
summary.lmw_est_aipw <- function(object, model = FALSE, ci = TRUE, alpha = .05, ...) {
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
  a <- cbind(diag(length(coef_levels)), diag(length(coef_levels)))

  means_order <- match(object$treat_levels, coef_levels)
  means <- drop(a %*% object$coef_aipw)[means_order]
  means_vcov <- (a %*% object$vcov_aipw %*% t(a))[means_order, means_order]

  contrasts <- utils::combn(object$treat_levels, 2, simplify = FALSE)

  a0 <- setNames(rep(0, length(object$treat_levels)), object$treat_levels)
  a <- do.call("rbind", lapply(contrasts, function(i) {
    a0[i] <- c(-1, 1)
    a0
  }))

  effects <- drop(a %*% means)
  effects_vcov <- a %*% means_vcov %*% t(a)
  effects_se <- sqrt(diag(effects_vcov))
  effects_tval <- effects/effects_se
  effects_pval <- 2 * pnorm(abs(effects_tval), lower.tail = FALSE)

  effects_mat <- cbind(Estimate = effects,
                       `Std. Error` = effects_se,
                       `z value` = effects_tval,
                       `Pr(>|z|)` = effects_pval)

  if (ci) {
    effects_ci <- matrix(nrow = nrow(effects_mat), ncol = 2)
    colnames(effects_ci) <- paste0(round(100*(1-alpha), 1), "% CI " , c("L", "U"))
    z.crit <- abs(qnorm(alpha/2))
    effects_ci[,1] <- effects - effects_se*z.crit
    effects_ci[,2] <- effects + effects_se*z.crit

    effects_mat <- cbind(effects_mat[,1:2,drop=FALSE], effects_ci, effects_mat[,-(1:2),drop=FALSE])
  }

  rownames(effects_mat) <- lapply(contrasts, function(i) {
    sprintf("E[Y%s-Y%s]", treat_level_inds[i[2]], treat_level_inds[i[1]])
  })

  means_mat <- NULL

  if (object$method == "MRI" && is.null(object$fixef)) {

    means_se <- sqrt(diag(means_vcov))
    means_tval <- means/means_se
    means_pval <- 2 * pnorm(abs(means_tval), lower.tail = FALSE)

    means_mat <- cbind(Estimate = means,
                       `Std. Error` = means_se,
                       `z value` = means_tval,
                       `Pr(>|z|)` = means_pval)

    if (ci) {
      means_ci <- matrix(nrow = nrow(means_mat), ncol = 2)
      colnames(means_ci) <- paste0(round(100*(1-alpha), 1), "% CI " , c("L", "U"))
      z.crit <- abs(qnorm(alpha/2))
      means_ci[,1] <- means - means_se*z.crit
      means_ci[,2] <- means + means_se*z.crit

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

  model_mat <- aliased <- NULL
  if (model) {
    model_se <- sqrt(diag(object$vcov))
    model_tval <- model_coefs/model_se

    model_mat <- cbind(Estimate = model_coefs,
                       `Std. Error` = model_se,
                       `z value` = model_tval,
                       `Pr(>|z|)` = 2 * pnorm(abs(model_tval), lower.tail = FALSE))
    aliased <- is.na(model_coefs)

    if (ci) {
      conf <- matrix(nrow = nrow(model_mat), ncol = 2)
      colnames(conf) <- paste0(round(100*(1-alpha), 1), "% CI " , c("L", "U"))
      z.crit <- abs(qnorm(alpha/2))
      conf[,1] <- model_coefs - model_se*z.crit
      conf[,2] <- model_coefs + model_se*z.crit

      model_mat <- cbind(model_mat[,1:2,drop=FALSE], conf, model_mat[,-(1:2),drop=FALSE])
    }
  }

  ans <- list(call = object$call,
              means = means_mat,
              coefficients = effects_mat,
              model.coefficients = model_mat,
              aliased = aliased,
              sigma = sigma,
              r.squared = r.squared,
              adj.r.squared = adj.r.squared,
              estimand = object$estimand,
              focal = object$focal,
              treat_levels = object$treat_levels,
              fixef_name = attr(object$fixef, "fixef_name"))

  class(ans) <- "summary.lmw_est"
  return(ans)
}
