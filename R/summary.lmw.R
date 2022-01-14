summary.lmw <- function(object, un = TRUE, addlvariables = NULL, standardize = TRUE, data = NULL, ...) {
  #Balance assessment, similar in structure to summary.matchit()

  if (is.null(object$covs)) {
    X <- matrix(nrow = length(object$treat), ncol = 0)
  }
  else {
    X <- covs_df_to_matrix(object$covs)
  }

  if (!is.null(addlvariables)) {

    data <- get_data(data, object)

    if (is.character(addlvariables)) {
      if (!is.null(data) && is.data.frame(data)) {
        if (all(addlvariables %in% names(data))) {
          addlvariables <- data[addlvariables]
        }
        else {
          stop("All variables in 'addlvariables' must be in 'data'.", call. = FALSE)
        }
      }
      else {
        stop("If 'addlvariables' is specified as a string, a data frame argument must be supplied to 'data'.", call. = FALSE)
      }
    }
    else if (inherits(addlvariables, "formula")) {
      vars.in.formula <- all.vars(addlvariables)
      if (!is.null(data) && is.data.frame(data)) data <- cbind(data[names(data) %in% vars.in.formula],
                                                               object$covs[names(data) %in% setdiff(vars.in.formula, names(data))])
      else data <- object$covs

      addlvariables <- covs_df_to_matrix(model.frame(addlvariables, data = data))
    }
    else if (!is.matrix(addlvariables) && !is.data.frame(addlvariables)) {
      stop("The argument to 'addlvariables' must be in one of the accepted forms. See ?summary.lmw for details.", call. = FALSE)
    }

    if (is.data.frame(addlvariables)) {
      addlvariables <- covs_df_to_matrix(addlvariables)
    }

    X <- cbind(X, addlvariables[, setdiff(colnames(addlvariables), colnames(X)), drop = FALSE])
  }

  X_target <- NULL
  if (!is.null(object$target)) {
    formula_without_treat <- remove_treat_from_formula(object$formula, attr(object$treat, "treat_name"))
    target <- process_target(object$target, formula_without_treat, object$covs)
    X_target <- covs_df_to_matrix(attr(target, "target_original"))
  }

  treat <- apply_contrast_to_treat(object$treat, object$contrast)
  focal <- if (!is.null(object$focal)) setNames(c("Control", "Treated"), levels(treat))[object$focal]
  levels(treat) <- c("Control", "Treated")

  weights <- object$weights
  s.weights <- if (is.null(object$s.weights)) rep(1, length(weights)) else object$s.weights

  kk <- ncol(X)

  if (kk > 0) {
    nam <- colnames(X)
    nam[startsWith(nam, "`") & endsWith(nam, "`")] <- substr(nam[startsWith(nam, "`") & endsWith(nam, "`")],
                                                             2, nchar(nam[startsWith(nam, "`") & endsWith(nam, "`")]) - 1)
  }

  sum.un <- sum.base.weighted <- sum.weighted <- NULL
  ## Summary Stats
  if (kk > 0) {
    if (un) {
      aa.un <- lapply(colnames(X), function(i) {
        balance_one_var(i, X, treat = treat, weights = s.weights, s.weights = s.weights,
                    standardize = standardize, focal = focal,
                    target = X_target)
      })

      sum.un <- do.call("rbind", aa.un)
      rownames(sum.un) <- nam

      if (!is.null(object$base.weights)) {
        aa.base.weighted <- lapply(colnames(X), function(i) {
          balance_one_var(i, X, treat = treat,
                      weights = s.weights*object$base.weights,
                      s.weights = s.weights,
                      standardize = standardize, focal = focal,
                      target = X_target)
        })

        sum.base.weighted <- do.call("rbind", aa.base.weighted)
        rownames(sum.base.weighted) <- nam
      }
    }

    aa.weighted <- lapply(colnames(X), function(i) {
      balance_one_var(i, X, treat = treat, weights = weights, s.weights = s.weights,
                  standardize = standardize, focal = focal,
                  target = X_target)
    })
    sum.weighted <- do.call("rbind", aa.weighted)
    rownames(sum.weighted) <- nam
  }

  ## Sample sizes
  nn_w <- nn(treat, weights, object$base.weights, s.weights)

  ## output
  res <- list(call = object$call, nn = nn_w, sum.un = sum.un,
              sum.base.weighted = sum.base.weighted,
              sum.weighted = sum.weighted,
              type = object$type,
              base.weights.origin = attr(object$base.weights, "origin"))
  class(res) <- "summary.lmw"
  return(res)
}

summary.lmw.multi <- function(object, un = TRUE, addlvariables = NULL, standardize = TRUE, data = NULL, contrast = NULL, ...) {
  #Balance assessment, similar in structure to summary.matchit()

  if (object$type == "URI") {
    if (hasName(match.call(), "contrast")) {
      warning("The 'contrast' argument is ignored when type = \"URI\" in the original call to lmw().", call. = FALSE)
    }
    contrast <- object$contrast
  }
  else if (!hasName(match.call(), "contrast")) {
    contrast <- object$contrast
  }

  if (is.null(object$covs)) {
    X <- matrix(nrow = length(object$treat), ncol = 0)
  }
  else {
    X <- covs_df_to_matrix(object$covs)
  }

  if (!is.null(addlvariables)) {

    data <- get_data(data, object)

    if (is.character(addlvariables)) {
      if (!is.null(data) && is.data.frame(data)) {
        if (all(addlvariables %in% names(data))) {
          addlvariables <- data[addlvariables]
        }
        else {
          stop("All variables in 'addlvariables' must be in 'data'.", call. = FALSE)
        }
      }
      else {
        stop("If 'addlvariables' is specified as a string, a data frame argument must be supplied to 'data'.", call. = FALSE)
      }
    }
    else if (inherits(addlvariables, "formula")) {
      vars.in.formula <- all.vars(addlvariables)
      if (!is.null(data) && is.data.frame(data)) data <- cbind(data[names(data) %in% vars.in.formula],
                                                               object$covs[names(data) %in% setdiff(vars.in.formula, names(data))])
      else data <- object$covs

      addlvariables <- covs_df_to_matrix(model.frame(addlvariables, data = data))
    }
    else if (!is.matrix(addlvariables) && !is.data.frame(addlvariables)) {
      stop("The argument to 'addlvariables' must be in one of the accepted forms. See ?summary.lmw for details.", call. = FALSE)
    }

    if (is.data.frame(addlvariables)) {
      addlvariables <- covs_df_to_matrix(addlvariables)
    }

    X <- cbind(X, addlvariables[, setdiff(colnames(addlvariables), colnames(X)), drop = FALSE])
  }

  X_target <- NULL
  if (!is.null(object$target)) {
    formula_without_treat <- remove_treat_from_formula(object$formula, attr(object$treat, "treat_name"))
    target <- process_target(object$target, formula_without_treat, object$covs)
    X_target <- covs_df_to_matrix(attr(target, "target_original"))
  }

  treat <- object$treat
  weights <- object$weights
  s.weights <- if (is.null(object$s.weights)) rep(1, length(weights)) else object$s.weights

  if (is.null(contrast)) {
    balance_fun <- balance_one_var.multi
    un_weights <- s.weights
  }
  else {
    contrast <- process_contrast(contrast, treat, object$type)
    balance_fun <- balance_one_var
    un_weights <- s.weights * (treat %in% contrast)
    treat <- apply_contrast_to_treat(treat, contrast)
    if (object$type == "MRI") weights <- weights * (treat %in% contrast)
  }

  kk <- ncol(X)

  if (kk > 0) {
    nam <- colnames(X)
    has_tics <- startsWith(nam, "`") & endsWith(nam, "`")
    nam[has_tics] <- substr(nam[has_tics], 2, nchar(nam[has_tics]) - 1)
  }

  sum.un <- sum.base.weighted <- sum.weighted <- NULL
  ## Summary Stats
  if (kk > 0) {
    if (un) {
      aa.un <- lapply(colnames(X), function(i) {
        balance_fun(i, X, treat = treat, weights = un_weights, s.weights = s.weights,
                    standardize = standardize, focal = object$focal,
                    target = X_target)
      })

      sum.un <- do.call("rbind", aa.un)
      rownames(sum.un) <- nam

      if (!is.null(object$base.weights)) {

        aa.base.weighted <- lapply(colnames(X), function(i) {
          balance_fun(i, X, treat = treat,
                      weights = un_weights*object$base.weights,
                      s.weights = s.weights,
                      standardize = standardize, focal = object$focal,
                      target = X_target)
        })

        sum.base.weighted <- do.call("rbind", aa.base.weighted)
        rownames(sum.base.weighted) <- nam
      }
    }

    aa.weighted <- lapply(colnames(X), function(i) {
      balance_fun(i, X, treat = treat, weights = weights, s.weights = s.weights,
                  standardize = standardize, focal = object$focal,
                  target = X_target)
    })
    sum.weighted <- do.call("rbind", aa.weighted)
    rownames(sum.weighted) <- nam
  }

  ## Sample sizes
  nn_w <- {
    if (is.null(contrast)) nn_multi(treat, weights, object$base.weights, s.weights)
    else nn(treat, weights, object$base.weights, s.weights)
  }

  ## output
  res <- list(call = object$call, nn = nn_w, sum.un = sum.un,
              sum.base.weighted = sum.base.weighted,
              sum.weighted = sum.weighted,
              type = object$type,
              base.weights.origin = attr(object$base.weights, "origin"))
  class(res) <- c("summary.lmw.multi"[is.null(contrast)], "summary.lmw")
  return(res)
}

print.summary.lmw <- function(x, digits = max(3, getOption("digits") - 3), ...){

  if (!is.null(x$call)) cat("\nCall:", deparse(x$call), sep = "\n")

  if (!is.null(x$sum.un)) {
    cat("\nSummary of Balance for Unweighted Data:\n")
    print.data.frame(round_df_char(x$sum.un[,-7, drop = FALSE], digits, pad = "0", na_vals = "."))
  }

  if (!is.null(x$sum.base.weighted)) {
    cat("\nSummary of Balance for Base Weighted Data:\n")
    print.data.frame(round_df_char(x$sum.base.weighted, digits, pad = "0", na_vals = "."))
  }

  if (!is.null(x$sum.weighted)) {
    cat("\nSummary of Balance for Weighted Data:\n")
    print.data.frame(round_df_char(x$sum.weighted, digits, pad = "0", na_vals = "."))
  }

  if (!is.null(x$nn)) {
    cat("\nEffective Sample Sizes:\n")
    print.data.frame(round_df_char(x$nn, 2, pad = " ", na_vals = "."))
  }
  cat("\n")
  invisible(x)
}

balance_one_var <- function(var, X, treat, weights, s.weights, standardize = TRUE, focal = NULL, target = NULL) {
  #weights must already have s.weights applied, which is true of regression weights from
  #lmw() but but not base.weights, so make sure they are multiplied by s.weights in the
  #function call.

  x <- X[,var]
  t1 <- levels(treat)[2] #treated level

  xsum <- rep(NA_real_, 6)
  if (standardize)
    names(xsum) <- c("SMD", paste("TSMD", levels(treat)[1:2]), "KS", paste("TKS", levels(treat)[1:2]))
  else
    names(xsum) <- c("MD", paste("TMD", levels(treat)[1:2]), "KS", paste("TKS", levels(treat)[1:2]))

  #If all variable values are the same, return 0s
  if (all(abs(x[-1] - x[1]) < sqrt(.Machine$double.eps))) {
    xsum[] <- 0
    return(xsum)
  }

  bin.var <- all(x == 0 | x == 1)

  too.small <- any(vapply(levels(treat), function(t) sum(s.weights[treat==t] != 0) < 2, logical(1L)))

  M <- {
    if (!is.null(focal)) mean_w(x, s.weights, treat==focal)
    else if (!is.null(target)) target[,var]
    else mean_w(x, s.weights)
  }

  M0 <- mean_w(x, weights, treat!=t1)
  M1 <- mean_w(x, weights, treat==t1)

  mdiff <- M1 - M0
  mdiff0 <- M0 - M
  mdiff1 <- M1 - M

  if (all(abs(c(mdiff, mdiff0, mdiff1)) < sqrt(.Machine$double.eps))) {
    xsum[1:3] <- 0
  }
  else if (standardize) {
    if (!too.small) {
      std <- {
        if (!is.null(focal)) sqrt(wvar(x[treat==focal], bin.var, s.weights[treat==focal]))
        else sqrt(mean(vapply(levels(treat), function(t) wvar(x[treat==t], bin.var, s.weights[treat==t]), numeric(1L))))
      }

      if (std < sqrt(.Machine$double.eps)) std <- sqrt(wvar(x, bin.var, s.weights)) #Avoid divide by zero

      xsum[1] <- mdiff/std
      xsum[2] <- mdiff0/std
      xsum[3] <- mdiff1/std
    }
  }
  else {
    xsum[1] <- mdiff
    xsum[2] <- mdiff0
    xsum[3] <- mdiff1
  }

  if (bin.var) {
    xsum[4] <- abs(mdiff)
    xsum[5] <- abs(mdiff0)
    xsum[6] <- abs(mdiff1)
  }
  else if (!too.small) {
    xsum[4] <- ks_w(x, treat, weights)
    if (!is.null(focal)) {
      xsum[5] <- tks_w(x[treat != t1], x[treat == focal], weights[treat != t1], s.weights[treat == focal])
      xsum[6] <- tks_w(x[treat == t1], x[treat == focal], weights[treat == t1], s.weights[treat == focal])
    }
    else if (!is.null(target)) {
      xsum[5] <- tks_w(x[treat != t1], target[,var], weights[treat != t1], 1)
      xsum[6] <- tks_w(x[treat == t1], target[,var], weights[treat == t1], 1)
    }
    else {
      xsum[5] <- tks_w(x[treat != t1], x, weights[treat != t1], s.weights)
      xsum[6] <- tks_w(x[treat == t1], x, weights[treat == t1], s.weights)
    }
  }

  xsum
}
balance_one_var.multi <- function(var, X, treat, weights = NULL, s.weights, standardize = TRUE, focal = NULL, target = NULL) {
  #weights must already have s.weights applied, which is true of regression weights from
  #lmw() but but not base.weights, so make sure they are multiplied by s.weights in the
  #function call.

  x <- X[,var]

  xsum <- rep(NA_real_, 6)
  if (standardize)
    names(xsum) <- c(paste("TSMD", levels(treat)), paste("TKS", levels(treat)))
  else
    names(xsum) <- c(paste("TMD", levels(treat)), paste("TKS", levels(treat)))

  #If all variable values are the same, return 0s
  if (all(abs(x[-1] - x[1]) < sqrt(.Machine$double.eps))) {
    xsum[] <- 0
    return(xsum)
  }

  bin.var <- all(x == 0 | x == 1)

  too.small <- any(vapply(levels(treat), function(t) sum(s.weights[treat==t] != 0) < 2, logical(1L)))

  M <- {
    if (!is.null(focal)) mean_w(x, s.weights, treat==focal)
    else if (!is.null(target)) target[,var]
    else mean_w(x, s.weights)
  }

  Mt <- vapply(levels(treat), function(t) mean_w(x, weights, treat==t), numeric(1L))

  mdifft <- Mt - M

  if (standardize) {
    if (!too.small) {
      std <- {
        if (!is.null(focal)) sqrt(wvar(x[treat==focal], bin.var, s.weights[treat==focal]))
        else sqrt(mean(vapply(levels(treat), function(t) wvar(x[treat==t], bin.var, s.weights[treat==t]), numeric(1L))))
      }

      if (std < sqrt(.Machine$double.eps)) std <- sqrt(wvar(x, bin.var, s.weights)) #Avoid divide by zero

      xsum[seq_len(nlevels(treat))] <- mdifft/std
    }
  }
  else {
    xsum[seq_len(nlevels(treat))] <- mdifft
  }

  if (bin.var) {
    xsum[nlevels(treat) + seq_len(nlevels(treat))] <- abs(mdifft)
  }
  else if (!too.small) {

    xsum[nlevels(treat) + seq_len(nlevels(treat))] <- {
      if (!is.null(focal)) {
        vapply(levels(treat), function(t) {
          tks_w(x[treat == t], x[treat == focal], weights[treat == t], s.weights[treat == focal])
        }, numeric(1L))
      }
      else if (!is.null(target)) {
        vapply(levels(treat), function(t) {
          tks_w(x[treat == t], target[,var], weights[treat == t], 1)
        }, numeric(1L))
      }
      else {
        vapply(levels(treat), function(t) {
          tks_w(x[treat == t], x, weights[treat == t], s.weights)
        }, numeric(1L))
      }
    }
  }

  xsum
}

ks_w <- function(x, treat, weights) {
  t1 <- levels(treat)[2]
  weights[treat==t1] <- weights[treat==t1]/sum(weights[treat==t1])
  weights[treat!=t1] <- weights[treat!=t1]/sum(weights[treat!=t1])

  ord <- order(x)
  x_ord <- x[ord]
  weights_ord <- weights[ord]
  treat_ord <- treat[ord]

  #Difference between ecdf of x for each group
  weights_ord_ <- weights_ord
  weights_ord_[treat_ord==t1] <- -weights_ord_[treat_ord==t1]
  ediff <- abs(cumsum(weights_ord_))[c(diff(x_ord) != 0, TRUE)]

  max(ediff)
}
tks_w <- function(x, x.target, weights, weights.target) {
  treat <- factor(c(rep(1, length(x)), rep(0, length(x.target))), levels = 0:1)
  x <- c(x, x.target)
  weights <- c(weights, weights.target)

  ks_w(x, treat, weights)
}

#(Weighted) variance that uses special formula for binary variables
wvar <- function(x, bin.var = NULL, w = NULL) {
  if (is.null(w)) w <- rep(1, length(x))
  if (is.null(bin.var)) bin.var <- all(x == 0 | x == 1)

  w <- w / sum(w) #weights normalized to sum to 1
  mx <- sum(w * x) #weighted mean

  if (bin.var) {
    mx*(1-mx)
  }
  else {
    #Reliability weights variance; same as cov.wt()
    sum(w * (x - mx)^2)/(1 - sum(w^2))
  }
}

#Compute sample sizes
nn <- function(treat, weights, base.weights, s.weights) {
  # weights <- weights * s.weights #s.weights already in weights

  t1 <- levels(treat)[2] #treated level
  t0 <- levels(treat)[1] #control level

  if (is.null(base.weights)) {
    n <- matrix(0, nrow=2, ncol=2, dimnames = list(c("All", "Weighted"),
                                                   levels(treat)[1:2]))

    n["All",] <- c(ESS(s.weights[treat==t0]),
                   ESS(s.weights[treat==t1]))
    n["Weighted",] <- c(ESS(weights[treat!=t1]),
                        ESS(weights[treat==t1]))
  }
  else {
    base.weights <- base.weights*s.weights
    n <- matrix(0, nrow=3, ncol=2, dimnames = list(c("All","Base weighted", "Weighted"),
                                                   levels(treat)[1:2]))

    n["All",] <- c(ESS(s.weights[treat==t0]),
                   ESS(s.weights[treat==t1]))
    n["Base weighted",] <- c(ESS((s.weights*base.weights)[treat==t0]),
                             ESS((s.weights*base.weights)[treat==t1]))
    n["Weighted",] <- c(ESS(weights[treat!=t1]),
                        ESS(weights[treat==t1]))
  }

  return(n)
}
nn_multi <- function(treat, weights, base.weights, s.weights) {
  if (is.null(s.weights)) s.weights <- rep(1, length(treat))
  # weights <- weights * s.weights #s.weights already in weights

  if (is.null(base.weights)) {
    n <- matrix(0, nrow=2, ncol=nlevels(treat), dimnames = list(c("All", "Weighted"),
                                                                levels(treat)))

    n["All",] <-     vapply(levels(treat), function(t) ESS(s.weights[treat==t]), numeric(1L))
    n["Weighted",] <- vapply(levels(treat), function(t) ESS(weights[treat==t]), numeric(1L))
  }
  else {
    base.weights <- base.weights*s.weights
    n <- matrix(0, nrow=3, ncol=nlevels(treat), dimnames = list(c("All","Base weighted", "Weighted"),
                                                                levels(treat)))

    n["All",] <-     vapply(levels(treat), function(t) ESS(s.weights[treat==t]), numeric(1L))
    n["Base weighted",] <- vapply(levels(treat), function(t) ESS(base.weights[treat==t]), numeric(1L))
    n["Weighted",] <- vapply(levels(treat), function(t) ESS(weights[treat==t]), numeric(1L))
  }

  return(n)
}
