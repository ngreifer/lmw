summary.lmw <- function(object, un = TRUE, addlvariables = NULL, standardize = TRUE, data = NULL, stat = "balance", ...) {
  #Balance assessment, similar in structure to summary.matchit()

  stat <- match_arg(stat, c("balance", "distribution"))

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

  X_target <- target.weights <- NULL
  if (!is.null(object$target)) {
    X_target <- covs_df_to_matrix(model.frame(remove_treat_from_formula(delete.response(terms(object$formula)), attr(object$treat, "treat_name")),
                                              data = object$target))
    target.weights <- attr(object$target, "target.weights")

    if (length(addl_diff <- setdiff(colnames(X), colnames(X_target))) > 0) {
      addl_target <- matrix(NA_real_, nrow = nrow(X_target), ncol = length(addl_diff),
                            dimnames = list(rownames(X_target), addl_diff))
      X_target <- cbind(X_target, addl_target)
    }
  }

  treat <- apply_contrast_to_treat(object$treat, object$contrast)
  focal <- if (!is.null(object$focal)) setNames(c("Control", "Treated"), levels(treat))[object$focal]
  levels(treat) <- c("Control", "Treated")

  weights <- object$weights
  s.weights <- if (is.null(object$s.weights)) rep(1, length(weights)) else object$s.weights

  kk <- ncol(X)

  bal.un <- bal.base.weighted <- bal.weighted <- NULL
  dist.un <- dist.base.weighted <- dist.weighted <- NULL

  if (kk > 0) {
    #Remove tics (``) from outside of variables names
    nam <- colnames(X)
    has_tics <- startsWith(nam, "`") & endsWith(nam, "`")
    nam[has_tics] <- substr(nam[has_tics], 2, nchar(nam[has_tics]) - 1)

    if (stat == "balance") {
      #Compute balance statistics - SMD, KS
      if (un) {
        aa.un <- lapply(colnames(X), function(i) {
          balance_one_var(X[,i], treat = treat, weights = s.weights, s.weights = s.weights,
                          standardize = standardize, focal = focal,
                          x_target = X_target[,i], target.weights = target.weights)
        })

        bal.un <- do.call("rbind", aa.un)
        rownames(bal.un) <- nam

        if (!is.null(object$base.weights)) {
          aa.base.weighted <- lapply(colnames(X), function(i) {
            balance_one_var(X[,i], treat = treat,
                            weights = s.weights*object$base.weights,
                            s.weights = s.weights,
                            standardize = standardize, focal = focal,
                            x_target = X_target[,i], target.weights = target.weights)
          })

          bal.base.weighted <- do.call("rbind", aa.base.weighted)
          rownames(bal.base.weighted) <- nam
        }
      }

      aa.weighted <- lapply(colnames(X), function(i) {
        balance_one_var(X[,i], treat = treat, weights = weights, s.weights = s.weights,
                        standardize = standardize, focal = focal,
                        x_target = X_target[,i], target.weights = target.weights)
      })
      bal.weighted <- do.call("rbind", aa.weighted)
      rownames(bal.weighted) <- nam
    }
    else if (stat == "distribution") {
      #Compute distribution statistics - mean, SD
      if (un) {
        aa.un <- lapply(colnames(X), function(i) {
          distribution_one_var(X[,i], treat = treat, weights = s.weights, s.weights = s.weights)
        })

        dist.un <- do.call("rbind", aa.un)
        colnames(dist.un)[1:2] <- c("Mean Overall", "SD Overall")
        rownames(dist.un) <- nam

        if (!is.null(object$base.weights)) {
          aa.base.weighted <- lapply(colnames(X), function(i) {
            distribution_one_var(X[,i], treat = treat,
                                 weights = s.weights*object$base.weights,
                                 s.weights = s.weights, focal = focal,
                                 x_target = X_target[,i], target.weights = target.weights)
          })

          dist.base.weighted <- do.call("rbind", aa.base.weighted)
          rownames(dist.base.weighted) <- nam
        }
      }

      aa.weighted <- lapply(colnames(X), function(i) {
        distribution_one_var(X[,i], treat = treat, weights = weights, s.weights = s.weights,
                             focal = focal, x_target = X_target[,i],
                             target.weights = target.weights)
      })
      dist.weighted <- do.call("rbind", aa.weighted)
      rownames(dist.weighted) <- nam
    }
  }

  ## Sample sizes
  nn_w <- nn(treat, weights, object$base.weights, s.weights)

  ## output
  res <- list(call = object$call,
              nn = nn_w,
              bal.un = bal.un,
              bal.base.weighted = bal.base.weighted,
              bal.weighted = bal.weighted,
              dist.un = dist.un,
              dist.base.weighted = dist.base.weighted,
              dist.weighted = dist.weighted,
              method = object$method,
              base.weights.origin = attr(object$base.weights, "origin"))
  class(res) <- "summary.lmw"
  return(res)
}

summary.lmw_multi <- function(object, un = TRUE, addlvariables = NULL, standardize = TRUE, data = NULL,
                              contrast = NULL, stat = "balance", ...) {
  #Balance assessment, similar in structure to summary.matchit()

  stat <- match_arg(stat, c("balance", "distribution"))

  if (object$method == "URI") {
    if (hasName(match.call(), "contrast")) {
      warning("The 'contrast' argument is ignored when method = \"URI\" in the original call to lmw().", call. = FALSE)
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

  X_target <- target.weights <- NULL
  if (!is.null(object$target)) {
    X_target <- covs_df_to_matrix(object$target)
    target.weights <- attr(object$target, "target.weights")
  }

  treat <- object$treat
  weights <- object$weights
  s.weights <- if (is.null(object$s.weights)) rep(1, length(weights)) else object$s.weights

  if (is.null(contrast)) {
    un_weights <- s.weights
    if (stat == "balance") {
      balance_fun <- balance_one_var.multi
    }
  }
  else {
    contrast <- process_contrast(contrast, treat, object$method)
    un_weights <- s.weights * (treat %in% contrast)
    treat <- apply_contrast_to_treat(treat, contrast)
    if (object$method == "MRI") weights <- weights * (treat %in% contrast)
    else treat[treat != levels(treat)[2]] <- levels(treat)[1]
    if (stat == "balance") {
      balance_fun <- balance_one_var
    }
  }

  kk <- ncol(X)

  bal.un <- bal.base.weighted <- bal.weighted <- NULL
  dist.un <- dist.base.weighted <- dist.weighted <- NULL

  ## Summary Stats
  if (kk > 0) {
    #Remove tics (``) from outside of variables names
    nam <- colnames(X)
    has_tics <- startsWith(nam, "`") & endsWith(nam, "`")
    nam[has_tics] <- substr(nam[has_tics], 2, nchar(nam[has_tics]) - 1)

    if (stat == "balance") {
      if (un) {
        aa.un <- lapply(colnames(X), function(i) {
          balance_fun(X[,i], treat = treat, weights = un_weights, s.weights = s.weights,
                      standardize = standardize, focal = object$focal,
                      x_target = X_target[,i], target.weights = target.weights)
        })

        bal.un <- do.call("rbind", aa.un)
        rownames(bal.un) <- nam

        if (!is.null(object$base.weights)) {

          aa.base.weighted <- lapply(colnames(X), function(i) {
            balance_fun(X[,i], treat = treat,
                        weights = un_weights*object$base.weights,
                        s.weights = s.weights,
                        standardize = standardize, focal = object$focal,
                        x_target = X_target[,i], target.weights = target.weights)
          })

          bal.base.weighted <- do.call("rbind", aa.base.weighted)
          rownames(bal.base.weighted) <- nam
        }
      }

      aa.weighted <- lapply(colnames(X), function(i) {
        balance_fun(X[,i], treat = treat, weights = weights, s.weights = s.weights,
                    standardize = standardize, focal = object$focal,
                    x_target = X_target[,i], target.weights = target.weights)
      })
      bal.weighted <- do.call("rbind", aa.weighted)
      rownames(bal.weighted) <- nam
    }
    else if (stat == "distribution") {
      #Compute distribution statistics - mean, SD
      if (un) {
        aa.un <- lapply(colnames(X), function(i) {
          distribution_one_var(X[,i], treat = treat, weights = un_weights, s.weights = s.weights,
                               contrast = contrast)
        })

        dist.un <- do.call("rbind", aa.un)
        colnames(dist.un)[1:2] <- c("Mean Overall", "SD Overall")
        rownames(dist.un) <- nam

        if (!is.null(object$base.weights)) {
          aa.base.weighted <- lapply(colnames(X), function(i) {
            distribution_one_var(X[,i], treat = treat,
                                 weights = s.weights*object$base.weights,
                                 s.weights = s.weights, focal = object$focal,
                                 x_target = X_target[,i], target.weights = target.weights,
                                 contrast = contrast)
          })

          dist.base.weighted <- do.call("rbind", aa.base.weighted)
          rownames(dist.base.weighted) <- nam
        }
      }

      aa.weighted <- lapply(colnames(X), function(i) {
        distribution_one_var(X[,i], treat = treat, weights = weights, s.weights = s.weights,
                             focal = object$focal, x_target = X_target[,i],
                             target.weights = target.weights,
                             contrast = contrast)
      })
      dist.weighted <- do.call("rbind", aa.weighted)
      rownames(dist.weighted) <- nam
    }
  }

  ## Sample sizes
  nn_w <- {
    if (is.null(contrast)) nn_multi(treat, weights, object$base.weights, s.weights)
    else nn(treat, weights, object$base.weights, s.weights)
  }

  ## output
  res <- list(call = object$call, nn = nn_w,
              bal.un = bal.un,
              bal.base.weighted = bal.base.weighted,
              bal.weighted = bal.weighted,
              dist.un = dist.un,
              dist.base.weighted = dist.base.weighted,
              dist.weighted = dist.weighted,
              method = object$method,
              base.weights.origin = attr(object$base.weights, "origin"))
  class(res) <- c("summary.lmw_multi"[is.null(contrast)], "summary.lmw")
  return(res)
}

print.summary.lmw <- function(x, digits = max(3, getOption("digits") - 4), ...){

  if (!is.null(x$call)) cat("\nCall:", deparse(x$call), sep = "\n")

  if (!is.null(x$bal.un)) {
    cat("\nSummary of Balance for Unweighted Data:\n")
    print.data.frame(round_df_char(x$bal.un, digits))
  }

  if (!is.null(x$bal.base.weighted)) {
    cat("\nSummary of Balance for Base Weighted Data:\n")
    print.data.frame(round_df_char(x$bal.base.weighted, digits))
  }

  if (!is.null(x$bal.weighted)) {
    cat("\nSummary of Balance for Weighted Data:\n")
    print.data.frame(round_df_char(x$bal.weighted, digits))
  }

  if (!is.null(x$dist.un)) {
    cat("\nDistribution Summary for Unweighted Data:\n")
    print.data.frame(add_peren_to_sd(round_df_char(x$dist.un, digits)))
  }

  if (!is.null(x$dist.base.weighted)) {
    cat("\nDistribution Summary for Base Weighted Data:\n")
    print.data.frame(add_peren_to_sd(round_df_char(x$dist.base.weighted, digits)))
  }

  if (!is.null(x$dist.weighted)) {
    cat("\nDistribution Summary for Weighted Data:\n")
    print.data.frame(add_peren_to_sd(round_df_char(x$dist.weighted, digits)))
  }

  if (!is.null(x$nn)) {
    cat("\nEffective Sample Sizes:\n")
    print.data.frame(round_df_char(x$nn, 2, pad = " "))
  }
  cat("\n")
  invisible(x)
}

balance_one_var <- function(x, treat, weights, s.weights, standardize = TRUE,
                            focal = NULL, x_target = NULL, target.weights = NULL) {
  #weights must already have s.weights applied, which is true of regression weights from
  #lmw() but but not base.weights, so make sure they are multiplied by s.weights in the
  #function call.

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
    else if (!is.null(x_target)) mean_w(x_target, target.weights)
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
        if (!is.null(focal)) sqrt(var_w(x, bin.var, s.weights, treat==focal))
        else if (!is.null(x_target) && length(x_target) > 1) sqrt(var_w(x_target, bin.var, target.weights))
        else sqrt(mean(vapply(levels(treat), function(t) var_w(x, bin.var, s.weights, treat==t), numeric(1L))))
      }

      if (std < sqrt(.Machine$double.eps)) std <- sqrt(var_w(x, bin.var, s.weights)) #Avoid divide by zero

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
    else if (!is.null(x_target)) {
      if (is.null(target.weights)) target.weights <- rep(1, length(x_target))
      xsum[5] <- tks_w(x[treat != t1], x_target, weights[treat != t1], target.weights)
      xsum[6] <- tks_w(x[treat == t1], x_target, weights[treat == t1], target.weights)
    }
    else {
      xsum[5] <- tks_w(x[treat != t1], x, weights[treat != t1], s.weights)
      xsum[6] <- tks_w(x[treat == t1], x, weights[treat == t1], s.weights)
    }
  }

  xsum
}
balance_one_var.multi <- function(x, treat, weights = NULL, s.weights, standardize = TRUE,
                                  focal = NULL, x_target = NULL, target.weights = NULL) {
  #weights must already have s.weights applied, which is true of regression weights from
  #lmw() but but not base.weights, so make sure they are multiplied by s.weights in the
  #function call.

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
    else if (!is.null(x_target)) mean_w(x_target, target.weights)
    else mean_w(x, s.weights)
  }

  Mt <- vapply(levels(treat), function(t) mean_w(x, weights, treat==t), numeric(1L))

  mdifft <- Mt - M

  if (standardize) {
    if (!too.small) {
      std <- {
        if (!is.null(focal)) sqrt(var_w(x, bin.var, s.weights, treat==focal))
        else if (!is.null(x_target) && length(x_target) > 1) sqrt(var_w(x_target, bin.var, target.weights))
        else sqrt(mean(vapply(levels(treat), function(t) var_w(x, bin.var, s.weights, treat==t), numeric(1L))))
      }

      if (std < sqrt(.Machine$double.eps)) std <- sqrt(var_w(x, bin.var, s.weights)) #Avoid divide by zero

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
      else if (!is.null(x_target)) {
        if (is.null(target.weights)) target.weights <- rep(1, length(x_target))
        vapply(levels(treat), function(t) {
          tks_w(x[treat == t], x_target, weights[treat == t], target.weights)
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

distribution_one_var <- function(x, treat, weights, s.weights, focal = NULL,
                                 x_target = NULL, target.weights = NULL, contrast = NULL) {
  #weights must already have s.weights applied, which is true of regression weights from
  #lmw() but but not base.weights, so make sure they are multiplied by s.weights in the
  #function call.
  bin.var <- all(x == 0 | x == 1)

  tlevs <- levels(treat)
  if (!is.null(contrast)) {
    tlevs <- intersect(tlevs, contrast)
  }

  xsum <- c(
    "Mean Target" = {
      if (!is.null(focal)) mean_w(x, s.weights, treat==focal)
      else if (!is.null(x_target)) mean_w(x_target, target.weights)
      else mean_w(x, s.weights)
    },
    "SD Target" = {
      if (!is.null(focal)) sqrt(var_w(x, bin.var, s.weights, treat==focal))
      else if (!is.null(x_target)) {
        if (length(x_target) > 1) sqrt(var_w(x_target, bin.var, target.weights))
        else NA_real_
      }
      else sqrt(var_w(x, bin.var, s.weights))
    },
    unlist(lapply(tlevs, function(t) {
      setNames(c(mean_w(x, weights, treat==t),
                 sqrt(var_w(x, bin.var, weights, treat==t))),
               paste(c("Mean", "SD"), t))
    }))
  )

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
tks_w <- function(x, x_target, weights, target.weights) {
  treat <- factor(c(rep(1, length(x)), rep(0, length(x_target))), levels = 0:1)
  x <- c(x, x_target)
  weights <- c(weights, target.weights)

  ks_w(x, treat, weights)
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
