summary.lmw <- function(object, un = TRUE, addlvariables = NULL, standardize = TRUE, data = NULL, ...) {
  #Balance assessment
  #Want to be similar to summary.matchit(), but use TASMD as primary metric

  if (is.null(object$covs)) {
    X <- matrix(nrow = length(object$treat), ncol = 0)
  }
  else {
    X <- covs_df_to_matrix(object$covs)
  }

  if (!is.null(addlvariables)) {
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

  treat <- object$treat
  weights <- object$weights
  s.weights <- if (is.null(object$s.weights)) rep(1, length(weights)) else object$s.weights

  kk <- ncol(X)

  if (kk > 0) {
    nam <- colnames(X)
    nam[startsWith(nam, "`") & endsWith(nam, "`")] <- substr(nam[startsWith(nam, "`") & endsWith(nam, "`")],
                                                             2, nchar(nam[startsWith(nam, "`") & endsWith(nam, "`")]) - 1)
  }

  ## Summary Stats
  if (kk > 0) {
    if (un) {
      aa.all <- setNames(lapply(seq_len(kk), function(i) balance_one_var(X[,i], treat = treat, weights = NULL, s.weights = s.weights,
                                                                         standardize = standardize, estimand = object$estimand)),
                         colnames(X))
      sum.un <- do.call("rbind", aa.all)
      dimnames(sum.un) <- list(nam, names(aa.all[[1]]))
    }

    aa.matched <- setNames(lapply(seq_len(kk), function(i) balance_one_var(X[,i], treat = treat, weights = weights, s.weights = s.weights,
                                                                           standardize = standardize, estimand = object$estimand)),
                           colnames(X))
    sum.weighted <- do.call("rbind", aa.matched)
    dimnames(sum.weighted) <- list(nam, names(aa.matched[[1]]))
  }

  ## Sample sizes
  nn_w <- nn(treat, weights, s.weights)

  ## output
  res <- list(call = object$call, nn = nn_w, sum.un = if (un) sum.un,
              sum.weighted = sum.weighted)
  class(res) <- "summary.lmw"
  return(res)
}

print.summary.lmw <- function(x, digits = max(3, getOption("digits") - 3), ...){

  if (!is.null(x$call)) cat("\nCall:", deparse(x$call), sep = "\n")

  if (!is.null(x$sum.un)) {
    cat("\nSummary of Balance for All Data:\n")
    print.data.frame(round_df_char(x$sum.un[,-7, drop = FALSE], digits, pad = "0", na_vals = "."))
    cat("\n")
  }

  if (!is.null(x$sum.weighted)) {
    cat("\nSummary of Balance for Weighted Data:\n")
    print.data.frame(round_df_char(x$sum.weighted, digits, pad = "0", na_vals = "."))
  }

  if (!is.null(x$nn)) {
    cat("\nEffective sample Sizes:\n")
    print.data.frame(round_df_char(x$nn, 2, pad = " ", na_vals = "."))
  }
  cat("\n")
  invisible(x)
}

balance_one_var <- function(x, treat, weights = NULL, s.weights, standardize = TRUE, estimand) {

  un <- is.null(weights)
  bin.var <- all(x == 0 | x == 1)

  xsum <- rep(NA_real_, 6)
  if (standardize)
    names(xsum) <- c("SMD", "TSMD Treated", "TSMD Control", "KS", "TKS Treated", "TKS Control")
  else
    names(xsum) <- c("MD", "TMD Treated", "TMD Control", "KS", "TKS Treated", "TKS Control")

  if (un) weights <- s.weights
  # else weights <- weights * s.weights #s.weights already multiplied in

  too.small <- sum(weights[treat==1] != 0) < 2 || sum(weights[treat==0] != 0) < 2

  M <- switch(estimand,
              "ATT" = mean_w(x, s.weights, treat==1),
              "ATC" = mean_w(x, s.weights, treat==0),
              "ATE" = mean_w(x, s.weights))

  M1 <- mean_w(x, weights, treat==1)
  M0 <- mean_w(x, weights, treat==0)

  mdiff <- M1 - M0
  mdiff1 <- M1 - M
  mdiff0 <- M - M0

  if (standardize) {
    if (!too.small) {
        std <- switch(estimand,
                      "ATT" = sqrt(wvar(x[treat==1], bin.var, s.weights[treat==1])),
                      "ATC" = sqrt(wvar(x[treat==0], bin.var, s.weights[treat==0])),
                      sqrt(.5*(wvar(x[treat==1], bin.var, s.weights[treat==1]) + wvar(x[treat==0], bin.var, s.weights[treat==0]))))

        if (std < sqrt(.Machine$double.eps)) std <- sqrt(wvar(x, bin.var, s.weights)) #Avoid divide by zero

      xsum[1] <- mdiff/std
      xsum[2] <- mdiff1/std
      xsum[3] <- mdiff0/std
    }
  }
  else {
    xsum[1] <- mdiff
    xsum[2] <- mdiff1
    xsum[3] <- mdiff0
  }

  if (bin.var) {
    xsum[4] <- abs(mdiff)
    xsum[5] <- abs(mdiff1)
    xsum[6] <- abs(mdiff0)
  }
  else if (!too.small) {
    xsum[4] <- ks_w(x, treat, weights)

    if (estimand == "ATE") {
      xsum[5] <- tks_w(x[treat == 1], x, weights[treat == 1], s.weights)
      xsum[6] <- tks_w(x[treat == 0], x, weights[treat == 0], s.weights)
    }
    else if (estimand == "ATT") {
      xsum[5] <- tks_w(x[treat == 1], x[treat == 1], weights[treat == 1], s.weights[treat == 1])
      xsum[6] <- tks_w(x[treat == 0], x[treat == 1], weights[treat == 0], s.weights[treat == 1])
    }
    else if (estimand == "ATC") {
      xsum[5] <- tks_w(x[treat == 1], x[treat == 0], weights[treat == 1], s.weights[treat == 0])
      xsum[6] <- tks_w(x[treat == 0], x[treat == 0], weights[treat == 0], s.weights[treat == 0])
    }
  }

  xsum
}

ks_w <- function(x, treat, weights) {
  for (i in unique(treat)) weights[treat==i] <- weights[treat==i]/sum(weights[treat==i])

  ord <- order(x)
  x_ord <- x[ord]
  weights_ord <- weights[ord]
  treat_ord <- treat[ord]

  #Difference between ecdf of x for each group
  weights_ord_ <- weights_ord
  weights_ord_[treat_ord==treat_ord[1]] <- -weights_ord_[treat_ord==treat_ord[1]]
  ediff <- abs(cumsum(weights_ord_))[c(diff(x_ord) != 0, TRUE)]

  max(ediff)
}
tks_w <- function(x, x.target, weights, weights.target) {
  treat <- c(rep(1, length(x)), rep(0, length(x.target)))
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
nn <- function(treat, weights, s.weights) {
  if (is.null(s.weights)) s.weights <- rep(1, length(treat))
  # weights <- weights * s.weights
  n <- matrix(0, ncol=2, nrow=2, dimnames = list(c("All", "Weighted"),
                                                 c("Control", "Treated")))

  #                      Control                Treated
  n["All",] <-     c(ESS(s.weights[treat==0]), ESS(s.weights[treat==1]))
  n["Weighted",] <- c(ESS(weights[treat==0]),  ESS(weights[treat==1]))

  return(n)
}
