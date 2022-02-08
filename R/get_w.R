get_w_from_X <- function(X, treat, method, base.weights = NULL, s.weights = NULL, dr.method = "WLS",
                         fixef = NULL) {

  if (is.null(s.weights)) s.weights <- rep(1, nrow(X))
  if (is.null(base.weights) || dr.method == "AIPW") w <- s.weights
  else w <- s.weights*base.weights

  if (method == "URI") {
    t <- X[,2]
  }

  if (!is.null(fixef)) {
    for (i in seq_len(ncol(X))) {
      X[,i] <- demean(X[,i], fixef, w)
    }
  }

  rw <- sqrt(w)
  qr_X <- qr(rw*X)
  p <- qr_X$rank

  XtX1 <- chol2inv(qr_X$qr[1:p, 1:p, drop = FALSE])

  #Remove linearly dependent columns
  X <- X[, qr_X$pivot[1:p], drop = FALSE]

  # XXtX1 <- (w*X) %*% XtX1

  if (is.null(s.weights)) s.weights <- rep(1, nrow(X))

  if (method == "URI") {
    #Treated group dummy always in second column
    weights <- drop(w * (X %*% XtX1[, qr_X$pivot[1:p] == 2, drop = FALSE]))
    weights[t == 0] <- -weights[t == 0]
  }
  else if (method == "MRI") {
    weights <- rep(0, length(treat))
    for (i in seq_len(nlevels(treat))) {
      in_t <- which(treat == levels(treat)[i])
      if (i %in% qr_X$pivot[1:p]) {
        weights[in_t] <- w[in_t] * (X[in_t,,drop = FALSE] %*% XtX1[, qr_X$pivot[1:p] == i, drop = FALSE])
      }
      else {
        weights[in_t] <- -w[in_t] * (X[in_t,,drop = FALSE] %*% XtX1[, 1, drop = FALSE])
      }
    }
  }

  if (!is.null(base.weights) && dr.method == "AIPW") {

    for (i in levels(t)) {
      base.weights[t == i] <- base.weights[t == i]/sum(base.weights[t == i])
    }

    # pred.weights <- base.weights %*% XXtX1 %*% t(X)
    # pred.weights <- base.weights - .lm.fit(rw*X, rw*base.weights)$residuals/rw
    # weights <- r.weights + base.weights - pred.weights

    weights <- weights + .lm.fit(rw*X, rw*base.weights)$residuals/rw
  }

  return(drop(weights))
}

get_w_from_X_iv <- function(X, treat, method, base.weights = NULL, s.weights = NULL) {
  #X should be from get_1st_stage_X_from_formula_iv()$X

  iv_names <- attr(X, "iv_names")

  #Remove linearly dependent columns
  qr_X <- qr(X)
  X <- X[, qr_X$pivot[seq_len(qr_X$rank)], drop = FALSE]

  iv_ind <- which(colnames(X) %in% iv_names)

  if (is.null(s.weights)) s.weights <- rep(1, nrow(X))
  if (is.null(base.weights)) base.weights <- rep(1, nrow(X))

  t <- as.numeric(treat == levels(treat)[1])

  rw <- sqrt(base.weights*s.weights)
  weights <- .lm.fit(rw*X[,-iv_ind, drop = FALSE], rw*t)$residuals -
    .lm.fit(rw*X, rw*t)$residuals

  for (i in levels(treat)) {
    weights[treat == i] <- weights[treat == i] / sum(weights[treat == i])
  }

  return(weights)
}
