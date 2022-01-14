get_w_from_X <- function(X, treat, type, base.weights = NULL, s.weights = NULL, dr.method = "WLS") {

  if (is.null(s.weights)) s.weights <- rep(1, nrow(X))
  if (is.null(base.weights) || dr.method == "AIPW") w <- s.weights
  else w <- s.weights*base.weights
  rw <- sqrt(w)
  qr_X <- qr(rw*X)
  p <- qr_X$rank

  XtX1 <- chol2inv(qr_X$qr[1:p, 1:p, drop = FALSE])

  #Remove linearly dependent columns
  X <- X[, qr_X$pivot[1:p], drop = FALSE]

  XXtX1 <- (w*X) %*% XtX1

  if (is.null(s.weights)) s.weights <- rep(1, nrow(X))

  if (is.null(base.weights) || dr.method == "WLS") {
    if (type == "URI") {
      #Treated group dummy always in second column
      weights <- XXtX1[,2]
      t <- X[,2]
      weights[t == 0] <- -weights[t == 0]
    }
    else if (type == "MRI") {
      weights <- rep(0, length(treat))
      for (i in seq_len(nlevels(treat))) {
        weights[treat == levels(treat)[i]] <- XXtX1[treat == levels(treat)[i],i]
      }
    }
  }
  else if (dr.method == "AIPW") {
    if (type == "URI") {
      #Treated group dummy always in second column
      r.weights <- XXtX1[,2]
      t <- X[,2]
      r.weights[t == 0] <- -r.weights[t == 0]
    }
    else if (type == "MRI") {
      r.weights <- rep(0, length(treat))
      for (i in seq_len(nlevels(treat))) {
        r.weights[treat == levels(treat)[i]] <- XXtX1[treat == levels(treat)[i],i]
      }
    }

    for (i in 0:1) {
      base.weights[t == i] <- base.weights[t == i]/sum(base.weights[t == i])
    }

    pred.weights <- base.weights %*% (XXtX1 %*% t(X))
    for (i in 0:1) {
      pred.weights[t == i] <- pred.weights[t == i]/sum(pred.weights[t == i])
    }

    weights <- r.weights + base.weights - pred.weights

  }

  return(drop(weights))
}

get_w_from_X_iv <- function(X, treat, type, base.weights = NULL, s.weights = NULL) {
  #X should be from get_1st_stage_X_from_formula_iv()
  #iv i nsecond column

  #Remove linearly dependent columns
  qr_X <- qr(X)
  X <- X[, qr_X$pivot[seq_len(qr_X$rank)], drop = FALSE]
  iv <- X[,2]
  X <- X[,-2, drop = FALSE]

  if (is.null(s.weights)) s.weights <- rep(1, nrow(X))
  if (is.null(base.weights)) base.weights <- rep(1, nrow(X))

  # wX <- base.weights*s.weights*X
  # Px <- tcrossprod(X %*% chol2inv(chol(crossprod(X, wX))), wX)
  # weights <- base.weights*s.weights*drop((diag(nrow(X)) - Px) %*% iv)

  rw <- sqrt(base.weights*s.weights)
  weights <- .lm.fit(rw*X, rw*iv)$residuals/rw  #faster, uses internal C routine

  tval1 <- levels(treat)[1]
  weights[treat != tval1] <- -weights[treat != tval1]

  for (i in levels(treat)) {
    weights[treat == i] <- weights[treat == i] / sum(weights[treat == i])
  }

  return(weights)
}
