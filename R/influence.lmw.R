influence.lmw <- function(model, outcome, data = NULL, ...) {

  data <- get_data(data, model)

  outcome <- do.call("get_outcome", list(substitute(outcome), data, model$formula))
  attributes(outcome) <- NULL

  #Get model matrix
  obj <- get_X_from_formula(model$formula, data = data, treat = model$treat,
                            type = model$type, estimand = model$estimand, target = model$target,
                            s.weights = model$s.weights, focal = model$focal)

  N <- length(model$weights)

  w <- NULL
  if (is.null(model$s.weights) && is.null(model$base.weights)) {
    fit <- lm.fit(x = obj$X, y = outcome)
    nm <- names(fit$residuals)
    r <- fit$residuals
  }
  else {
    if (is.null(model$s.weights)) model$s.weights <- rep(1, length(outcome))
    if (is.null(model$base.weights)) model$base.weights <- rep(1, length(outcome))
    w <- model$s.weights * model$base.weights

    fit <- lm.wfit(x = obj$X, y = outcome, w = w)
    nm <- names(fit$residuals)
    r <- fit$residuals * fit$weights
  }

  lmw_hat <- hat_fast(obj$X, w)

  SIC <- (N - 1) * abs(r * model$weights/(1 - lmw_hat))
  #Note: URI formula used; same results with MRI formula when used with MRI weights

  return(list(hat = setNames(lmw_hat, nm),
              wt.res = setNames(r, nm),
              sic = setNames(SIC/max(SIC), nm)))
}

influence.lmw_est <- function(model, ...) {
stop("NOT READY YET")
  X <- model.matrix(model)
  N <- nrow(X)
  nm <- names(model$residuals)

  w <- model$weights
  if (is.null(w)) {
    r <- model$residuals
  }
  else {
    r <- model$residuals * w
  }

  lmw_hat <- hatvalues(model)

  lmw_weights <- get_w_from_X(X, w)

  SIC <- (N - 1) * abs(r * lmw_weights/(1 - lmw_hat))
  #Note: URI formula used; same results with MRI formula when used with MRI weights

  return(list(hat = setNames(lmw_hat, nm),
              wt.res = setNames(r, nm),
              sic = setNames(SIC/max(SIC), nm)))

}
