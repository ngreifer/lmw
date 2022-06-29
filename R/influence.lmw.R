influence.lmw <- function(model, outcome, data = NULL, ...) {

  data <- get_data(data, model)

  #Get model matrix
  obj <- get_X_from_formula(model$formula, data = data, treat = model$treat,
                            method = model$method, estimand = model$estimand, target = model$target,
                            s.weights = model$s.weights, target.weights = attr(model$target, "target.weights"),
                            focal = model$focal)

  outcome <- do.call("get_outcome", list(substitute(outcome), data, model$formula,
                                         obj$X))
  outcome_name <- attr(outcome, "outcome_name")

  attributes(outcome) <- NULL

  #Fit regression model; note use lm.wfit() instead of lm() because
  #we already have the model matrix (obj$X)

  if (is.null(model$s.weights)) model$s.weights <- rep(1, length(outcome))
  if (is.null(model$base.weights)) model$base.weights <- rep(1, length(outcome))
  w <- model$s.weights * model$base.weights
  pos_w <- which(w > 0)

  if (!is.null(model$fixef)) {
    for (i in seq_len(ncol(obj$X))[-1]) {
      obj$X[pos_w,i] <- demean(obj$X[pos_w,i], model$fixef[pos_w], w[pos_w])
    }
    outcome[pos_w] <- demean(outcome[pos_w], model$fixef[pos_w], w[pos_w])
  }

  fit <- lm.wfit(x = obj$X[pos_w,,drop = FALSE], y = outcome[pos_w], w = w[pos_w])

  n <- nrow(obj$X)

  wr <- rep(0, n)
  wr[pos_w] <- fit$residuals * w[pos_w]

  nm <- rownames(obj$X)
  if (is.null(nm)) nm <- seq_len(n)

  lmw_hat <- hat_fast(obj$X, w)

  SIC <- (n - 1) * abs(wr * model$weights/(1 - lmw_hat))
  #Note: URI formula used; same results with MRI formula when used with MRI weights

  return(list(hat = setNames(lmw_hat, nm),
              wt.res = setNames(wr, nm),
              sic = setNames(SIC/max(SIC), nm)))
}

influence.lmw_est <- function(model, ...) {

  n <- sum(model$weights > 0)

  wr <- model$residuals * model$weights

  nm <- rownames(model$model.matrix)
  if (is.null(nm)) nm <- seq_len(n)

  lmw_hat <- hatvalues(model)

  SIC <- (n - 1) * abs(wr * model$lmw.weights/(1 - lmw_hat))
  #Note: URI formula used; same results with MRI formula when used with MRI weights

  return(list(hat = setNames(lmw_hat, nm),
              wt.res = setNames(wr, nm),
              sic = setNames(SIC/max(SIC), nm)))
}

