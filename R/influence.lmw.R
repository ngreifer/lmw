influence.lmw <- function(model, outcome, data = NULL, ...) {

  data <- get_data(data, model)

  outcome <- do.call("get_outcome", list(substitute(outcome), data, model$formula))
  outcome_name <- attr(outcome, "outcome_name")

  attributes(outcome) <- NULL

  #Get model matrix
  obj <- get_X_from_formula(model$formula, data = data, treat = model$treat,
                            method = model$method, estimand = model$estimand, target = model$target,
                            s.weights = model$s.weights, target.weights = attr(model$target, "target.weights"),
                            focal = model$focal)

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

  wr <- fit$residuals * w

  nm <- names(fit$residuals)

  lmw_hat <- hat_fast(obj$X, w)

  N <- length(pos_w)

  SIC <- (N - 1) * abs(wr * w/(1 - lmw_hat))
  #Note: URI formula used; same results with MRI formula when used with MRI weights

  return(list(hat = setNames(lmw_hat, nm),
              wt.res = setNames(wr, nm),
              sic = setNames(SIC/max(SIC), nm)))
}
