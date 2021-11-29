#Compute weights from formula
lmw <- function(formula, data = NULL, type = "URI", estimand = "ATE", treat = NULL, target = NULL, base.weights = NULL, s.weights = NULL) {
  call <- match.call()

  type <- match_arg(type, c("URI", "MRI"))
  estimand <- match_arg(estimand, c("ATE", "ATT", "ATC", "CATE"))

  rhs_formula <- delete.response(terms(formula, data = data))

  if (type == "URI") obj <- get_X_from_formula_URI(rhs_formula, data, treat, estimand, target, s.weights)
  else if (type == "MRI") obj <- get_X_from_formula_MRI(rhs_formula, data, treat, estimand, target, s.weights)

  weights <- get_w_from_X(obj$X, base.weights, s.weights)

  out <- list(treat = obj$t,
              weights = weights,
              covs = obj$mf,
              estimand = estimand,
              type = type,
              base.weights = base.weights,
              s.weights = s.weights,
              call = call,
              formula = formula)

  class(out) <- "lmw"
  return(out)
}

print.lmw <- function(x, ...) {
  cat("An lmw object\n")

}
