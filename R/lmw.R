#Compute weights from formula
lmw <- function(formula, data = NULL, method = "URI", estimand = "ATE", treat = NULL, target = NULL, base.weights = NULL,
                s.weights = NULL, dr.method = "WLS", obj = NULL, contrast = NULL, focal = NULL) {
  call <- match.call()

  method <- match_arg(method, c("URI", "MRI"))

  estimand <- process_estimand(estimand, target, obj)

  data <- process_data(data, obj)

  base.weights <- process_base.weights(base.weights, obj)

  s.weights <- process_s.weights(s.weights, obj)

  dr.method <- process_dr.method(dr.method, base.weights, method)

  treat_name <- process_treat_name(treat, formula, method, obj)

  #treat changes meaning from treatment name to treatment vector
  treat <- process_treat(treat_name, data)

  contrast <- process_contrast(contrast, treat, method)

  #treat_contrast has levels re-ordered so contrast is first
  treat_contrast <- apply_contrast_to_treat(treat, contrast)

  focal <- process_focal(focal, treat_contrast, estimand)

  X_obj <- get_X_from_formula(formula, data, treat_contrast, method, estimand, target, s.weights, focal)

  weights <- get_w_from_X(X_obj$X, treat_contrast, method, base.weights, s.weights, dr.method)

  out <- list(treat = treat,
              weights = weights,
              covs = X_obj$mf,
              estimand = estimand,
              method = method,
              base.weights = base.weights,
              s.weights = s.weights,
              dr.method = dr.method,
              call = call,
              formula = formula,
              target = attr(X_obj$target, "target_original"),
              contrast = contrast,
              focal = focal)

  class(out) <- c("lmw.multi"[nlevels(treat) > 2], "lmw")
  return(out)
}

print.lmw <- function(x, ...) {
  cat("An lmw object\n")
  cat(sprintf(" - treatment: %s (%s levels)\n", attr(x$treat, "treat_name"), nlevels(x$treat)))
  cat(sprintf(" - method: %s\n", switch(x$method, "URI" = "URI (uni-regression imputation)", "MRI" = " MRI (multi-regression imputation)")))
  cat(sprintf(" - number of obs.: %s\n", length(x$treat)))
  cat(sprintf(" - sampling weights: %s\n", if (is.null(x$s.weights)) "none" else "present"))
  cat(sprintf(" - base weights: %s\n",
              if (is.null(x$base.weights)) "none"
              else if (is.null(attr(x$base.weights, "origin"))) "present"
              else switch(attr(x$base.weights, "origin"),
                          "matchit"  = "matching weights from MatchIt",
                          "weightit" = "weights from WeightIt",
                          "present")))
  cat(sprintf(" - target estimand: %s\n", x$estimand))
  if (!is.null(x$covs)) {
    cat(sprintf(" - covariates: %s\n",
        if (length(names(x$covs)) > 30) paste(c(names(x$covs)[1:30], sprintf("and %s more", ncol(x$covs) - 30)), collapse = ", ")
        else paste(names(x$covs), collapse = ", ")))
  }
}
