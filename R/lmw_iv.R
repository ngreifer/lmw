#Compute weights from formula
lmw_iv <- function(formula, data = NULL, estimand = "ATE", method = "URI", treat = NULL, iv, base.weights = NULL,
                s.weights = NULL, obj = NULL, fixef = NULL, target = NULL, target.weights = NULL, contrast = NULL, focal = NULL) {
  call <- match.call()

  if (missing(iv)) {
    stop("An argument to 'iv' specifying the instrumental variable(s) is required.", call. = FALSE)
  }

  method <- match_arg(method, c("URI", "MRI"))

  estimand <- process_estimand(estimand, target, obj)

  data <- process_data(data, obj)

  base.weights <- do.call("process_base.weights", list(substitute(base.weights), data, obj))

  s.weights <- do.call("process_s.weights", list(substitute(s.weights), data, obj))

  treat_name <- process_treat_name(treat, formula, data, method, obj)

  fixef <- process_fixef(fixef, formula, data, treat_name)

  #treat changes meaning from treatment name to treatment vector
  treat <- process_treat(treat_name, data, multi.ok = FALSE)

  check_lengths(treat, data, s.weights, base.weights, fixef)

  contrast <- process_contrast(contrast, treat, method)

  #treat_contrast has levels re-ordered so contrast is first
  treat_contrast <- apply_contrast_to_treat(treat, contrast)

  focal <- process_focal(focal, treat_contrast, estimand)

  iv <- process_iv(iv, formula, data)

  X_obj <- get_1st_stage_X_from_formula_iv(formula, data, treat_contrast, iv,
                                           method, estimand, target, s.weights,
                                           target.weights, focal)

  weights <- get_w_from_X_iv(X_obj$X, X_obj$A, treat_contrast, method, base.weights, s.weights,
                             fixef)

  out <- list(treat = treat,
              weights = weights,
              covs = X_obj$mf,
              estimand = estimand,
              method = method,
              base.weights = base.weights,
              s.weights = s.weights,
              call = call,
              fixef = fixef,
              formula = formula,
              iv = iv,
              target = attr(X_obj$target, "target_original"),
              contrast = contrast,
              focal = focal)

  class(out) <- c("lmw_iv", "lmw")
  return(out)
}

print.lmw_iv <- function(x, ...) {
  cat("An lmw_iv object\n")
  cat(sprintf(" - treatment: %s (%s levels)\n", attr(x$treat, "treat_name"), nlevels(x$treat)))
  cat(sprintf(" - instrument%s: %s\n", if (length(attr(terms(x$iv), "term.labels")) > 1) "s" else "",
              paste(attr(terms(x$iv), "term.labels"), collapse = ", ")))
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
  if (!is.null(x$fixef)) {
    cat(sprintf(" - fixed effect: %s\n",
                attr(x$fixef, "fixef_name")))
  }
}
