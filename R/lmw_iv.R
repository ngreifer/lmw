#Compute weights from formula
lmw_iv <- function(formula, data = NULL, type = "URI", estimand = "ATE", treat = NULL, iv, target = NULL, base.weights = NULL,
                s.weights = NULL, obj = NULL, contrast = NULL, focal = NULL) {
  call <- match.call()

  type <- match_arg(type, c("URI"))

  estimand <- process_estimand(estimand, target, obj)

  data <- process_data(data, obj)

  base.weights <- process_base.weights(base.weights, obj)

  s.weights <- process_s.weights(s.weights, obj)

  treat_name <- process_treat_name(treat, formula, type, obj)

  #treat changes meaning from treatment name to treatment vector
  treat <- process_treat(treat_name, data, multi.ok = FALSE)

  contrast <- process_contrast(contrast, treat, type)

  #treat_contrast has levels re-ordered so contrast is first
  treat_contrast <- apply_contrast_to_treat(treat, contrast)

  focal <- process_focal(focal, treat_contrast, estimand)

  iv_name <- process_iv_name(iv, formula)

  iv <- process_iv(iv_name, data)

  X_obj <- get_1st_stage_X_from_formula_iv(formula, data, treat_contrast, iv,
                                           type, estimand, target, s.weights, focal)

  weights <- get_w_from_X_iv(X_obj$X, treat_contrast, type, base.weights, s.weights)

  out <- list(treat = treat,
              iv = iv,
              weights = weights,
              covs = X_obj$mf,
              estimand = estimand,
              type = type,
              base.weights = base.weights,
              s.weights = s.weights,
              call = call,
              formula = formula,
              target = attr(X_obj$target, "target_original"),
              contrast = contrast,
              focal = focal)

  class(out) <- c("lmw_iv", "lmw")
  return(out)
}

print.lmw_iv <- function(x, ...) {
  cat("An lmw_iv object\n")
  cat(sprintf(" - treatment: %s (%s levels)\n", attr(x$treat, "treat_name"), nlevels(x$treat)))
  cat(sprintf(" - instrument: %s\n", attr(x$iv, "iv_name")))
  cat(sprintf(" - type: %s\n", switch(x$type, "URI" = "URI (uni-regression imputation)", "MRI" = " MRI (multi-regression imputation)")))
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
