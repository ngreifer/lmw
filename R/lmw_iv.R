#' Compute instrumental variable regression-implied weights
#'
#' @description
#' Computes the weights implied by an instrumental variable (IV) model that
#' would estimate a weighted difference in outcome means equal to the treatment
#' effect resulting from the supplied model fit with two-stage least squares.
#'
#' @details
#' `lmw_iv()` computes weights that make the weighted difference in
#' outcome means between the treatment groups equal to the two-stage least
#' squares (2SLS) estimate of the treatment effect. `formula` corresponds
#' to the second-stage (reduced form) model, with the treatment replaced by its
#' fitted values resulting from the first stage model. The first stage is fit
#' by replacing the treatment in the supplied `formula` with the IVs named
#' in `iv` and using the treatment as the outcome. The treatment is
#' assumed to be endogenous and the supplied instrumental variables assumed to
#' be instruments conditional on the other covariates, which are assumed to to
#' be exogenous.
#'
#' When any treatment-by-covariate interactions are present in `formula`
#' or when `method = "MRI"`, covariates are centered at specific values to
#' ensure the resulting weights correspond to the desired estimand as supplied
#' to the `estimand` argument. For the ATE, the covariates are centered at
#' their means in the full sample. For the ATT and ATC, the covariates are
#' centered at their means in the treatment or control group (i.e., the
#' `focal` group), respectively. For the CATE, the covariates are centered
#' according to the argument supplied to `target` (see below). Note that
#' when covariate-by-covariate interactions are present, they will be centered
#' after computing the interaction rather than the interaction being computed
#' on the centered covariates unless `estimand = "CATE"`, in which case
#' the covariates will be centered at the values specified in `target`
#' prior to involvement in interactions. Note that the resulting effect estimate
#' does not actually correspond to the estimand supplied unless all effect
#' heterogeneity is due to the included covariates.
#'
#' When treatment-by-covariate interactions are included in `formula`,
#' additional instruments will be formed as the product of the supplied IVs and
#' the interacting covariates. When `method = "MRI"`, instruments will be
#' formed as the product of the supplied IVs and each of the covariates. All
#' treatment-by-covariate interactions are considered endogenous.
#'
#' \subsection{Base weights and sampling weights}{
#'
#' Base weights (`base.weights`) and sampling weights (`s.weights`)
#' are similar in that they both involve combining weights with an outcome
#' regression model. However, they differ in a few ways. Sampling weights are
#' primarily used to adjust the target population; when the outcome model is
#' fit, it is fit using weighted least squares, and when target balance is
#' assessed, it is assessed using the sampling weighted population as the
#' target population. Centering of covariates in the outcome model is done
#' using the sampling weighted covariate means. Base weights are primarily used
#' to offer a second level of balancing beyond the implied regression weights,
#' i.e., to fit the 2SLS models in the base-weighted sample. Base weights do
#' not change the target population, so when target balance is assessed, it is
#' assessed using the unweighted population as the target population.
#'
#' Some forms of weights both change the target population and provide an extra
#' layer of balancing, like propensity score weights that target estimands
#' other than the ATT, ATC, or ATE (e.g., overlap weights), or matching weights
#' where the target population is defined by the matching (e.g., matching with
#' a caliper, cardinality matching, or coarsened exact matching). Because these
#' weights change the target population, they should be supplied to
#' `s.weights` to ensure covariates are appropriately centered. In
#' `lmw_iv()`, whether weights are supplied to `base.weights` or
#' `s.weights` will not matter for the estimation of the weights but will
#' affect the target population in [balance assessment][summary.lmw].
#'
#' When both `base.weights` and `s.weights` are supplied, e.g., when
#' the base weights are the result of a propensity score model fit with
#' sampling weights, it is assumed the base weights do not incorporate the
#' sampling weights; that is, it is assumed that to estimate a treatment effect
#' *without* regression adjustment, the base weights and the sampling
#' weights would have to be multiplied together. This is true, for example, for
#' the weights in a `matchit` or `weightit` object (see below) but
#' not for weights in the output of `MatchIt::match.data()` unless called
#' with `include.s.weights = FALSE` or weights resulting from
#' `CBPS::CBPS()`.
#' }
#'
#' \subsection{2SLS after using \pkg{MatchIt} or \pkg{WeightIt}}{
#' Instrumental variable regression weights can be computed in a matched or weighted sample
#' by supplying a `matchit` or `weightit` object (from \pkg{MatchIt}
#' or \pkg{WeightIt}, respectively) to the `obj` argument of `lmw()`.
#' The estimand, base weights, and sampling weights (if any) will be taken from
#' the supplied object and used in the calculation of the implied regression
#' weights, unless these have been supplied separately to `lmw_iv()`. The
#' `weights` component of the supplied object containing the matching or
#' balancing weights will be passed to `base.weights` and the
#' `s.weights` component will be passed to `s.weights`. Arguments
#' supplied to `lmw_iv()` will take precedence over the corresponding
#' components in the `obj` object.
#' }
#'
#' \subsection{Multi-category treatments}{
#' Multi-category treatments are not
#' currently supported for `lmw_iv()`.
#' }
#'
#' \subsection{Fixed effects}{
#' A fixed effects variable can be supplied to the
#' `fixef` argument. This is equivalent to adding the fixed effects
#' variable as an exogenous predictor that does not interact with the
#' treatment, IV, or any other covariate. The difference is that computation is
#' much faster when the fixed effect has many levels because demeaning is used
#' rather than including the fixed effect variable as a collection of dummy
#' variables. When using URI, the weights will be the same regardless of
#' whether the fixed effect variable is included as a covariate or supplied to
#' `fixef`; when using MRI, results will differ because the fixed effect
#' variable does not interact with treatment. The fixed effects variable will
#' not appear in the [summary.lmw()] output (but can be added using
#' `addlvariables` argument) or in the model output of [lmw_est()] or
#' [summary.lmw_est()]. Because it does not interact with the
#' treatment, the distribution of the fixed effect variable may not correspond
#' to the target population, so caution should be used if it is expected the
#' treatment effect varies across levels of this variable (in which case it
#' should be included as a predictor). Currently only one fixed effect variable
#' is allowed.
#' }
#'
#' @param formula a one-sided formula with the treatment and covariates on the
#' right-hand side corresponding to the second-stage (reduced form) outcome
#' regression model to be fit. If an outcome variable is supplied on the
#' left-hand side, it will be ignored. This model should not include an IV. See
#' Details for how this formula is interpreted in light of other options.
#' @param data a data frame containing the variables named in `formula`,
#' `treat`, and `iv`.
#' @param estimand the estimand of interest, which determines how covariates
#' are centered. Should be one of `"ATE"` for the average treatment
#' effect, `"ATT"` for the average treatment effect in the treated,
#' `"ATC"` for the average treatment effect in the control, or
#' `"CATE"` for the conditional average treatment effect. When
#' `estimand = "CATE"`, an argument to `target` must be supplied.
#' This argument also affects what [summary.lmw()] considers to be
#' the target population. Default is `"ATE"` unless `obj` is
#' specified, in which case it takes its value from the supplied object.
#' @param method the method used to estimate the weights; either `"URI"`
#' (the default) for uni-regression imputation weights, where a single model is
#' fit to the whole dataset, or `"MRI"` for multi-regression imputation,
#' where the covariates fully interact with the treatment. This affects the
#' interpretation of `formula`. See Details.
#' @param treat the name of the treatment variable in `data`. If
#' unspecified, the first variable present in `formula` will be taken as
#' the treatment variable with a message. Currently, only binary treatments are
#' supported. See Details.
#' @param iv a character vector or one-sided formula containing the names of
#' the IVs in `data`. These variables should not appear in `formula`.
#' Multiple IVs are allowed. See Details. This argument is required.
#' @param base.weights a vector of base weights. See Details. If omitted and
#' `obj` is specified, the weights from the supplied object will be used.
#' Can be supplied as a numeric vector, a string containing the name of the
#' variable in `data` containing the base weights, or the unquoted name of
#' the variable in `data` containing the base weights.
#' @param s.weights a vector of sampling weights. See Details. If omitted and
#' `obj` is specified, the sampling weights from the supplied object will
#' be used. Can be supplied as a numeric vector, a string containing the name
#' of the variable in `data` containing the sampling weights, or the
#' unquoted name of the variable in `data` containing the sampling
#' weights.
#' @param obj a `matchit` or `weightit` object corresponding to the
#' matched or weighted sample in which the implied IV regression would take
#' place. See Details.
#' @param fixef optional; a string or one-sided formula containing the name of
#' the fixed effects variable in `data`. See Details.
#' @param target a list or data frame containing the target values for each
#' covariate included in `formula`. Ignored with a warning when
#' `estimand` is not `"CATE"`.
#' @param target.weights a vector of sampling weights to be applied to
#' `target` when supplied as a data frame. Ignored with a warning when
#' `estimand` is not `"CATE"`.
#' @param contrast ignored.
#' @param focal the level of the treatment variable to be considered "focal"
#' (i.e., the "treated" level when `estimand = "ATT"` or the control level
#' when `estimand = "ATC"`). Ignored when `estimand` is `"ATE"`
#' or `"CATE"`. For binary treatments, this generally does not need to be
#' supplied.
#'
#' @return An `lmw_iv` object, which inherits from `lmw` objects and
#' contains the following components:
#' \item{treat}{the treatment variable,
#' given as a factor.}
#' \item{weights}{the computed implied regression weights.}
#' \item{covs}{a data frame containing the covariates included the model
#' formula.}
#' \item{estimand}{the requested estimand.}
#' \item{method}{the method
#' used to estimate the weights (`"URI"` or `"MRI"`).}
#' \item{base.weights}{the weights supplied to `base.weights`.}
#' \item{s.weights}{the weights supplied to `s.weights`.}
#' \item{call}{the
#' original call to `lmw_iv()`.}
#' \item{fixef}{the fixed effects variable
#' if supplied to `fixef`.}
#' \item{formula}{the model formula.}
#' \item{iv}{the instrumental variables, given as a one-sided formula.}
#' \item{target}{the supplied covariate target values when `estimand =
#' "CATE"`, after some initial processing.}
#' \item{contrast}{the contrasted
#' treatment groups.}
#' \item{focal}{the focal treatment levels when
#' `estimand` is `"ATT"` or `"ATC"`.}
#'
#' All functions that lack a specific `lmw_iv` method work with
#' `lmw_iv` objects as they do for `lmw` objects, such as
#' [summary.lmw()], [plot.lmw()], etc.
#'
#' @seealso [summary.lmw()] for summarizing balance and
#' representativeness; [plot.lmw()] for plotting features of the
#' weights; [lmw_est()] for estimating treatment effects from
#' `lmw_iv` objects; [influence.lmw()] for influence measures;
#' `ivreg()` in the \pkg{ivreg} package for fitting 2SLS models.
#'
#' @references Chattopadhyay, A., & Zubizarreta, J. R. (2023). On the implied weights of linear regression for causal inference. *Biometrika*, 110(3), 615–629. \doi{doi.org/10.1093/biomet/asac058}
#'
#' @examples
#' # URI for the ATT using instrument `Ins`
#' lmw.out <- lmw_iv(~ treat + age + education + race +
#'                     re74, data = lalonde,
#'                   estimand = "ATT", method = "URI",
#'                   treat = "treat", iv = ~Ins)
#' lmw.out
#' summary(lmw.out, addlvariables = ~married + re75)

#' @export lmw_iv
lmw_iv <- function(formula, data = NULL, estimand = "ATE", method = "URI", treat = NULL, iv, base.weights = NULL,
                s.weights = NULL, obj = NULL, fixef = NULL, target = NULL, target.weights = NULL, contrast = NULL, focal = NULL) {
  call <- match.call()

  if (missing(iv)) {
    chk::err("an argument to `iv` specifying the instrumental variable(s) is required")
  }

  chk::chk_string(method)
  method <- toupper(method)
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

  focal <- process_focal(focal, treat_contrast, estimand, obj)

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
  out
}

#' @exportS3Method print lmw_iv
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

  invisible(x)
}
