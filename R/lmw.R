#' Compute linear regression-implied weights
#'
#' @description
#' Computes the weights implied by a linear outcome regression model that would
#' estimate a weighted difference in outcome means equal to the
#' covariate-adjusted treatment effect resulting from the supplied regression
#' model.
#'
#' @details
#' `formula` is interpreted differently depending on whether `method`
#' is `"URI"` or `"MRI"`. When `method = "URI"`, the formula is
#' taken literally as the right-hand side of the outcome model formula. The
#' only difference is that the covariates will be centered based on the
#' argument to `estimand` (see below). When `method = "MRI"`, all
#' references to the treatment are removed (i.e., covariate interactions with
#' treatment become covariate main effects if not already present), and the new
#' formula is taken as the right-hand side of the model formula fit within each
#' treatment group. This is equivalent to allowing all covariates to have both
#' main effects and interactions with treatment after centering the covariates
#' based on the argument to `estimand`. Allowing the treatment to interact
#' with all covariates with `method = "URI"` is equivalent to specifying
#' `method = "MRI"`, and, for binary treatments, the returned weights will
#' be the same when `fixef = NULL`.
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
#' prior to involvement in interactions.
#'
#' \subsection{Estimating a CATE}{
#'
#' When `estimand = "CATE"`, `target` can be supplied either as a
#' single target profile (i.e., a list or a data frame with one row) or as a
#' target dataset, potentially with its own sampling weights, which are
#' supplied to `target.weights`. The variables included in `target`
#' must correspond to all the named *covariates* in `formula`; for
#' example, if `formula = ~ X1 + log(X1) + X2 + X1:X2`, values in
#' `target` must be given for `X1` and `X2`, but not
#' `log(X1)` or `X1:X2`. To choose a target profile value for a
#' factor corresponding to a proportion (e.g., a target value of .5 for a
#' variable like `sex` indicating a target population with a 50-50 sex
#' split), the factor variable must be split into a numeric variable
#' beforehand, e.g., using [model.matrix()] or
#' `cobalt::splitfactor()`. `target` values cannot be given to
#' variables specified using `$`, `[[]]`, or `[]` (e.g.,
#' `data$X1`), so an error will be thrown if they are used in
#' `formula`. When a target dataset is supplied, covariates will be
#' centered at their means in the (`target.weights`-weighted) target
#' dataset.
#' }
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
#' to offer a second level of balancing beyond the implied regression weights;
#' they can be incorporated into the effect estimate either using weighted
#' least squares or using the augmented inverse probability weighting (AIPW)
#' estimator. Base weights do not change the target population, so when target
#' balance is assessed, it is assessed using the unweighted population as the
#' target population.
#'
#' Some forms of weights both change the target population and provide an extra
#' layer of balancing, like propensity score weights that target estimands
#' other than the ATT, ATC, or ATE (e.g., overlap weights), or matching weights
#' where the target population is defined by the matching (e.g., matching with
#' a caliper, cardinality matching, or coarsened exact matching). Because these
#' weights change the target population, they should be supplied to
#' `s.weights` to ensure covariates are appropriately centered. When there
#' are no treatment-by-covariate interactions and `method = "URI"`,
#' whether weights are supplied to `base.weights` or `s.weights` will
#' not matter for the estimation of the weights but will affect the target
#' population in [balance assessment][summary.lmw].
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
#' \subsection{Regression after using \pkg{MatchIt} or \pkg{WeightIt}}{
#' Regression weights can be computed in a matched or weighted sample by
#' supplying a `matchit` or `weightit` object (from \pkg{MatchIt} or
#' \pkg{WeightIt}, respectively) to the `obj` argument of `lmw()`.
#' The estimand, base weights, and sampling weights (if any) will be taken from
#' the supplied object and used in the calculation of the implied regression
#' weights, unless these have been supplied separately to `lmw()`. The
#' `weights` component of the supplied object containing the matching or
#' balancing weights will be passed to `base.weights` and the
#' `s.weights` component will be passed to `s.weights`. Arguments
#' supplied to `lmw()` will take precedence over the corresponding
#' components in the `obj` object.
#' }
#'
#' \subsection{Multi-category treatments}{
#' There are a few differences when the
#' treatment has multiple (i.e., more than 2) categories. If `estimand` is
#' `"ATT"` or `"ATC"`, an argument should be supplied to `focal`
#' identifying which group is the treated or control (i.e., "focal") group,
#' respectively.
#'
#' The key difference, though, is when `method = "URI"`, because in this
#' case the contrast between each pair of treatment groups has its own weights
#' and its own implied target population. Because `lmw()` only produces
#' one set of weights, an argument must be supplied to `contrast`
#' identifying which groups are to be used as the contrast for computing the
#' weights. In addition, to compute the treatment effect corresponding to the
#' chosen contrast as a weighted difference in outcome means, the difference
#' must be taken between the weighted mean of the non-reference group and the
#' weighted mean of *all other groups combined*, rather than simply the
#' weighted mean of the reference group.
#'
#' The implication of this is that contrast statistics computed in the weighted
#' sample involve all units, even those not in the contrasted groups, whereas
#' statistics computed in the unweighted sample only involve units in the
#' contrasted groups. See [summary.lmw()] for more information on
#' assessing balance using the regression weights for multi-category
#' treatments. Given these complications, it is generally best to use
#' `method = "MRI"` with multi-category treatments.
#' }
#'
#' \subsection{Fixed effects}{
#' A fixed effects variable can be supplied to the
#' `fixef` argument. This is equivalent to adding the fixed effects
#' variable as a predictor that does not interact with the treatment or any
#' other covariate. The difference is that computation is much faster when the
#' fixed effect has many levels because demeaning is used rather than including
#' the fixed effect variable as a collection of dummy variables. When using
#' URI, the weights will be the same regardless of whether the fixed effect
#' variable is included as a covariate or supplied to `fixef`; when using
#' MRI, results will differ because the fixed effect variable does not interact
#' with treatment. The fixed effects variable will not appear in the
#' [summary.lmw()] output (but can be added using `addlvariables`
#' argument) or in the model output of [lmw_est()] or
#' [summary.lmw_est()]. Because it does not interact with the
#' treatment, the distribution of the fixed effect variable may not correspond
#' to the target population, so caution should be used if it is expected the
#' treatment effect varies across levels of this variable (in which case it
#' should be included as a predictor). Currently only one fixed effect variable
#' is allowed.
#' }
#'
#' @param formula a one-sided formula with the treatment and covariates on the
#' right-hand side corresponding to the outcome regression model to be fit. The
#' outcome variable is not involved in computing the weights and does not need
#' to be specified. See Details for how this formula is interpreted in light of
#' other options.
#' @param data a data frame containing the variables named in `formula`
#' and `treat`.
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
#' where the model is fit separately in the treatment groups. This affects the
#' interpretation of `formula`. See Details.
#' @param treat the name of the treatment variable in `data`. If
#' unspecified, the first variable present in `formula` will be taken as
#' the treatment variable with a message. See Details.
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
#' @param dr.method the method used to incorporate the `base.weights` into
#' a doubly-robust estimator. Can be one of `"WLS"` for weighted least
#' squares or `"AIPW"` for augmented inverse probability weighting.
#' Ignored when `base.weights` is `NULL`.
#' @param obj a `matchit` or `weightit` object corresponding to the
#' matched or weighted sample in which the implied outcome regression would
#' take place. See Details.
#' @param fixef optional; a string or one-sided formula containing the name of
#' the fixed effects variable in `data`. See Details. Cannot be used with
#' `dr.method = "AIPW"`.
#' @param target a list or data frame containing the target values for each
#' covariate included in `formula`. Ignored with a warning when
#' `estimand` is not `"CATE"`. See Details.
#' @param target.weights a vector of sampling weights to be applied to
#' `target` when supplied as a data frame. Ignored with a warning when
#' `estimand` is not `"CATE"`. See Details.
#' @param contrast for multi-category treatments with `method = "URI"`, a
#' vector containing the names or indices of the two treatment levels to be
#' contrasted (since in this case the weights depend on the specific contrast).
#' See Details.
#' @param focal the level of the treatment variable to be considered "focal"
#' (i.e., the "treated" level when `estimand = "ATT"` or the control level
#' when `estimand = "ATC"`). Ignored when `estimand` is `"ATE"`
#' or `"CATE"`. Otherwise, if unspecified, the second value of
#' `contrast` will be considered focal when `estimand = "ATT"` and
#' the first value of `contrast` will be considered focal when
#' `estimand = "ATC"`. For binary treatments, this generally does not need
#' to be supplied. See Details.
#'
#' @return An `lmw` object, which contains the following components:
#' \item{treat}{the treatment variable, given as a factor.}
#' \item{weights}{the
#' computed implied regression weights.}
#' \item{covs}{a data frame containing
#' the covariates included the model formula.}
#' \item{estimand}{the requested
#' estimand.}
#' \item{method}{the method used to estimate the weights
#' (`"URI"` or `"MRI"`).}
#' \item{base.weights}{the weights supplied to
#' `base.weights`.}
#' \item{s.weights}{the weights supplied to
#' `s.weights`.}
#' \item{dr.method}{when `base.weights` are supplied,
#' the method for computing the doubly-robust weights.}
#' \item{call}{the
#' original call to `lmw()`.}
#' \item{fixef}{the fixed effects variable if
#' supplied to `fixef`.}
#' \item{formula}{the model formula.}
#' \item{target}{the supplied target profile or dataset when `estimand =
#' "CATE"`, after some initial processing. The `"target.weights"`
#' attribute contains the `target.weights` if supplied.}
#' \item{contrast}{the contrasted treatment groups.}
#' \item{focal}{the focal
#' treatment level when `estimand` is `"ATT"` or `"ATC"`.}
#'
#' @seealso [summary.lmw()] for summarizing balance and
#' representativeness; [plot.lmw()] for plotting features of the
#' weights; [lmw_est()] for estimating treatment effects from
#' `lmw` objects; [influence.lmw()] for influence measures;
#' [lm()] for fitting standard regression models.
#'
#' @references Chattopadhyay, A., & Zubizarreta, J. R. (2022). On the implied
#' weights of linear regression for causal inference. *Biometrika*, asac058. \doi{10.1093/biomet/asac058}
#'
#' @examples
#' data("lalonde")
#'
#' # URI regression for ATT
#' lmw.out1 <- lmw(~ treat + age + education + race + married +
#'                   nodegree + re74 + re75, data = lalonde,
#'                 estimand = "ATT", method = "URI",
#'                 treat = "treat")
#' lmw.out1
#' summary(lmw.out1)
#'
#' # MRI regression for ATT
#' lmw.out2 <- lmw(~ treat + age + education + race + married +
#'                   nodegree + re74 + re75, data = lalonde,
#'                 estimand = "ATT", method = "MRI",
#'                 treat = "treat")
#' lmw.out2
#' summary(lmw.out2)
#'
#' @examplesIf requireNamespace("MatchIt", quietly = TRUE)
#' # MRI regression for ATT after propensity score matching
#' m.out <- MatchIt::matchit(treat ~ age + education + race +
#'                             married + nodegree + re74 + re75,
#'                           data = lalonde, method = "nearest",
#'                           estimand = "ATT")
#' lmw.out3 <- lmw(~ treat + age + education + race + married +
#'                   nodegree + re74 + re75, data = lalonde,
#'                 method = "MRI", treat = "treat", obj = m.out)
#' lmw.out3
#' summary(lmw.out3)
#'
#' @examples
#' # MRI regression for CATE with given target profile
#' target.prof <- list(age = 25, education = 11, race = "black",
#'                     married = 0, nodegree = 1, re74 = 0,
#'                     re75 = 0)
#' lmw.out4 <- lmw(~ treat + age + education + race + married +
#'                   nodegree + re74 + re75, data = lalonde,
#'                 estimand = "CATE", method = "MRI",
#'                 treat = "treat", target = target.prof)
#' lmw.out4
#' summary(lmw.out4)
#'
#' # MRI regression for CATE with given target dataset (in
#' # this case, will give the same as with estimand = "ATT")
#' target.data <- subset(lalonde, treat == 1)
#' lmw.out4 <- lmw(~ treat + age + education + race + married +
#'                   nodegree + re74 + re75, data = lalonde,
#'                 estimand = "CATE", method = "MRI",
#'                 treat = "treat", target = target.data)
#' lmw.out4
#' summary(lmw.out4)
#'
#' # URI regression with fixed effects for 'race'
#' lmw.out5 <- lmw(~ treat + age + education + married +
#'                   nodegree + re74 + re75, data = lalonde,
#'                 method = "URI", treat = "treat",
#'                 fixef = ~race)
#' lmw.out5
#'
#' # Produces the same weights as when included as a covariate
#' all.equal(lmw.out1$weights, lmw.out5$weights)
#'
#' # MRI for a multi-category treatment, ATT with 1 as the focal
#' # group
#' lmw.out6 <- lmw(~ treat_multi + age + education + race + married +
#'                   nodegree + re74 + re75, data = lalonde,
#'                 estimand = "ATT", method = "MRI",
#'                 treat = "treat_multi", focal = "1")
#' lmw.out6
#' summary(lmw.out6)
#'
#' # URI for a multi-category treatment; need to specify
#' # contrast because only two groups can be compared at
#' # a time
#' lmw.out7 <- lmw(~ treat_multi + age + education + race + married +
#'                   nodegree + re74 + re75, data = lalonde,
#'                 estimand = "ATE", method = "URI",
#'                 treat = "treat_multi", contrast = c("2", "3"))
#' lmw.out7
#' summary(lmw.out7)
#'

#' @export lmw
lmw <- function(formula, data = NULL, estimand = "ATE", method = "URI", treat = NULL, base.weights = NULL,
                s.weights = NULL, dr.method = "WLS", obj = NULL, fixef = NULL, target = NULL, target.weights = NULL,
                contrast = NULL, focal = NULL) {
  call <- match.call()

  method <- match_arg(method, c("URI", "MRI"))

  estimand <- process_estimand(estimand, target, obj)

  data <- process_data(data, obj)

  base.weights <- do.call("process_base.weights", list(substitute(base.weights), data, obj))

  s.weights <- do.call("process_s.weights", list(substitute(s.weights), data, obj))

  dr.method <- process_dr.method(dr.method, base.weights, method, estimand)

  treat_name <- process_treat_name(treat, formula, data, method, obj)

  fixef <- process_fixef(fixef, formula, data, treat_name)

  #treat changes meaning from treatment name to treatment vector
  treat <- process_treat(treat_name, data)

  check_lengths(treat, data, s.weights, base.weights, fixef)

  contrast <- process_contrast(contrast, treat, method)

  #treat_contrast has levels re-ordered so contrast is first
  treat_contrast <- apply_contrast_to_treat(treat, contrast)

  focal <- process_focal(focal, treat_contrast, estimand)

  X_obj <- get_X_from_formula(formula, data, treat_contrast, method, estimand,
                              target, s.weights, target.weights, focal)

  weights <- get_w_from_X(X_obj$X, treat_contrast, method, base.weights, s.weights, dr.method, fixef)

  out <- list(treat = treat,
              weights = weights,
              covs = X_obj$mf,
              estimand = estimand,
              method = method,
              base.weights = base.weights,
              s.weights = s.weights,
              dr.method = dr.method,
              call = call,
              fixef = fixef,
              formula = formula,
              target = attr(X_obj$target, "target_original"),
              contrast = contrast,
              focal = focal)

  class(out) <- c("lmw_aipw"[!is.null(dr.method) && dr.method == "AIPW"],
                  "lmw_multi"[nlevels(treat) > 2], "lmw")
  return(out)
}

#' @exportS3Method print lmw
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
  if (!is.null(x$dr.method)) {
    cat(sprintf(" - doubly-robust method: %s\n", switch(x$dr.method,
                                                        "WLS" = "weighted least squares (WLS)",
                                                        "AIPW" = "augmented inverse probability weighting (AIPW)")))
  }
  if (inherits(x, "lmw_multi") && !is.null(x$focal)) {
    cat(sprintf(" - target estimand: %s (focal = \"%s\")\n",
                x$estimand, x$focal))
  }
  else {
    cat(sprintf(" - target estimand: %s\n", x$estimand))
  }
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
