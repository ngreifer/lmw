#' Estimate a treatment effect from a linear model
#'
#' @description
#' `lmw_est()` fits the outcome regression corresponding to the model used
#' to compute the weights in the supplied `lmw` object and returns the
#' model coefficients and their covariance matrix. Use
#' [summary.lmw_est()] to compute and view the treatment effect and
#' potential outcome mean estimates and their standard errors.
#'
#' @details
#' `lmw_est()` uses [lm.fit()] or [lm.wfit()] to fit
#' the outcome regression model (and first stage model for `lmw_iv`
#' objects) and returns the output of these functions augmented with other
#' components related to the estimation of the weights. Unlike with
#' `lm.[w]fit()`, the covariance matrix of the parameter estimates is also
#' included in the output.
#'
#' For `lmw` objects, the model fit is that supplied to the `formula`
#' input to `lmw()` except that it is fit in a dataset appropriately
#' centered to ensure the estimand corresponds with the one requested. When
#' `method = "MRI"` in the call to `lmw()`, the model is fit as an
#' interaction between the treatment and all the (centered) terms in the model
#' formula. The results will be similar to those from using [lm()] on
#' this model and supplied data except that the covariates are centered
#' beforehand. The product of the sampling weights and base weights supplied to
#' `lmw()`, if any, will be supplied to `lm.wfit()` to fit the model
#' using weighted least squares.
#'
#' For `lmw_aipw` objects, the model is fit as above except that base
#' weights are not included in the model fitting and are instead used to
#' compute additional augmentation terms that are added to the estimated
#' potential outcome means from the outcome regression. The variance-covariance
#' matrix is computed using M-estimation; this corresponds to the HC0 robust
#' covariance matrix for the model parameters with the base weights treated as
#' fixed, which yields conservative standard errors for the ATE. Inference is
#' only approximate for the ATT and ATC.
#'
#' For `lmw_iv` objects, the first stage model is constructed by removing
#' the treatment from the supplied model formula, adding the instrumental
#' variable as a main effect, and using the treatment variable as the outcome.
#' For the second stage (reduced form) model, the fitted values of the
#' treatment from the first stage model are used in place of the treatment in
#' the outcome model. The results are similar to those from using
#' `ivreg::ivreg()`, and the coefficients estimates will be the same
#' except for the intercept due to the centering of covariates.
#'
#' Although some coefficients in the model may be interpretable as treatment
#' effect estimates, [summary.lmw_est()] should be used to view and
#' extract the treatment effect and potential outcome mean estimates, standard
#' errors, and other model statistics. The output of `lmw_est()` should
#' rarely be used except to be supplied to `summary()`.
#'
#' @param x an `lmw` or `lmw_iv` object; the output of a call to
#' [lmw()] or [lmw_iv()].
#' @param outcome the name of the outcome variable. Can be supplied as a string
#' containing the name of the outcome variable or as the outcome variable
#' itself. If not supplied, the outcome variable in the `formula` supplied
#' to `lmw()` or `lmw_iv()`, if any, will be used.
#' @param data an optional data frame containing the outcome variable named in
#' `outcome` and the cluster variable(s) when `cluster` is supplied
#' as a `formula`.
#' @param robust whether to compute the robust covariance matrix for the model
#' coefficients. Allowable values include those allowed for the `type`
#' argument of [sandwich::vcovHC()] or [sandwich::vcovCL()]
#' when `cluster` is specified. Can also be specified as `TRUE` (the
#' default), which means `"HC3"` or `"HC1"` when `cluster` is
#' specified, or `FALSE`, which means `"const"` (i.e., the standard
#' non-robust covariance). When `cluster` is specified, `robust` will
#' be set to `TRUE` if `FALSE`. When AIPW is used, `robust` is
#' ignored; the HC0 robust covariance matrix is used.
#' @param cluster the clustering variable(s) for computing a cluster-robust
#' covariance matrix. See [sandwich::vcovCL()]. If supplied as a
#' `formula`, the clustering variables must be present in the original
#' dataset used to compute the weights or `data`. When AIPW is used,
#' `cluster` is ignored.
#' @param \dots other arguments passed to [sandwich::vcovHC()] or
#' [sandwich::vcovCL()].
#'
#' @return An `lmw_est` object with the following components:
#' \item{coefficients, residuals, fitted.values, effects, weights, rank,
#' df.residual, qr}{for `lmw` objects, the output of the
#' [lm.fit()] or [lm.wfit()] call used to fit the outcome
#' model. For `lmw_iv` objects, the output of the [lm.fit()] or
#' [lm.wfit()] call used to fit the the second stage model, with
#' `residuals` corresponding to the residuals computed when substituting
#' the true treatment variable in place of the fitted treatment values in the
#' model.}
#' \item{model.matrix}{the model matrix (supplied to the `x`
#' argument of `lm.fit`).}
#' \item{vcov}{the estimated covariance matrix of
#' the parameter estimates as produced by [sandwich::vcovHC()] or
#' [sandwich::vcovCL()].}
#' \item{lmw.weights}{the implied regression
#' weights computed by `lmw_est()`.}
#' \item{call}{the call to
#' `lmw_est()`.}
#' \item{estimand}{the requested estimand.}
#' \item{focal}{the
#' focal treatment level when `estimand` is `"ATT"` or `"ATC"`.}
#' \item{method}{the method used to estimate the weights (`"URI"` or
#' `"MRI"`).}
#' \item{robust}{the type standard error used.}
#' \item{outcome}{the name of the outcome variable.}
#' \item{treat_levels}{the
#' levels of the treatment.}
#'
#' When AIPW is used, the object will be of class `lmw_est_aipw`, which
#' inherits from `lmw_est`, and contains the additional components:
#' \item{coef_aipw}{the model-predicted potential outcome means (`mu`) and
#' the augmentation terms (`aug`).}
#' \item{vcov_aipw}{the covariance matrix
#' of the quantities in `coef_aipw`.}
#'
#' When weights are included in the estimation (i.e., `base.weights` or
#' `s.weights` supplied to `lmw()` or `lmw_iv()`), any units
#' will weights equal to zero will be removed from the data prior to model
#' fitting.
#'
#' Methods exist for `lmw_est` objects for [model.matrix()],
#' [vcov()], [hatvalues()], [sandwich::bread()],
#' and [sandwich::estfun()], all of which are used internally to
#' compute the parameter estimate covariance matrix. The first two simply
#' extract the corresponding component from the `lmw_est` object and the
#' last three imitate the corresponding methods for `lm` objects (or
#' `ivreg` objects for `lmw_iv` inputs). Other regression-related
#' functions, such as [coef()], [residuals()], and
#' [fitted()], use the default methods and should work correctly with
#' `lmw_est` objects.
#'
#' Note that when fixed effects are supplied through the `fixef` argument
#' to `lmw()` or `lmw_iv()`, standard error estimates computed using
#' functions outside \pkg{lmw} may not be accurate due to issues relating to
#' degrees of freedom. In particular, this affects conventional and HC1-robust
#' standard errors. Otherwise, `sandwich::vcovHC()` can be used to compute
#' standard errors (setting `type = "const"` for conventional standard
#' errors), though `sandwich::vcovCL()` may not work as expected and
#' should not be used. To calculate cluster-robust standard errors, supply an
#' argument to `cluster` in `lmw_est()`.
#'
#' @note `lmw_est()` uses non-standard evaluation to interpret its
#' `outcome` argument. For programmers who wish to use `lmw_est()`
#' inside other functions, an effective way to pass the name of an arbitrary
#' outcome (e.g., `y` passed as a string) is to use [do.call()],
#' for example: \preformatted{fun <- function(model, outcome, data) {
#' do.call("lmw_est", list(model, outcome, data)) } } When using
#' `lmw_est()` inside [lapply()] or `purrr::map` to loop
#' over outcomes, this syntax must be used as well.
#'
#' @seealso [summary.lmw_est()] for viewing and extracting the
#' treatment effect and potential outcome mean estimates, standard errors, and
#' other model statistics; [lmw()] or [lmw_iv()] for
#' estimating the weights that correspond to the model estimated by
#' `lmw_est()`; [lm()] and [lm.fit()] for fitting the
#' corresponding model; `ivreg()` in the \pkg{ivreg} package for fitting
#' 2SLS models; [influence.lmw_est()] for influence measures
#'
#' @examples
#' data("lalonde")
#'
#' # MRI regression for ATT
#' lmw.out1 <- lmw(~ treat + age + education + race + married +
#'                   nodegree + re74 + re75, data = lalonde,
#'                   estimand = "ATT", method = "MRI",
#'                   treat = "treat")
#'
#' lmw.fit1 <- lmw_est(lmw.out1, outcome = "re78")
#' lmw.fit1
#'
#' summary(lmw.fit1)
#'
#' @examplesIf requireNamespace("MatchIt", quietly = TRUE)
#' # MRI regression for ATT after propensity score matching
#' m.out <- MatchIt::matchit(treat ~ age + education + race +
#'                             married + nodegree + re74 + re75,
#'                           data = lalonde, method = "nearest",
#'                           estimand = "ATT")
#' lmw.out2 <- lmw(~ treat + age + education + race + married +
#'                   nodegree + re74 + re75, data = lalonde,
#'                 method = "MRI", treat = "treat", obj = m.out)
#'
#' ## Using a cluster-robust SE with subclass (pair membership)
#' ## as the cluster variable
#' lmw.fit2 <- lmw_est(lmw.out2, outcome = "re78", cluster = ~subclass)
#' lmw.fit2
#'
#' summary(lmw.fit2)
#'
#' @examples
#' # AIPW for ATE with MRI regression after propensity score
#' # weighting
#' ps <- glm(treat ~ age + education + race + married + nodegree +
#'             re74 + re75, data = lalonde,
#'             family = binomial)$fitted
#' ipw <- ifelse(lalonde$treat == 1, 1/ps, 1/(1-ps))
#'
#' lmw.out3 <- lmw(re78 ~ treat + age + education + race + married +
#'                   nodegree + re74 + re75, data = lalonde,
#'                 method = "MRI", treat = "treat",
#'                 base.weights = ipw, dr.method = "AIPW")
#' lmw.fit3 <- lmw_est(lmw.out3)
#' lmw.fit3
#'
#' summary(lmw.fit3)
#'
#' # MRI for multi-category treatment ATE
#' lmw.out3 <- lmw(~ treat_multi + age + education + race + married +
#'                   nodegree + re74 + re75, data = lalonde,
#'                 estimand = "ATE", method = "MRI",
#'                 treat = "treat_multi")
#' lmw.fit3 <- lmw_est(lmw.out3, outcome = "re78")
#' lmw.fit3
#'
#' summary(lmw.fit3)

#' @export lmw_est
lmw_est <- function(x, ...) {
  UseMethod("lmw_est")
}

#' @exportS3Method lmw_est lmw
#' @rdname lmw_est
lmw_est.lmw <- function(x, outcome, data = NULL, robust = TRUE, cluster = NULL, ...) {

  call <- match.call()

  data <- get_data(data, x)

  #Get model matrix
  obj <- get_X_from_formula(x$formula, data = data, treat = x$treat,
                            method = x$method, estimand = x$estimand, target = x$target,
                            s.weights = x$s.weights, target.weights = attr(x$target, "target.weights"),
                            focal = x$focal)

  outcome <- do.call("get_outcome", list(substitute(outcome), data, x$formula,
                                         obj$X))
  outcome_name <- attr(outcome, "outcome_name")

  attributes(outcome) <- NULL

  #Fit regression model; note use lm.[w]fit() instead of lm() because
  #we already have the model matrix (obj$X)

  if (is.null(x$s.weights)) x$s.weights <- rep(1, length(outcome))
  if (is.null(x$base.weights)) x$base.weights <- rep(1, length(outcome))
  w <- x$s.weights * x$base.weights
  pos_w <- which(w > 0)

  if (!is.null(x$fixef)) {
    for (i in seq_len(ncol(obj$X))[-1]) {
      obj$X[,i] <- demean(obj$X[,i], x$fixef, w)
    }
    outcome <- demean(outcome, x$fixef, w)
  }

  fit <- lm.wfit(x = obj$X, y = outcome, w = w)

  # non_pos_w <- which(w <= 0)
  # fit$na.action <- setNames(non_pos_w, rownames(obj$X[non_pos_w]))
  # class(fit$na.action) <- "omit"

  class(fit) <- c("lmw_est")

  fit$model.matrix <- obj$X

  if (!is.null(x$fixef)) {
    fit$df.residual <- fit$df.residual - length(unique(x$fixef[pos_w])) + 1
    fit$fixef <- droplevels(x$fixef)
  }

  if (isTRUE(robust)) {
    if (is.null(cluster)) robust <- "HC3"
    else robust <- "HC1"
  }
  else if (isFALSE(robust)) {
    if (is.null(cluster)) robust <- "const"
    else {
      robust <- "HC1"
      chk::wrn("setting `robust = \"HC1\"` because `cluster` is non-`NULL`")
    }
  }
  else if (!is.character(robust) || length(robust) != 1 ||
           !robust %in% eval(formals(sandwich::vcovHC.default)$type)) {
    chk::err("`robust` must be `TRUE`, `FALSE`, or one of the allowable inputs to the `type` argument of `sandwich::vcovHC()`")
  }

  #Subset model outputs to those with positive weights
  #for compatibility with vcovHC
  fit_sub <- subset_fit(fit)

  if (robust == "const") { #Regular OLS vcov; faster than sandwich
    fit$vcov <- {
      if (is.null(w)) bread(fit_sub)/length(residuals(fit_sub)) * sum(residuals(fit_sub)^2)/fit_sub$df.residual
      else            bread(fit_sub)/length(residuals(fit_sub)) * sum(fit_sub$weights * residuals(fit_sub)^2)/fit_sub$df.residual
    }
  }
  else if (is.null(cluster)) {
    fit$vcov <- sandwich::vcovHC(fit_sub, type = robust, ...)
  }
  else {
    if (inherits(cluster, "formula")) {
      cluster <- model.frame(cluster,
                             data = data,
                             na.action = na.pass)
    }
    else {
      cluster <- as.data.frame(cluster)
    }

    if (nrow(cluster) == nrow(data)) cluster <- cluster[pos_w,, drop = FALSE]
    else if (nrow(cluster) != length(pos_w)) {
      chk::err("`cluster` must have the same number of rows as the original data set")
    }

    withCallingHandlers({
      fit$vcov <- sandwich::vcovCL(fit_sub, type = robust, cluster = cluster, ...)
    },
    warning = function(w) {
      if (conditionMessage(w) != "clustered HC2/HC3 are only applicable to (generalized) linear regression models")
        chk::wrn(w, tidy = FALSE)
      invokeRestart("muffleWarning")
    })
  }

  #Need to correct df for HC1 when fixed effects are present
  if (!is.null(x$fixef) && robust %in% c("const", "HC1")) {
    n <- length(pos_w)
    fit$vcov <- fit$vcov * (n - ncol(fit$model.matrix))/fit$df.residual
  }

  #For cluster SE, adjust df as the min number of clusters across clustering variables
  if (!is.null(cluster)) {
    fit$df.residual <- min(vapply(cluster, function(cl) length(unique(cl)), numeric(1L))) - 1
  }

  fit$lmw.weights <- x$weights
  fit$call <- call
  fit$estimand <- x$estimand
  fit$focal <- x$focal
  fit$method <- x$method
  fit$robust <- robust
  fit$outcome <- outcome_name
  fit$treat_levels <- levels(x$treat)

  fit
}

#' @exportS3Method print lmw_est
print.lmw_est <- function(x, ...) {
  cat(sprintf("An %s object\n", class(x)[1]))
  cat(" - outcome:", x$outcome, "\n")
  cat(" - standard errors:", if (utils::hasName(x$call, "cluster") && !is.null(x$call[["cluster"]])) "cluster",
      if (x$robust == "const") "usual" else sprintf("robust (%s)", x$robust), "\n")
  cat(" - estimand:", x$estimand, "\n")
  cat(" - method:", x$method, "\n")
  if (!is.null(x$fixef)) cat(" - fixed effects:", attr(x$fixef, "fixef_name"), "\n")
  cat("\n")
  cat("Use summary() to examine estimates, standard errors, p-values, and confidence intervals.", "\n")
  invisible(x)
}

subset_fit <- function(fit) {
  if (is.null(fit$weights)) return(fit)
  pos_w <- fit$weights > 0
  fit$residuals <- fit$residuals[pos_w]
  fit$fitted.values <- fit$fitted.values[pos_w]
  fit$weights <- fit$weights[pos_w]
  fit$model.matrix <- fit$model.matrix[pos_w,,drop = FALSE]
  if (!is.null(fit$fixef)) fit$fixef <- fit$fixef[pos_w]

  fit
}

#####
# lmw_est.formula <- function(x, data = NULL, estimand = "ATE", method = "URI", treat = NULL, base.weights = NULL,
#                             s.weights = NULL, dr.method = "WLS", obj = NULL, target = NULL, target.weights = NULL,
#                             contrast = NULL, focal = NULL,
#                             outcome, robust = TRUE, cluster = NULL, ...) {
#   call <- match.call()
#
#   formula <- x
#
#   method <- match_arg(method, c("URI", "MRI"))
#
#   estimand <- process_estimand(estimand, target, obj)
#
#   data <- process_data(data, obj)
#
#   base.weights <- process_base.weights(base.weights, obj)
#
#   s.weights <- process_s.weights(s.weights, obj)
#
#   dr.method <- process_dr.method(dr.method, base.weights, method)
#
#   treat_name <- process_treat_name(treat, formula, data, method, obj)
#
#   #treat changes meaning from treatment name to treatment vector
#   treat <- process_treat(treat_name, data)
#
#   contrast <- process_contrast(contrast, treat, method)
#
#   #treat_contrast has levels re-ordered so contrast is first
#   treat_contrast <- apply_contrast_to_treat(treat, contrast)
#
#   focal <- process_focal(focal, treat_contrast, estimand)
#
#   # data <- get_data(data, x)
#
#   outcome <- do.call("get_outcome", list(substitute(outcome), data, formula))
#   outcome_name <- attr(outcome, "outcome_name")
#
#   attributes(outcome) <- NULL
#
#   #Get model matrix; note: use original treat, unlike
#   obj <- get_X_from_formula(formula, data, treat, method, estimand,
#                             target, s.weights, target.weights, focal)
#
#   #Fit regression model; note use lm.[w]fit() instead of lm() because
#   #we already have the model matrix (obj$X)
#   w <- NULL
#   rn <- rownames(obj$X)
#   if (is.null(s.weights) && is.null(base.weights)) {
#     fit <- lm.fit(x = obj$X, y = outcome)
#     pos_w <- seq_along(outcome)
#   }
#   else {
#     if (is.null(s.weights)) s.weights <- rep(1, length(outcome))
#     if (is.null(base.weights)) base.weights <- rep(1, length(outcome))
#     w <- s.weights * base.weights
#     pos_w <- which(w > 0)
#     obj$X <- obj$X[pos_w,, drop = FALSE]
#
#     fit <- lm.wfit(x = obj$X, y = outcome[pos_w], w = w[pos_w])
#
#     non_pos_w <- which(w <= 0)
#     fit$na.action <- setNames(non_pos_w, rn[non_pos_w])
#     class(fit$na.action) <- "omit"
#   }
#
#   class(fit) <- "lmw_est"
#
#   fit$model.matrix <- obj$X
#
#   if (isTRUE(robust)) {
#     if (is.null(cluster)) robust <- "HC3"
#     else robust <- "HC1"
#   }
#   else if (isFALSE(robust)) {
#     if (is.null(cluster)) robust <- "const"
#     else {
#       robust <- "HC1"
#       warning("Setting robust = \"HC1\" because 'cluster' is non-NULL.", call. = FALSE)
#     }
#   }
#   else if (!is.character(robust) || length(robust) != 1 ||
#            !robust %in% eval(formals(sandwich::vcovHC.default)$type)) {
#     stop("'robust' must be TRUE, FALSE, or one of the allowable inputs to the 'type' argument of sandwich::vcovHC().", call. = FALSE)
#   }
#
#   if (is.null(cluster)) {
#     fit$vcov <- sandwich::vcovHC(fit, type = robust, ...)
#   }
#   else {
#     if (inherits(cluster, "formula")) {
#       cluster <- model.frame(cluster,
#                              data = data[pos_w,, drop = FALSE],
#                              na.action = na.pass)
#     }
#     else {
#       cluster <- as.data.frame(cluster)
#     }
#
#     if (nrow(cluster) == nrow(data)) cluster <- cluster[pos_w,, drop = FALSE]
#     else if (nrow(cluster) != length(pos_w)) {
#       stop("'cluster' must have the same number of rows as the original data set.", call. = FALSE)
#     }
#     fit$vcov <- sandwich::vcovCL(fit, type = robust, cluster = cluster, ...)
#   }
#
#   fit$call <- call
#   fit$estimand <- estimand
#   fit$focal <- focal
#   fit$method <- method
#   fit$robust <- robust
#   fit$outcome <- outcome_name
#   fit$treat_levels <- levels(treat)
#
#   fit
# }
