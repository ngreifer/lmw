process_base.weights <- function(base.weights = NULL, data = NULL, obj = NULL) {

  base.weights_sub <- substitute(base.weights)
  base.weights_char <- deparse1(base.weights_sub)
  base.weights <- try(eval(base.weights_sub, data, parent.frame(2)), silent = TRUE)
  if (inherits(base.weights, "try-error")) {
    cond <- conditionMessage(attr(base.weights, "condition"))
    if (startsWith(cond, "object") && endsWith(cond, "not found")) {
      if (is.null(data)) {
        chk::err(sprintf("the base weights variable '%s' cannot be found in the environment. Please supply an argument to `data` containing the base weights",
                     base.weights_char))
      }
      else {
        chk::err(sprintf("the base weights variable '%s' must be present in the supplied dataset or environment",
                     base.weights_char))
      }
    }
    else {
      chk::err(cond, tidy = FALSE)
    }
  }
  else if (length(base.weights) == 0) {
    if (!is.null(obj) &&
        (inherits(obj, "matchit") || inherits(obj, "weightit")) &&
        !is.null(obj$weights) && is.numeric(obj$weights) &&
        abs(diff(range(obj$weights))) > sqrt(.Machine$double.eps)) {
      base.weights <- obj$weights
      attr(base.weights, "origin") <- if (inherits(obj, "matchit")) "matchit" else "weightit"
    }
  }
  else if (is.character(base.weights) && length(base.weights) == 1) {
    if (is.null(data)) {
      chk::err("a dataset must be present when `base.weights` is supplied as a string. Please supply an argument to `data` containing the base weights")
    }
    base.weights_char <- base.weights
    base.weights <- try(eval(str2expression(base.weights_char), data), silent = TRUE)
    if (length(base.weights) == 0 || inherits(base.weights, "try-error")) {
      chk::err("the base weights variable must be present in the dataset")
    }
  }

  if (length(base.weights) == 0) {
    return(NULL)
  }
  if (!is.numeric(base.weights)) {
    chk::err("the base weights variable must be numeric")
  }

  base.weights
}

process_s.weights <- function(s.weights = NULL, data = NULL, obj = NULL) {

  s.weights_sub <- substitute(s.weights)
  s.weights_char <- deparse1(s.weights_sub)
  s.weights <- try(eval(s.weights_sub, data, parent.frame(2)), silent = TRUE)
  if (inherits(s.weights, "try-error")) {
    cond <- conditionMessage(attr(s.weights, "condition"))
    if (startsWith(cond, "object") && endsWith(cond, "not found")) {
      if (is.null(data)) {
        chk::err(sprintf("the sampling weights variable '%s' cannot be found in the environment. Please supply an argument to `data` containing the sampling weights", s.weights_char))
      }
      else {
        chk::err(sprintf("the sampling weights variable '%s' must be present in the supplied dataset or environment",
                         s.weights_char))
      }
    }
    else {
      chk::err(cond, tidy = FALSE)
    }
  }
  else if (length(s.weights) == 0) {
    if (!is.null(obj) &&
        (inherits(obj, "matchit") || inherits(obj, "weightit")) &&
        !is.null(obj$s.weights) && is.numeric(obj$s.weights) &&
        abs(diff(range(obj$s.weights))) > sqrt(.Machine$double.eps)) {
      s.weights <- obj$s.weights
      attr(s.weights, "origin") <- if (inherits(obj, "matchit")) "matchit" else "weightit"
    }
  }
  else if (is.character(s.weights) && length(s.weights) == 1) {
    if (is.null(data)) {
      chk::err("a dataset must be present when `s.weights` is supplied as a string. Please supply an argument to `data` containing the sampling weights")
    }
    s.weights_char <- s.weights
    s.weights <- try(eval(str2expression(s.weights_char), data), silent = TRUE)
    if (length(s.weights) == 0 || inherits(s.weights, "try-error")) {
      chk::err("the sampling weights variable must be present in the dataset")
    }
  }

  if (length(s.weights) == 0) {
    return(NULL)
  }
  if (!is.numeric(s.weights)) {
    chk::err("the sampling weights variable must be numeric")
  }

  s.weights
}

process_dr.method <- function(dr.method, base.weights, method, estimand) {
  if (is.null(base.weights)) return(NULL)

  chk::chk_string(dr.method)
  dr.method <- toupper(dr.method)
  if (dr.method == "IPWRA") dr.method <- "WLS"
  # dr.method <- match_arg(dr.method, c("WLS"))
  # dr.method <- match_arg(dr.method, c("WLS", "AIPW"[method == "MRI"]))
  dr.method <- match_arg(dr.method, c("WLS", "AIPW"))

  if (estimand == "CATE" && dr.method == "AIPW") {
    chk::err("the CATE cannot be used with AIPW")
  }

  dr.method
}

process_treat <- function(treat_name, data, multi.ok = TRUE) {

  treat <- model.response(model.frame(reformulate("0", treat_name),
                                      data = data, na.action = "na.pass"))

  if (anyNA(treat)) {
    chk::err("missing values are not allowed in the treatment")
  }

  unique_treat <- unique(treat)

  if (length(unique_treat) == 2) {
    if (is.factor(treat)) treat <- factor(treat, levels = levels(treat)[levels(treat) %in% unique_treat])
    else if (!is.numeric(treat) || all(treat == 0 | treat == 1)) {
      treat <- factor(treat, levels = sort(unique_treat))
    }
    else {
      chk::err("if the treatment is not a 0/1 variable, it must be a factor variable")
    }
  }
  else if (multi.ok) {
    if (is.character(treat)) treat <- factor(treat, levels = sort(unique_treat))
    else if (is.factor(treat)) treat <- factor(treat, levels = levels(treat)[levels(treat) %in% unique_treat])
    else {
      chk::err("the treatment must be a factor variable if it takes on more than two values")
    }
  }
  else {
    chk::err("the treatment must be binary")
  }

  attr(treat, "treat_name") <- treat_name

  treat
}

process_estimand <- function(estimand, target, obj) {
  #Get lmw() call to see if estimand arg was specified; if not, it may
  #be replaced by estimand in obj
  m <- match.call(sys.function(sys.parent()),
                  sys.call(sys.parent()))

  estimand.supplied <- utils::hasName(m, "estimand")

  chk::chk_string(estimand)
  estimand <- toupper(estimand)
  estimand <- match_arg(estimand, c("ATE", "ATT", "ATC", "CATE"))

  if (estimand == "CATE") {
    if (is.null(target)) {
      chk::err("`target` must be specified when `estimand = \"CATE\"`")
    }
  }
  else if (!is.null(target)) {
    if (estimand.supplied) {
      chk::wrn("setting `estimand` to \"CATE\" because `target` was supplied")
    }
    estimand <- "CATE"
  }

  if (!is.null(obj) && (inherits(obj, "matchit") || inherits(obj, "weightit"))) {
    if (!estimand.supplied) {
      estimand <- obj$estimand
    }
    else if (!identical(estimand, obj$estimand)) {
      chk::wrn(sprintf("`estimand` (\"%s\") does not agree with the estimand specified in the supplied %s object (\"%s\") Using \"%s\".",
                      estimand, if (inherits(obj, "matchit")) "matchit" else "weightit", obj$estimand, estimand))
    }
  }

  estimand
}

process_mf <- function(mf) {
  for (i in seq_len(ncol(mf))) {
    if (is.character(mf[[i]])) mf[[i]] <- factor(mf[[i]])
    else if (any(!is.finite(mf[[i]]))) {
      chk::err("non-finite values are not allowed in the covariates")
    }
  }
  mf
}

process_data <- function(data = NULL, obj = NULL) {

  null.data <- is.null(data)
  if (!null.data) {
    if (is.matrix(data)) {
      data <- as.data.frame(data)
    }
    else if (!is.data.frame(data)) {
      chk::err("`data` must be a data.frame object")
    }
  }

  obj.data <- NULL
  if (inherits(obj, "matchit")) {
    if (!requireNamespace("MatchIt", quietly = TRUE)) {
      if (null.data) {
        chk::wrn("The 'MatchIt' package should be installed when a matchit object is supplied to `obj`")
      }
    }
    else {
      obj.data <- MatchIt::match.data(obj, drop.unmatched = FALSE, include.s.weights = FALSE)
    }
  }

  if (!null.data) {
    if (!is.null(obj.data)) {
      data <- cbind(data, obj.data[setdiff(names(obj.data), names(data))])
    }
  }
  else {
    data <- obj.data #NULL if no obj
  }

  data
}

process_treat_name <- function(treat, formula, data, method, obj) {
  tt.factors <- attr(terms(formula, data = data), "factors")

  obj_treat <- NULL
  if (inherits(obj, "matchit")) {
    obj_treat <- deparse1(as.list(obj$formula)[[2]])
  }

  if (is.null(treat)) {
    if (is.null(obj_treat)) {
      if (!any(colSums(tt.factors != 0) == 1)) {
        chk::err("please supply an argument to `treat` to identify the treatment variable")
      }

      #Use first variable in formula that doesn't involve an interaction
      treat <- colnames(tt.factors)[which(colSums(tt.factors != 0) == 1)[1]]
      chk::msg(sprintf("using \"%s\" as the treatment variable. If this is incorrect or to suppress this message, please supply an argument to 'treat' to identify the treatment variable", treat))
    }
    else {
      if (method == "URI" && !obj_treat %in% rownames(tt.factors)) {
        chk::err(sprintf("the treatment variable in the supplied matchit object (%s) does not align with any variables in `formula`",
                     obj_treat))
      }
      treat <- obj_treat
    }
  }
  else {
    if (length(treat) != 1 || !is.character(treat)) {
      chk::err("`treat` must be a string naming the treatment variable")
    }
    if (method == "URI" && !treat %in% rownames(tt.factors)) {
      chk::err(sprintf("the supplied treatment variable (\"%s\") does not align with any variables in `formula`",
                   treat))
    }
  }

  treat
}

process_contrast <- function(contrast = NULL, treat, method) {

  treat_f <- if (is.factor(treat)) droplevels(treat) else as.factor(treat)
  t_levels <- levels(treat_f)

  if (is.null(contrast)) {
    if (is.null(contrast) && length(t_levels) > 2 && method == "URI") {
      chk::err("`contrast` must be specified when the treatment has more than two levels and `method = \"URI\"`")
    }
    return(NULL)
  }

  if (is.numeric(contrast)) {
    if (can_str2num(treat)) {
      contrast <- as.character(contrast)
    }
    else if (all(contrast %in% seq_along(t_levels))) {
      contrast <- t_levels[contrast]
    }
    else {
      chk::err("`contrast` must contain the names or indices of treatment levels to be contrasted")
    }
  }

  if (is.factor(contrast)) {
    contrast <- as.character(contrast)
  }

  #contrast is now character
  if (!all(contrast %in% t_levels)) {
    chk::err("`contrast` must contain the names or indices of treatment levels to be contrasted")
  }

  if (length(contrast) == 1) {
    if (length(t_levels) == 2) {
      contrast <- c(contrast, t_levels[t_levels != contrast])
    }
    else if (contrast != t_levels[1]) {
      contrast <- c(contrast, t_levels[1])
    }
    else {
      chk::err("if `contrast` is a single value, it cannot be the reference value of the treatment")
    }
  }
  else if (length(contrast) != 2) {
    chk::err("`contrast` cannot have length greater than 2")
  }

  contrast
}

process_focal <- function(focal = NULL, treat, estimand, obj = NULL) {
  if (estimand %in% c("ATE", "CATE")) {
    if (!is.null(focal)) chk::wrn(sprintf("`focal` is ignored when `estimand = \"%s\"`", estimand))
    return(NULL)
  }

  if (!is.null(focal)) {
    if (length(focal) != 1) {
      chk::err("`focal` must be of length 1")
    }

    if (!as.character(focal) %in% levels(treat)) {
      chk::err("`focal` must be the name of a value of the treatment variable")
    }

    if (!is.null(obj) && !is.null(obj$focal) && !identical(focal, obj$focal)) {
      chk::wrn("`focal` does not align with the `focal` component of the `obj` input")
    }
  }
  else {
    if (is.null(obj$focal) || is.null(obj$focal)) {
      focal <- switch(estimand, "ATT" = levels(treat)[2], levels(treat)[1])
      if (nlevels(treat) > 2 || !can_str2num(unique(treat, nmax = 2)))
        chk::msg(sprintf("using \"%s\" as the focal (%s) group. If this is incorrect or to suppress this message, please supply an argument to `focal` to identify the focal treatment level",
                         focal, switch(estimand, "ATT" = "treated", "control")))
    }
    else {
      focal <- obj$focal

      if (length(focal) != 1) {
        chk::err("`focal` must be of length 1")
      }

      if (!as.character(focal) %in% levels(treat)) {
        chk::err("`focal` must be the name of a value of the treatment variable")
      }
    }
  }

  as.character(focal)
}

apply_contrast_to_treat <- function(treat, contrast = NULL) {
  if (is.null(contrast)) return(treat)

  attrs <- attributes(treat)
  treat <- factor(treat, levels = c(rev(contrast), levels(treat)[!levels(treat) %in% contrast]))
  for (i in setdiff(names(attrs), "levels")) attr(treat, i) <- attrs[[i]]

  treat
}

#Extract data from an lmw object
get_data <- function(data, x) {
  #x is a lmw object

  if (is.null(data)) {
    f_env <- environment(x$formula)
    data <- try(eval(x$call$data, envir = f_env), silent = TRUE)

    if (length(data) == 0 || inherits(data, "try-error") || length(dim(data)) != 2 || nrow(data) != length(x[["treat"]])) {
      data <- try(eval(x$call$data, envir = parent.frame(2)), silent = TRUE)
    }

    obj <- try(eval(x$call$obj, envir = f_env), silent = TRUE)
    if (length(obj) == 0 || inherits(obj, "try-error") || (!inherits(obj, "matchit") && !inherits(obj, "weightit"))) {
      obj <- try(eval(x$call$obj, envir = parent.frame(2)), silent = TRUE)
    }

    data <- process_data(data, obj)

    if (length(data) != 0 && (length(dim(data)) != 2 || nrow(data) != length(x[["treat"]]))) {
      chk::err("a valid dataset could not be found. Please supply an argument to `data` containing the original dataset used to estimate the weights")
    }
  }
  else {
    if (!is.data.frame(data)) {
      if (!is.matrix(data)) chk::err("`data` must be a data frame")
      data <- as.data.frame.matrix(data)
    }

    if (nrow(data) != length(x$treat)) {
      chk::err("`data` must have as many rows as there were units in the original call to `lmw()`")
    }

    obj <- try(eval(x$call$obj, envir = f_env), silent = TRUE)
    if (length(obj) == 0 || inherits(obj, "try-error") || (!inherits(obj, "matchit") && !inherits(obj, "weightit"))) {
      obj <- try(eval(x$call$obj, envir = parent.frame(2)), silent = TRUE)
    }

    data <- process_data(data, obj)
  }

  data
}

#Get outcome from 'outcome' or formula
# 'outcome' can be a string or the name of a variable in 'data' or the environment
# containing the outcome
# Uses NSE, so must be called with this idiom when used inside functions:
#    do.call("get_outcome", list(substitute(outcome_arg), data_arg, formula_arg))
get_outcome <- function(outcome, data = NULL, formula, X) {

  tt <- terms(formula, data = data)

  if (missing(outcome)) {
    if (attr(tt, "response") == 0) {
      chk::err("`outcome` must be supplied")
    }

    outcome_char <- deparse1(attr(tt, "variables")[[2]])

    mf <- try(eval(model.frame(tt, data = data)), silent = TRUE)
    if (inherits(mf, "try-error")) {
      cond <- attr(mf, "condition")$message
      if (startsWith(cond, "object") && endsWith(cond, "not found")) {
        chk::err(sprintf("the outcome variable '%s' must be present in the supplied dataset or environment", outcome_char))
      }
    }
    outcome <- model.response(mf)
  }
  else {
    outcome_sub <- substitute(outcome)
    outcome_char <- deparse1(outcome_sub)
    outcome <- try(eval(outcome_sub, envir = data), silent = TRUE)
    if (inherits(outcome, "try-error")) {
      cond <- attr(outcome, "condition")$message
      if (startsWith(cond, "object") && endsWith(cond, "not found")) {
        if (is.null(data)) {
          chk::err(sprintf("the outcome variable '%s' cannot be found in the environment. Please supply an argument to `data` containing the original dataset used to estimate the weights", outcome_char))
        }
        else {
          chk::err(sprintf("the outcome variable '%s' must be present in the supplied dataset or environment", outcome_char))
        }
      }
    }
    if (is.character(outcome) && length(outcome) == 1) {
      if (is.null(data)) {
        chk::err("a dataset must be present when `outcome` is supplied as a string. Please supply an argument to `data` containing the original dataset used to estimate the weights")
      }
      outcome_char <- outcome
      outcome <- try(eval(str2expression(outcome_char), data), silent = TRUE)
      if (length(outcome) == 0 || inherits(outcome, "try-error")) {
        chk::err("the outcome variable must be present in the dataset")
      }
    }
  }

  if (length(outcome) == 0) {
    chk::err("the outcome variable cannot be `NULL`")
  }
  if (!is.numeric(outcome) && !is.logical(outcome)) {
    chk::err("the outcome variable must be numeric")
  }
  if (length(outcome) != nrow(X)) {
    chk::err("the outcome variable must have length equal to the number of units in the dataset")
  }

  attr(outcome, "outcome_name") <- outcome_char
  outcome
}

process_target <- function(target, formula, mf, target.weights = NULL) {
  if (any(c("$", "[", "[[") %in% all.names(formula))) {
    chk::err("subsetting operations ($, [.], [[.]]) are not allowed in the model formula when `target` is specified")
  }
  if (!is.list(target)) {
    chk::err("`target` must be a list of covariate-value pairs or a data frame containing the target population")
  }
  if (!is.data.frame(target) && !all(lengths(target) == 1L)) {
    chk::err("all entries in `target` must have lengths of 1 when supplied as a list")
  }
  if (!is.null(target.weights)) {
    if (!is.data.frame(target) || nrow(target) == 1) {
      chk::wrn("`target.weights` is ignored when `target` is a target profile")
    }
    if (!is.numeric(target.weights) || length(target.weights) != nrow(target)) {
      chk::err("`target.weights` must be a numeric vector with length equal to the number of rows of the target dataset")
    }
    target.weights <- as.numeric(target.weights)
  }

  vars_in_formula <- all.vars(formula)
  vars_in_target <- names(target)
  vars_in_formula_not_in_target <- setdiff(vars_in_formula, vars_in_target)
  if (length(vars_in_formula_not_in_target) > 0) {
    chk::err(paste0("All covariates in the model formula must be present in `target`; variable(s) not present:\n\t",
                paste(vars_in_formula_not_in_target, collapse = ", ")), tidy = FALSE)
  }

  vars_in_target_not_in_formula <- setdiff(vars_in_target, vars_in_formula)
  if (length(vars_in_target_not_in_formula) > 0) {
    if (!is.data.frame(target)) {
      chk::wrn(paste0("The following value(s) in `target` will be ignored:\n\t",
                     paste(vars_in_target_not_in_formula, collapse = ", ")), tidy = FALSE)
    }
    target <- target[vars_in_target %in% vars_in_formula]
  }

  target <- process_mf(as.data.frame(target))
  .checkMFClasses(vapply(mf, .MFclass, character(1L)), target)

  for (i in names(target)[vapply(target, is.factor, logical(1L))]) {
    target[[i]] <- factor(as.character(target[[i]]), levels = levels(mf[[i]]))
  }

  target_mf <- model.frame(formula, data = target)
  target_mm <- model.matrix(formula, data = target_mf)[,-1, drop = FALSE]

  out <- {
    if (nrow(target_mm) == 1) drop(target_mm)
    else colMeans_w(target_mm, target.weights)
  }

  names(out) <- colnames(target_mm)

  attr(target, "target.weights") <- target.weights
  attr(out, "target_original") <- target

  out
}

process_iv <- function(iv, formula, data = NULL) {
  if (length(iv) == 0) {
    chk::err("an argument to `iv` specifying the instrumental variable(s) is required")
  }

  tt.factors <- attr(terms(formula, data = data), "factors")

  if (is.character(iv)) {
    if (is.null(data) || !is.data.frame(data)) {
      chk::err("if `iv` is specified as a string, a data frame argument must be supplied to `data`")
    }

    if (!all(iv %in% names(data))) {
      chk::err("all variables in `iv` must be in `data`")
    }

    iv_f <- reformulate(iv)
  }
  else if (inherits(iv, "formula")) {
    iv_f <- iv
  }
  else {
    chk::err("`iv` must be a one-sided formula or character vector of instrumental variables")
  }

  iv.factors <- attr(terms(iv_f, data = data), "term.labels")

  if (any(iv.factors %in% rownames(tt.factors))) {
    chk::err("the instrumental variable(s) should not be present in the model formula")
  }

  iv_f
}

process_fixef <- function(fixef, formula, data = NULL, treat_name) {
  if (length(fixef) == 0) return(NULL)

  tt.factors <- attr(terms(formula, data = data), "factors")

  if (is.character(fixef)) {
    if (is.null(data) || !is.data.frame(data)) {
      chk::err("if `fixef` is specified as a string, a data frame argument must be supplied to `data`")
    }

    if (!all(fixef %in% names(data))) {
      chk::err("all variables in `fixef` must be in `data`")
    }

    fixef_f <- reformulate(fixef)
  }
  else if (inherits(fixef, "formula")) {
    fixef_f <- fixef
  }
  else {
    chk::err("`fixef` must be a one-sided formula or string naming the fixed effect variable")
  }

  fixef_name <- attr(terms(fixef_f, data = data), "term.labels")
  if (length(fixef_name) > 1) {
    chk::err("only one fixed effect variable may be supplied")
  }

  if (any(fixef_name == treat_name)) {
    chk::err("the fixed effect variable cannot be the same as the treatment variable")
  }
  if (any(fixef_name %in% rownames(tt.factors))) {
    chk::err("the fixed effect variable should not be present in the model formula")
  }

  fixef_mf <- model.frame(fixef_f, data = data, na.action = "na.pass")

  if (anyNA(fixef_mf)) {
    chk::err("missing values are not allowed in the fixed effect variable")
  }

  fixef <- factor(fixef_mf[[fixef_name]])
  attr(fixef, "fixef_name") <- fixef_name

  fixef
}

check_lengths <- function(treat, ...) {
  arg_names <- c("treatment", unlist(lapply(substitute(list(...))[-1], deparse1)))
  args <- c(list(treat), list(...))

  lengths <- setNames(vapply(args, function(x) {
    if (length(dim(x)) == 2) NROW(x)
    else length(x)
  }, integer(1L)), arg_names)

  lengths <- lengths[lengths > 0L]
  arg_names <- names(lengths)

  if (!all(lengths == lengths[1])) {
    error_df <- data.frame(sort(unique(lengths)))
    rownames(error_df) <- format(sapply(error_df[[1]], function(n) {
      paste0(paste(arg_names[lengths == n], collapse = ", "), ":")
    }), justify = "right")
    names(error_df) <- NULL

    msg <- paste(
      sprintf("%s must have the same length. Variable lengths:",
              word_list(c("the treatment", add_quotes(arg_names[-1], "`")))),
      format(do.call("paste", lapply(rownames(error_df), function(i) {
        paste0("\n", i, " ", error_df[i, 1])
      })), justfify = "right")
    )

    chk::err(msg, tidy = FALSE)
  }

  invisible(treat)
}
