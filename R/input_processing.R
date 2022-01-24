process_base.weights <- function(base.weights, obj) {
  if (is.null(base.weights)) {
    if (!is.null(obj) &&
        (inherits(obj, "matchit") || inherits(obj, "weightit")) &&
        !is.null(obj$weights) && is.numeric(obj$weights) &&
        abs(diff(range(obj$weights))) > sqrt(.Machine$double.eps)) {
      base.weights <- obj$weights
      attr(base.weights, "origin") <- if (inherits(obj, "matchit")) "matchit" else "weightit"
    }
  }
  else if (!is.numeric(base.weights)) {
    stop("'base.weights' must be a numeric vector.", call. = FALSE)
  }

  return(base.weights)
}

process_s.weights <- function(s.weights, obj) {
  if (is.null(s.weights)) {
    if (!is.null(obj) &&
        (inherits(obj, "matchit") || inherits(obj, "weightit")) &&
        !is.null(obj$s.weights) && is.numeric(obj$s.weights) &&
        abs(diff(range(obj$s.weights))) > sqrt(.Machine$double.eps)) {
      s.weights <- obj$s.weights
      attr(s.weights, "origin") <- if (inherits(obj, "matchit")) "matchit" else "weightit"
    }
  }
  else if (!is.numeric(s.weights)) {
    stop("'s.weights' must be a numeric vector.", call. = FALSE)
  }

  return(s.weights)
}

process_dr.method <- function(dr.method, base.weights, method) {
  if (is.null(base.weights)) return(NULL)
  if (length(dr.method) != 1 || !is.character(dr.method)) {
    stop("'dr.method' must be a string.", call. = FALSE)
  }
  dr.method <- toupper(dr.method)
  dr.method <- match_arg(dr.method, c("WLS"))
  # dr.method <- match_arg(dr.method, c("WLS", "AIPW"[method == "MRI"]))

  dr.method
}

process_treat <- function(treat_name, data, multi.ok = TRUE) {

  treat <- model.response(model.frame(reformulate("0", treat_name),
                                      data = data, na.action = "na.pass"))

  if (anyNA(treat)) stop("Missing values are not allowed in the treatment.", call. = FALSE)

  unique_treat <- unique(treat)

  if (length(unique_treat) == 2) {
    if (is.factor(treat)) treat <- factor(treat, levels = levels(treat)[levels(treat) %in% unique_treat])
    else if (is.numeric(treat) && !all(treat == 0 | treat == 1)) {
      stop("If the treatment is not a 0/1 variable, it must be a factor variable.", call. = FALSE)
    }
    else treat <- factor(treat, levels = sort(unique_treat))
  }
  else if (multi.ok) {
    if (is.character(treat)) treat <- factor(treat, levels = sort(unique_treat))
    else if (is.factor(treat)) treat <- factor(treat, levels = levels(treat)[levels(treat) %in% unique_treat])
    else {
      stop("The treatment must be a factor variable if it takes on more than two values.", call. = FALSE)
    }
  }
  else {
    stop("The treatment must be binary.", call. = FALSE)
  }

  attr(treat, "treat_name") <- treat_name

  return(treat)
}

process_estimand <- function(estimand, target, obj) {
  #Get lmw() call to see if estimand arg was specified; if not, it may
  #be replaced by estimand in obj
  m <- match.call(sys.function(sys.parent()),
                  sys.call(sys.parent()))

  estimand.supplied <- hasName(m, "estimand")

  if (!is.character(estimand) || length(estimand) != 1) {
    stop("'estimand' must be a string of length 1.", call. = FALSE)
  }
  estimand <- toupper(estimand)
  estimand <- match_arg(estimand, c("ATE", "ATT", "ATC", "CATE"))

  if (estimand == "CATE") {
    if (is.null(target)) {
      stop("'target' must be specified when estimand = \"CATE\".", call. = FALSE)
    }
  }
  else if (!is.null(target)) {
    if (estimand.supplied) {
      warning("Setting 'estimand' to \"CATE\" because 'target' was supplied.", call. = FALSE)
    }
    estimand <- "CATE"
  }

  if (!is.null(obj) && inherits(obj, "matchit") || inherits(obj, "weightit")) {
    if (!estimand.supplied) {
      estimand <- obj$estimand
    }
    else if (!identical(estimand, obj$estimand)) {
      warning(sprintf("'estimand' (\"%s\") does not agree with the estimand specified in the supplied %s object (\"%s\"). Using \"%s\".",
                      estimand, if (inherits(obj, "matchit")) "matchit" else "weightit", obj$estimand, estimand),
              call. = FALSE)
    }
  }

  estimand
}

process_mf <- function(mf) {
  for (i in seq_len(ncol(mf))) {
    if (anyNA(mf[[i]])) stop("Missing values are not allowed in the covariates.", call. = FALSE)
    if (is.character(mf[[i]])) mf[[i]] <- factor(mf[[i]])
    else if (any(!is.finite(mf[[i]]))) stop("Non-finite values are not allowed in the covariates.", call. = FALSE)
  }
  mf
}

process_data <- function(data = NULL, obj = NULL) {

  obj.data <- NULL
  if (inherits(obj, "matchit")) {
    if (!requireNamespace("MatchIt", quietly = TRUE)) {
      warning("The 'MatchIt' package should be installed when a matchit object is supplied to 'obj'.", call. = FALSE)
    }
    else {
      obj.data <- MatchIt::match.data(obj, drop.unmatched = FALSE, include.s.weights = FALSE)
    }
  }

  if (is.null(data)) {
    if (!is.null(obj.data)) {
      data <- obj.data
    }
  }
  else {
    if (is.matrix(data)) {
      data <- as.data.frame(data)
    }
    if (!is.data.frame(data)) {
      stop("'data' must be a data.frame object.", call. = FALSE)
    }
    if (!is.null(obj.data)) {
      data <- cbind(data, obj.data[setdiff(names(obj.data), names(data))])
    }
  }
  return(data)
}

process_treat_name <- function(treat, formula, method, obj) {
  tt.factors <- attr(terms(formula), "factors")

  obj_treat <- NULL
  if (inherits(obj, "matchit")) {
    obj_treat <- deparse1(as.list(obj$formula)[[2]])
  }

  if (is.null(treat)) {
    if (is.null(obj_treat)) {
      if (!any(colSums(tt.factors != 0) == 1)) {
        stop("Please supply an argument to 'treat' to identify the treatment variable.", call. = FALSE)
      }

      #Use first variable in formula that doesn't involve an interaction
      treat <- colnames(tt.factors)[which(colSums(tt.factors != 0) == 1)[1]]
      message(paste0("Using \"", treat, "\" as the treatment variable. If this is incorrect or to suppress this message, please supply an argument to 'treat' to identify the treatment variable."))
    }
    else {
      if (method == "URI" && !obj_treat %in% rownames(tt.factors)) {
        stop(sprintf("The treatment variable in the supplied matchit object (%s) does not align with any variables in 'formula'.",
                     obj_treat), call. = FALSE)
      }
      treat <- obj_treat
    }
  }
  else {
    if (length(treat) != 1 || !is.character(treat)) {
      stop("'treat' must be a string naming the treatment variable.", call. = FALSE)
    }
    if (method == "URI" && !treat %in% rownames(tt.factors)) {
      stop(sprintf("The supplied treatment variable (\"%s\") does not align with any variables in 'formula'.",
                   treat), call. = FALSE)
    }
  }
  return(treat)
}

process_contrast <- function(contrast = NULL, treat, method) {

  treat_f <- if (is.factor(treat)) droplevels(treat) else as.factor(treat)
  t_levels <- levels(treat_f)

  if (is.null(contrast)) {
    if (is.null(contrast) && length(t_levels) > 2 && method == "URI") {
      stop("'contrast' must be specified when the treatment has more than two levels and method = \"URI\".", call. = FALSE)
    }
    return(NULL)
  }

  if (is.numeric(contrast)) {
    if (can_str2num(treat)) {
      contrast <- as.character(contrast)
    }
    else if (!all(contrast %in% seq_along(t_levels))) {
      stop("'contrast' must contain the names or indices of treatment levels to be contrasted.", call. = FALSE)
    }
    else {
      contrast <- t_levels[contrast]
    }
  }
  if (is.factor(contrast)) {
    contrast <- as.character(contrast)
  }

  #contrast is now character
  if (!all(contrast %in% t_levels)) {
    stop("'contrast' must contain the names or indices of treatment levels to be contrasted.", call. = FALSE)
  }

  if (length(contrast) == 1) {
    if (length(t_levels) == 2) {
      contrast <- c(contrast, t_levels[t_levels != contrast])
    }
    else if (contrast == t_levels[1]) {
      stop("If 'contrast' is a single value, it cannot be the reference value of the treatment.", call. = FALSE)
    }
    else {
      contrast <- c(contrast, t_levels[1])
    }
  }
  else if (length(contrast) != 2) {
    stop("'contrast' cannot have length greater than 2.", call. = FALSE)
  }

  return(contrast)
}

process_focal <- function(focal = NULL, treat, estimand) {
  if (estimand %in% c("ATE", "CATE")) {
    if (!is.null(focal)) warning(sprintf("'focal' is ignored when estimand = \"%s\".", estimand), call. = FALSE)
    return(NULL)
  }

  if (is.null(focal)) {
    focal <- switch(estimand, "ATT" = levels(treat)[2], levels(treat)[1])
    if (nlevels(treat) > 2 || !can_str2num(unique(treat, nmax = 2)))
      message(sprintf("Using \"%s\" as the focal (%s) group. If this is incorrect or to suppress this message, please supply an argument to 'focal' to identify the focal treatment level.",
                      focal, switch(estimand, "ATT" = "treated", "control")))
  }
  else if (length(focal) != 1) {
    stop("'focal' must be of length 1.", call. = FALSE)
  }
  else if (!as.character(focal) %in% levels(treat)) {
    stop("'focal' must be the name of a value of the treatment variable.", call. = FALSE)
  }

  return(as.character(focal))
}

apply_contrast_to_treat <- function(treat, contrast = NULL) {
  if (is.null(contrast)) return(treat)

  #Re-order levels of treat; use treat[] to retain attributes (e.g., "treat_name")
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
      stop("A valid dataset could not be found. Please supply an argument to 'data' containing the original dataset used to estimate the weights.", call. = FALSE)
    }
  }
  else {
    if (!is.data.frame(data)) {
      if (is.matrix(data)) data <- as.data.frame.matrix(data)
      else stop("'data' must be a data frame.", call. = FALSE)
    }
    if (nrow(data) != length(x$treat)) {
      stop("'data' must have as many rows as there were units in the original call to lmw().", call. = FALSE)
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
#    do.call("get_outcome", list(substitute(outcome_arg), data_arg, formula_args))
get_outcome <- function(outcome, data = NULL, formula) {

  tt <- terms(formula, data = data)

  if (missing(outcome)) {
    if (attr(tt, "response") == 0) {
      stop("'outcome' must be supplied.", call. = FALSE)
    }

    outcome_char <- deparse1(attr(terms(tt), "variables")[[2]])

    mf <- try(eval(model.frame(tt, data = data)), silent = TRUE)
    if (inherits(mf, "try-error")) {
      cond <- attr(mf, "condition")$message
      if (startsWith(cond, "object") && endsWith(cond, "not found")) {
        stop(sprintf("The outcome variable '%s' must be present in the supplied dataset or environment.", outcome_char), call. = FALSE)
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
          stop(sprintf("The outcome variable '%s' cannot befound in the environment. Please supply an argument to 'data' containing the original dataset used to estimate the weights.", outcome_char),
               call. = FALSE)
        }
        else {
          stop(sprintf("The outcome variable '%s' must be present in the supplied dataset or environment.", outcome_char), call. = FALSE)
        }
      }
    }
    if (is.character(outcome) && length(outcome) == 1) {
      if (is.null(data)) {
        stop("A dataset must be present when 'outcome' is supplied as a string. Please supply an argument to 'data' containing the original dataset used to estimate the weights.", call. = FALSE)
      }
      outcome_char <- outcome
      outcome <- try(eval(str2expression(outcome_char), data), silent = TRUE)
      if (length(outcome) == 0 || inherits(outcome, "try-error")) {
        stop("The outcome variable must be present in the dataset.", call. = FALSE)
      }
    }
  }

  if (length(outcome) == 0) {
    stop("The outcome variable cannot be NULL.", call. = FALSE)
  }
  if (!is.numeric(outcome) && !is.logical(outcome)) {
    print(outcome)
    stop("The outcome variable must be numeric.", call. = FALSE)
  }
  if (length(outcome) != nrow(data)) {
    stop("The outcome variable must have length equal to the number of units in the dataset.", call. = FALSE)
  }

  attr(outcome, "outcome_name") <- outcome_char
  outcome
}

process_target <- function(target, formula, mf) {
  if (any(c("$", "[", "[[") %in% all.names(formula))) {
    stop("Subsetting operations ($, [.], [[.]]) are not allowed in the model formula when 'target' is specified.", call. = FALSE)
  }
  if (!is.list(target) && !is.data.frame(target)) {
    stop("'target' must be a list of covariate-value pairs.", call. = FALSE)
  }
  if (!all(lengths(target) == 1L)) {
    stop("All entries in 'target' must have lengths of 1.", call. = FALSE)
  }
  vars_in_formula <- all.vars(formula)
  vars_in_target <- names(target)
  vars_in_formula_not_in_target <- setdiff(vars_in_formula, vars_in_target)
  if (length(vars_in_formula_not_in_target) > 0) {
    stop(paste0("All covariates in the model formula must be present in 'target'; variable(s) not present:\n\t",
                paste(vars_in_formula_not_in_target, collapse = ", ")), call. = FALSE)
  }
  vars_in_target_not_in_formula <- setdiff(vars_in_target, vars_in_formula)
  if (length(vars_in_target_not_in_formula) > 0) {
    warning(paste0("The following value(s) in 'target' will be ignored:\n\t",
                   paste(vars_in_target_not_in_formula, collapse = ", ")), call. = FALSE)
    target <- target[vars_in_target %in% vars_in_formula]
  }

  target <- process_mf(as.data.frame(target))
  .checkMFClasses(vapply(mf, .MFclass, character(1L)), target)

  for (i in names(target)[vapply(target, is.factor, logical(1L))]) {
    target[[i]] <- factor(as.character(target[[i]]), levels = levels(mf[[i]]))
  }

  target_mf <- model.frame(formula, data = target)
  target_mm <- model.matrix(formula, data = target_mf)[,-1, drop = FALSE]

  out <- setNames(drop(target_mm), colnames(target_mm))
  attr(out, "target_original") <- target

  out
}

process_iv_name <- function(iv, formula) {
  tt.factors <- attr(terms(formula), "factors")

  if (missing(iv) || is.null(iv)) {
    stop("Please supply an argument to 'iv' to identify the instrumental variable.", call. = FALSE)
  }

  if (length(iv) != 1 || !is.character(iv)) {
    stop("'iv' must be a string naming the instrumental variable.", call. = FALSE)
  }
  if (iv %in% rownames(tt.factors)) {
    stop("The instrumental variable should not be present in the model formula.", call. = FALSE)
  }

  return(iv)
}

process_iv <- function(iv_name, data) {
  iv <- model.response(model.frame(reformulate("0", iv_name),
                                   data = data, na.action = "na.pass"))

  if (anyNA(iv)) stop("Missing values are not allowed in the instrumental variable.", call. = FALSE)

  if (is.logical(iv)) iv <- as.numeric(iv)
  if (!is.numeric(iv)) {
    stop("The instrumental variable must be numeric, not a factor or character.", call. = FALSE)
  }

  attr(iv, "iv_name") <- iv_name

  return(iv)
}
