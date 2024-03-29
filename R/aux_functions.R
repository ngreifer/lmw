#More informative and cleaner version of base::match.arg. From WeightIt with edits.
match_arg <- function(arg, choices, several.ok = FALSE) {
  #Replaces match.arg() but gives cleaner error message and processing
  #of arg.
  if (missing(arg))
    stop("No argument was supplied to match_arg.", call. = FALSE)
  arg.name <- deparse1(substitute(arg))

  if (missing(choices)) {
    formal.args <- formals(sys.function(sysP <- sys.parent()))
    choices <- eval(formal.args[[as.character(substitute(arg))]],
                    envir = sys.frame(sysP))
  }

  if (is.null(arg)) return(choices[1L])
  if (!is.character(arg))
    chk::err(sprintf("the argument to `%s` must be NULL or a character vector", arg.name))

  if (!several.ok) {
    if (identical(arg, choices)) return(arg[1L])
    if (length(arg) > 1L) {
      chk::err(sprintf("the argument to `%s` must be of length 1", arg.name))
    }
  }
  else if (length(arg) == 0) {
    chk::err(sprintf("the argument to `%s` must be of length >= 1", arg.name))
  }

  i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
  if (all(i == 0L))
    chk::err(sprintf("the argument to `%s` should be %s%s",
                     arg.name,
                     ngettext(length(choices), "", if (several.ok) "at least one of " else "one of "),
                     word_list(choices, and.or = "or", quotes = 2)))

  i <- i[i > 0L]

  choices[i]
}

#Function to turn a vector into a string with "," and "and" or "or" for clean messages. 'and.or'
#controls whether words are separated by "and" or "or"; 'is.are' controls whether the list is
#followed by "is" or "are" (to avoid manually figuring out if plural); quotes controls whether
#quotes should be placed around words in string. From WeightIt.
word_list <- function(word.list = NULL, and.or = c("and", "or"), is.are = FALSE, quotes = FALSE) {
  #When given a vector of strings, creates a string of the form "a and b"
  #or "a, b, and c"
  #If is.are, adds "is" or "are" appropriately
  L <- length(word.list)
  word.list <- add_quotes(word.list, quotes)

  if (L == 0) {
    out <- ""
    attr(out, "plural") <- FALSE
  }
  else {
    word.list <- word.list[!word.list %in% c(NA_character_, "")]
    L <- length(word.list)
    if (L == 0) {
      out <- ""
      attr(out, "plural") <- FALSE
    }
    else if (L == 1) {
      out <- word.list
      if (is.are) out <- paste(out, "is")
      attr(out, "plural") <- FALSE
    }
    else {
      and.or <- match_arg(and.or, c("and", "or"))
      if (L == 2) {
        out <- paste(word.list, collapse = paste0(" ", and.or," "))
      }
      else {
        out <- paste(paste(word.list[seq_len(L-1)], collapse = ", "),
                     word.list[L], sep = paste0(", ", and.or," "))

      }
      if (is.are) out <- paste(out, "are")
      attr(out, "plural") <- TRUE
    }

  }

  out
}

#Add quotation marks around a string.
add_quotes <- function(x, quotes = 2L) {
  if (!isFALSE(quotes)) {
    if (is.character(quotes)) x <- paste0(quotes, x, quotes)
    else if (isTRUE(quotes) || as.integer(quotes) == 2L) x <- paste0("\"", x, "\"")
    else if (as.integer(quotes) == 1L) x <- paste0("\'", x, "\'")
    else stop("'quotes' must be a string, boolean, 1, or 2.")
  }
  x
}

#Effective sample size
ESS <- function(w) {
  # sum(abs(w))^2/sum(w^2)
  sum(w)^2/sum(w^2)
}

#Weighted colMeans
colMeans_w <- function(mat, w = NULL, subset = NULL) {
  if (length(subset) != 0) {
    mat <- mat[subset,,drop = FALSE]
    if (length(w) != 0) w <- w[subset]
  }

  if (length(w) == 0) return(colSums(mat)/nrow(mat))

  colSums(mat * w)/sum(w)
}

#Weighted mean (faster than weighted.mean())
mean_w <- function(x, w = NULL, subset = NULL) {
  if (length(subset) != 0) {
    x <- x[subset]
    if (length(w) != 0) w <- w[subset]
  }

  if (length(w) == 0) return(sum(x)/length(x))

  sum(x * w)/sum(w)
}

#(Weighted) variance that uses special formula for binary variables
var_w <- function(x, bin.var = NULL, w = NULL, subset = NULL) {
  if (is.null(bin.var)) bin.var <- all(x == 0 | x == 1)
  if (length(subset) != 0) {
    x <- x[subset]
    if (length(w) != 0) w <- w[subset]
  }

  if (is.null(w)) w <- rep(1, length(x))

  w <- w / sum(w) #weights normalized to sum to 1
  mx <- sum(w * x) #weighted mean

  if (bin.var) {
    v <- mx*(1-mx)
  }
  else {
    #Reliability weights variance; same as cov.wt()
    v <- sum(w * (x - mx)^2)/(1 - sum(w^2))
  }
  abs(v)
}

#Determine whether a character vector can be coerced to numeric
can_str2num <- function(x) {
  nas <- is.na(x)
  suppressWarnings(x_num <- as.numeric(as.character(x[!nas])))

  !anyNA(x_num)
}

#Cleanly coerces a character vector to numeric; best to use after can_str2num()
str2num <- function(x) {
  nas <- is.na(x)
  suppressWarnings(x_num <- as.numeric(as.character(x)))
  x_num[nas] <- NA

  x_num
}

#Clean printing of data frames with numeric and NA elements.
round_df_char <- function(df, digits, pad = "0", na_vals = ".") {
  #Digits is passed to round(). pad is used to replace trailing zeros so decimal
  #lines up. Should be "0" or " "; "" (the empty string) un-aligns decimals.
  #na_vals is what NA should print as.

  if (NROW(df) == 0 || NCOL(df) == 0) return(df)
  if (!is.data.frame(df)) {
    df <- as.data.frame.matrix(df, stringsAsFactors = FALSE)
  }
  rn <- rownames(df)
  cn <- colnames(df)

  infs <- array(FALSE, dim = dim(df))
  # o.negs <- array(FALSE, dim = dim(df))

  nas <- is.na(df)
  nums <- vapply(df, is.numeric, logical(1))
  infs[,nums] <- vapply(which(nums), function(i) !nas[,i] & !is.finite(df[[i]]), logical(NROW(df)))

  for (i in which(!nums)) {
    if (can_str2num(df[[i]])) {
      df[[i]] <- str2num(df[[i]])
      nums[i] <- TRUE
    }
  }

  # o.negs[,nums] <- !nas[,nums] & df[nums] < 0 & round(df[nums], digits) == 0
  df[nums] <- round(df[nums], digits = digits)

  for (i in which(nums)) {
    df[[i]] <- format(df[[i]], scientific = FALSE, justify = "none", trim = TRUE,
                      drop0trailing = !identical(as.character(pad), "0"))

    if (!identical(as.character(pad), "0") && any(grepl(".", df[[i]], fixed = TRUE))) {
      s <- strsplit(df[[i]], ".", fixed = TRUE)
      lengths <- lengths(s)
      digits.r.of.. <- rep(0, NROW(df))
      digits.r.of..[lengths > 1] <- nchar(vapply(s[lengths > 1], `[[`, character(1L), 2))
      max.dig <- max(digits.r.of..)

      dots <- ifelse(lengths > 1, "", if (as.character(pad) != "") "." else pad)
      pads <- vapply(max.dig - digits.r.of.., function(n) paste(rep(pad, n), collapse = ""), character(1L))

      df[[i]] <- paste0(df[[i]], dots, pads)
    }
  }

  # df[o.negs] <- paste0("-", df[o.negs]) #Requested to remove to prevent -0

  # Insert NA placeholders
  df[nas] <- na_vals
  df[infs] <- "N/A"

  if (length(rn) > 0) rownames(df) <- rn
  if (length(cn) > 0) names(df) <- cn

  attr(df, "na_vals") <- na_vals
  df
}

#Adds perentheses around a number in SD columns; e.g., 5.46 -> (5.46)
add_peren_to_sd <- function(df) {
  for (i in names(df)) {
    if (startsWith(i, "SD") && !all(df[[i]] == attr(df, "na_vals"))) {
      df[[i]][df[[i]] != attr(df, "na_vals")] <- sprintf("(%s)",  df[[i]][df[[i]] != attr(df, "na_vals")])
    }
  }
  df
}

#Transform number to subscript
num2sub <- function(x) {
  x <- as.character(x)

  chartr("0123456789",
         "\u2080\u2081\u2082\u2083\u2084\u2085\u2086\u2087\u2088\u2089",
         x)
}

#Get covariates from data frame; for use in summary()
covs_df_to_matrix <- function(covs) {

  if (NCOL(covs) == 0) {
    return(as.matrix(covs))
  }
  fnames <- colnames(covs)
  fnames[!startsWith(fnames, "`")] <- paste0("`", fnames[!startsWith(fnames, "`")], "`")
  formula <- reformulate(fnames)

  mf <- model.frame(terms(formula, data = covs), covs,
                    na.action = na.pass)

  chars.in.mf <- vapply(mf, is.character, logical(1L))
  mf[chars.in.mf] <- lapply(mf[chars.in.mf], factor)

  X <- model.matrix(formula, data = mf,
                    contrasts.arg = lapply(Filter(is.factor, mf),
                                           contrasts, contrasts = FALSE))
  assign <- attr(X, "assign")[-1]
  X <- X[,-1, drop=FALSE]
  attr(X, "assign") <- assign

  X
}

#Quickly compute diagonal of hat matrix without having to compute
#full projection matrix. Uses a special formula when a fixed effect (f)
#is present to simplify calculation. Assumes X first column is an
#intercept.
hat_fast <- function(X, w = NULL, f = NULL) {
  if (is.null(f)) {
    QR <- {
      if (is.null(w)) qr.default(X)
      else qr.default(sqrt(w) * X)
    }

    Q <- qr.qy(QR, diag(1, nrow = nrow(QR$qr), ncol = QR$rank))

    return(rowSums(Q * Q))
  }

  #Fixed effects block version
  fmm <- do.call("cbind", lapply(levels(f), function(i) as.numeric(f == i)))

  if (!is.null(w)) {
    rw <- sqrt(w)
    diag_h_f <- hat_fast(fmm, w)
  }
  else {
    rw <- 1
    diag_h_f <- 1/tabulate(f)[as.integer(f)]
  }

  diag_h_X <- hat_fast(.lm.fit(rw*fmm, rw*X[,-1, drop = FALSE])$residuals/rw, w)

  diag_h_f + diag_h_X
}

treat_name_from_coefs <- function(coef_names, treat_levels) {
  shortest_name <- coef_names[which.min(nchar(coef_names))]
  for (i in 1:nchar(shortest_name)) {
    treat <- substring(shortest_name, 1, i)
    if (sum(paste0(treat, treat_levels) %in% coef_names) == length(coef_names)) {
      return(treat)
    }
  }

  return("")
}

treat_levels_from_coefs <- function(coef_names, treat_levels, treat_name = NULL) {
  if (is.null(treat_name)) {
    treat_name <- treat_name_from_coefs(coef_names, treat_levels)
  }

  coef_levels <- sub(treat_name, "", coef_names, fixed = TRUE)

  c(setdiff(treat_levels, coef_levels), coef_levels)
}

#Group mean centers a variable x for a factor f. For
#use with fixed effects.
demean <- function(x, f, w = NULL) {
  f <- as.factor(f)
  for (i in levels(f)) {
    x[f == i] <- x[f == i] - mean_w(x, w, subset = f == i)
  }
  x
}

#Check whether function is being called within another specific function
called_from <- function(...) {
  calls <- sys.calls()
  any(unlist(list(...)) %in% unlist(lapply(calls, function(x) deparse1(x[[1]]))))
}
