#More informative and cleaner version of base::match.arg. From WeightIt.
match_arg <- function(arg, choices, several.ok = FALSE) {
  #Replaces match.arg() but gives cleaner error message and processing
  #of arg.
  if (missing(arg))
    stop("No argument was supplied to match_arg.", call. = FALSE)
  arg.name <- paste(deparse(substitute(arg), width.cutoff = 500L), collapse = " ")

  if (missing(choices)) {
    formal.args <- formals(sys.function(sysP <- sys.parent()))
    choices <- eval(formal.args[[as.character(substitute(arg))]],
                    envir = sys.frame(sysP))
  }

  if (is.null(arg))
    return(choices[1L])
  else if (!is.character(arg))
    stop(paste0("The argument to '", arg.name, "' must be NULL or a character vector"), call. = FALSE)
  if (!several.ok) {
    if (identical(arg, choices))
      return(arg[1L])
    if (length(arg) > 1L)
      stop(paste0("The argument to '", arg.name, "' must be of length 1"), call. = FALSE)
  }
  else if (length(arg) == 0)
    stop(paste0("The argument to '", arg.name, "' must be of length >= 1"), call. = FALSE)

  i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
  if (all(i == 0L))
    stop(paste0("The argument to '", arg.name, "' should be ", if (length(choices) > 1) {if (several.ok) "at least one of " else "one of "} else "",
                word_list(choices, and.or = "or", quotes = 2), "."),
         call. = FALSE)
  i <- i[i > 0L]
  if (!several.ok && length(i) > 1)
    stop("There is more than one match in 'match_arg'")
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
      and.or <- match_arg(and.or)
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
  return(out)
}

#Add quotation marks around a string.
add_quotes <- function(x, quotes = 2L) {
  if (!isFALSE(quotes)) {
    if (isTRUE(quotes) || as.integer(quotes) == 2L) x <- paste0("\"", x, "\"")
    else if (as.integer(quotes) == 1L) x <- paste0("\'", x, "\'")
    else stop("'quotes' must be boolean, 1, or 2.")
  }
  x
}

#Effective sample size (using |w|)
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
  else return(colSums(mat * w)/sum(w))
}

#Weighted mean
mean_w <- function(x, w = NULL, subset = NULL) {
  if (length(subset) != 0) {
    x <- x[subset]
    if (length(w) != 0) w <- w[subset]
  }

  if (length(w) == 0) return(sum(x)/length(x))
  else return(sum(x * w)/sum(w))
}

#Determine whether a character vector can be coerced to numeric
can_str2num <- function(x) {
  nas <- is.na(x)
  suppressWarnings(x_num <- as.numeric(as.character(x[!nas])))
  return(!anyNA(x_num))
}

#Cleanly coerces a character vector to numeric; best to use after can_str2num()
str2num <- function(x) {
  nas <- is.na(x)
  suppressWarnings(x_num <- as.numeric(as.character(x)))
  x_num[nas] <- NA
  return(x_num)
}

#Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

#Clean printing of data frames with numeric and NA elements.
round_df_char <- function(df, digits, pad = "0", na_vals = "") {
  #Digits is passed to round(). pad is used to replace trailing zeros so decimal
  #lines up. Should be "0" or " "; "" (the empty string) un-aligns decimals.
  #na_vals is what NA should print as.

  if (NROW(df) == 0 || NCOL(df) == 0) return(df)
  if (!is.data.frame(df)) df <- as.data.frame.matrix(df, stringsAsFactors = FALSE)
  rn <- rownames(df)
  cn <- colnames(df)

  infs <- o.negs <- array(FALSE, dim = dim(df))
  nas <- is.na(df)
  nums <- vapply(df, is.numeric, logical(1))
  infs[,nums] <- vapply(which(nums), function(i) !nas[,i] & !is.finite(df[[i]]), logical(NROW(df)))

  for (i in which(!nums)) {
    if (can_str2num(df[[i]])) {
      df[[i]] <- str2num(df[[i]])
      nums[i] <- TRUE
    }
  }

  o.negs[,nums] <- !nas[,nums] & df[nums] < 0 & round(df[nums], digits) == 0
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

  df[o.negs] <- paste0("-", df[o.negs])

  # Insert NA placeholders
  df[nas] <- na_vals
  df[infs] <- "N/A"

  if (length(rn) > 0) rownames(df) <- rn
  if (length(cn) > 0) names(df) <- cn

  return(df)
}

#Transform number to subscript
num2sub <- function(x) {
  x <- as.character(x)

  chartr("0123456789",
         "\u2080\u2081\u2082\u2083\u2084\u2085\u2086\u2087\u2088\u2089",
         x)
}

#Get covariates (RHS) vars from formula
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

  return(X)
}

#Turn a vector into a 0/1 vector. 'zero' and 'one' can be supplied to make it clear which is
#which; otherwise, a guess is used. From WeightIt.
binarize <- function(variable, zero = NULL, one = NULL) {
  if (length(unique(variable)) > 2) stop(paste0("Cannot binarize ", deparse1(substitute(variable)), ": more than two levels."))
  if (is.character(variable) || is.factor(variable)) {
    variable <- factor(variable, nmax = 2)
    unique.vals <- levels(variable)
  }
  else {
    unique.vals <- unique(variable, nmax = 2)
  }

  if (is.null(zero)) {
    if (is.null(one)) {
      if (can_str2num(unique.vals)) {
        variable.numeric <- str2num(variable)
      }
      else {
        variable.numeric <- as.numeric(variable)
      }

      if (0 %in% variable.numeric) zero <- 0
      else zero <- min(variable.numeric, na.rm = TRUE)

      return(setNames(as.integer(variable.numeric != zero), names(variable)))
    }
    else {
      if (one %in% unique.vals) return(setNames(as.integer(variable == one), names(variable)))
      else stop("The argument to 'one' is not the name of a level of variable.", call. = FALSE)
    }
  }
  else {
    if (zero %in% unique.vals) return(setNames(as.integer(variable != zero), names(variable)))
    else stop("The argument to 'zero' is not the name of a level of variable.", call. = FALSE)
  }
}

hat_fast <- function(X, w = NULL) {
  if (is.null(w)) QR <- qr.default(X)
  else QR <- qr.default(sqrt(w) * X)

  Q <- qr.qy(QR, diag(1, nrow = nrow(QR$qr), ncol = QR$rank))
  rowSums(Q * Q)
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
  return(c(setdiff(treat_levels, coef_levels), coef_levels))
}
