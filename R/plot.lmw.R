plot.lmw <- function(x, type = "weights", ...) {
  type <- match_arg(type, c("weights", "extrapolation", "influence"))

  if (type == "weights") {
    weights_plot(x, ...)
  }
  else if (type == "extrapolation") {
    extrapolation_plot(x, ...)
  }
  else if (type == "influence") {
    influence_plot(x, ...)
  }

  invisible(x)
}

weights_plot <- function(x, rug = TRUE, mean = TRUE, ess = TRUE, ...) {
  .pardefault <- par(no.readonly = TRUE)
  on.exit(par(.pardefault))

  t <- x$treat
  w <- x$weights
  sw <- if (is.null(x$s.weights)) rep(1, length(t)) else x$s.weights

  tlevs <- if (is.null(x$contrast)) levels(t) else x$contrast

  layout(do.call("rbind", as.list(seq_along(tlevs))))

  par(mar = c(2.75, 3, 1.75, 1),
      mgp = c(1.5, 0.5, 0))

  for (i in tlevs) {
    if (length(tlevs) == 2 && i != tlevs[1] && x$method == "URI") {
      in_i <- which(t != tlevs[1])
    }
    else {
      in_i <- which(t == i)
    }

    wi <- w[in_i]

    n <- length(wi)

    mw <- sum(wi)/n

    if (ess) {
      e <- ESS(wi)
      e_un <- ESS(sw[t == i])
    }

    wi_range <- range(wi)
    if (abs(diff(wi_range)) < sqrt(.Machine$double.eps)) {
      withCallingHandlers(
        dens <- density(wi, adjust = 1e-4, ...),
        warning = function(w){
          if (grepl("collapsing to unique 'x' values", w$message))
            invokeRestart("muffleWarning")
        })
      # Note: without the above lines, a warning is thrown; see
      # https://stackoverflow.com/a/56861070/6348551
      xlim <- c(0, 2*wi_range[2])
    }
    else {
      dens <- density(wi, ...)
      x_range <- range(dens$x)
      xlim <- {
        if (ess) c(x_range[1],
                   x_range[2] + .1 * diff(x_range))
        else x_range
      }
    }

    ylim <- c(0, max(dens$y)*1.05)

    plot(dens,
         ylab = "Density",
         xlab = "Weight",
         yaxt = "n",
         cex.lab = 0.8,
         cex.axis = 0.8,
         xlim = xlim,
         ylim = ylim,
         yaxs = "i",
         main = "")
    title(main = sprintf("Distributon of Weights (%s)", {
      if (length(tlevs) == 2 && all(tlevs %in% c("0", "1"))) switch(i, "0" = "Control", "1" = "Treated")
      else i
    }),
    cex.main = 0.9,
    font.main = 1,
    line = 0.5)
    if (rug) {
      rug(wi,
          col = "black",
          ticksize = 0.1,
          lwd = 0.3,
          side = 1)
    }
    if (mean) {
      rug(mw, col = "red", ticksize = 0.3, lwd = 1)
    }

    if (ess) {
      legend("topright",
             legend = sprintf("N = %s\nESS = %s", round(e_un, 1), round(e, 1)),
             bty = "n",
             cex = 0.8)
    }
  }
}

extrapolation_plot <- function(x, var, data = NULL, ...) {
  .pardefault <- par(no.readonly = TRUE)
  on.exit(par(.pardefault))

  par(mar = c(2.75, 3, 1.75, 1),
      mgp = c(1.5, 0.5, 0))

  # set.seed(1234)

  t <- x$treat
  w <- x$weights

  tlevs <- if (is.null(x$contrast)) levels(t) else x$contrast

  data <- get_data(data, x)

  if (missing(var)) {
    stop("The variables for which extrapolation is to be assessed must be named in the 'var' argument.", call. = FALSE)
  }
  if (is.character(var)) {
    if (!is.null(data) && is.data.frame(data)) {
      if (all(var %in% names(data))) {
        covs <- covs_df_to_matrix(data[var])
      }
      else {
        stop("All variables in 'var' must be in 'data'.", call. = FALSE)
      }
    }
    else {
      stop("If 'var' is specified as a string, a data frame argument must be supplied to 'data'.", call. = FALSE)
    }
  }
  else if (inherits(var, "formula")) {
    vars.in.formula <- all.vars(var)
    if (!is.null(data) && is.data.frame(data)) data <- cbind(data[names(data) %in% vars.in.formula],
                                                             x$covs[names(data) %in% setdiff(vars.in.formula, names(data))])
    else data <- x$covs

    covs <- covs_df_to_matrix(model.frame(var, data = data))
  }

  v <- as.data.frame(covs)

  if (!is.null(x$target)) {
    X_target <- covs_df_to_matrix(model.frame(remove_treat_from_formula(delete.response(terms(x$formula)), attr(x$treat, "treat_name")),
                                              data = x$target))
    target_means <- colMeans_w(X_target, attr(x$target, "target.weights"))

    if (length(addl_diff <- setdiff(colnames(covs), colnames(X_target))) > 0) {
      target_means <- c(target_means, setNames(rep(NA_real_, length(addl_diff)), addl_diff))
    }
  }

  par(mfrow = c(1, ncol(v)))

  cex <- numeric(length(w))
  for (i in tlevs) {
    if (length(tlevs) == 2 && i != tlevs[1] && x$method == "URI") {
      in_i <- which(t != tlevs[1])
    }
    else {
      in_i <- which(t == i)
    }
    # sqrt(): makes it so size of weigth corresponds to area, not radius
    # length(in_i)*: makes it so
    # cex[in_i] <- sqrt(sum(abs(w[in_i]) > 1e-8)*abs(w[in_i] / sum(w[in_i]))) #Total weights sum to n on each side
    cex[in_i] <- 30*sqrt(abs(w[in_i])/sum(w[in_i])) #Total weights sum to 1 on each side, negative weights count as negative
    # cex[in_i] <- 30*sqrt(abs(w[in_i])/sum(abs(w[in_i]))) #Total abs(weights) sum to 1 on each side, same amount of ink on both sides
  }

  if (length(v) <= 2) cex.text <- .8
  else cex.text <- .8/.66

  #Transparency of points depends on number of points logarithmically
  alpha <- vapply(tlevs, function(i) min(1, .6/log10(sum(abs(w[t == i]) > 1e-8))),
                  numeric(1L))

  col <- character(length(t))
  for (i in tlevs) {
    col[t == i & w >= 0] <- adjustcolor("black", alpha.f = alpha[i])
    col[t == i & w < 0] <- adjustcolor("red", alpha.f = alpha[i])
  }

  for (j in names(v)) {
    vj <- v[[j]]

    means_vj <- vapply(tlevs, function(i) {
      if (length(tlevs) == 2 && i != tlevs[1] && x$method == "URI") {
        in_i <- which(t != tlevs[1])
      }
      else {
        in_i <- which(t == i)
      }
      mean_w(vj, w, in_i)
    }, numeric(1L))

    target_vj <- {
      if (!is.null(x$focal)) mean_w(vj, x$s.weights, t==x$focal)
      else if (!is.null(x$target)) target_means[j]
      else mean_w(vj, x$s.weights)
    }

    if (length(unique(vj)) <= 2) vj <- vj + runif(length(vj), -.02, .02)

    range_vj <- range(vj)
    xlim <- c(min(range(vj)[1], target_vj, means_vj) - .025 * diff(range_vj),
              max(range(vj)[2], target_vj, means_vj) + .025 * diff(range_vj))

    jitter <- runif(length(t), -.1, .1)

    #Data points
    plot(x = vj,
         y = .5 + 2*(length(tlevs) - as.numeric(t)) + 1*(w >= 0) + jitter,
         pch = 20,
         col = col,
         ylim = c(0, 2*length(tlevs)),
         xlim = xlim,
         xlab = j,
         yaxt = "n",
         yaxs = "i",
         cex = cex,
         ylab = NA,
         main = NA,
         frame.plot = TRUE,
         cex.axis = cex.text,
         cex.lab = cex.text)

    #Lines seperating treatment groups
    abline(h = 2 * seq_along(tlevs)[-length(tlevs)])

    #Vertical lines for observed weighted means
    for (i in seq_along(tlevs)) {
      segments(x0 = means_vj[i],
               y0 = 2*(length(tlevs) - i + 1),
               y1 = 2*(length(tlevs) - i))
    }

    #Xs for target means
    points(x = rep(target_vj, length(tlevs)),
           y = 2*seq_along(tlevs) - 1,
           pch = 4, cex = 4, lwd = 0.5)

    axis(2, at = 2*rev(seq_along(tlevs)) - 1,
         labels = {
           if (length(tlevs) == 2 && all(tlevs %in% c("0", "1"))) {
             c("0" = "Control", "1" = "Treated")[tlevs]
           }
           else tlevs
         },
         tick = FALSE, cex.axis = cex.text)
  }

}

influence_plot <- function(x, outcome, data = NULL, id.n = 3, ...) {
  call <- match.call(expand.dots = FALSE)

  .pardefault <- par(no.readonly = TRUE)
  on.exit(par(.pardefault))

  par(mar = c(5.1, 4.1, 1.75, 1),
      mgp = c(3, 1, 0))

  inf <- do.call("influence", list(x, substitute(outcome), data))

  SIC <- inf$sic

  ## Scaled SIC
  SIC_std <- SIC/max(SIC)

  plot(SIC_std,
       type = c("h"),
       lty = c("solid"),
       xlab = "Index",
       ylab = "Scaled SIC",
       cex.lab = 0.8,
       cex.axis = 0.8,
       col = grey(0),
       axes = FALSE)

  labels.id <- seq_along(SIC)

  show.r <- order(SIC_std, decreasing = TRUE)[1:id.n]
  y.id <- show.r
  y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
  labpos <- c(4, 2)[(y.id > mean(par("usr")[1:2])) + 1L]
  text(y.id, SIC_std[show.r], labels.id[show.r], cex = 0.75, xpd = TRUE,
       pos = labpos, offset = 0.25)

  N <- length(SIC_std)
  indices <- c(1, pretty(seq_len(N))[-1])
  indices <- indices[-length(indices)]
  indices[length(indices)] <- N
  axis(side = 1,
       at = indices,
       cex.axis = .8)
  axis(side = 2,
       at = seq(0, 1, 0.2),
       cex.axis = .8)
}

plot.lmw_est <- function(x, type = "influence", ...) {
  type <- match_arg(type, c("influence", "lm"))

  if (type == "influence") {
    influence_plot(x, ...)
  }
  else if (type == "lm") {
    class(x) <- c(class(x), "lm")
    utils::getS3method("plot", "lm")(x, ...)
  }

  invisible(x)
}
