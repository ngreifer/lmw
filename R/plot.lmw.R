#' Plots diagnosing regression-implied weights
#'
#' @description
#' Produces plots to diagnose properties of the weights, including their
#' distribution, to what degree the distribution of covariates involves
#' extrapolation in the weighted sample, and how much influence each unit has
#' on the effect estimate.
#'
#' @details
#' When `type = "weights"`, `plot.lmw()` produces a density plot of
#' the weights within each treatment group. By construction, these weights will
#' have a mean of 1. Some weights may be negative. The effective sample size
#' (ESS) and original sample size (N) will be displayed in the upper right
#' corner of the plot when `ess = TRUE`.
#'
#' When `type = "extrapolation"`, `plot.lmw()` produces a plot of the
#' distribution of weights and covariates for each treatment group. Each dot
#' represents a unit, with values arranged on the x-axis according to their
#' covariate value and the size of the dots corresponding to the magnitude of
#' the weight. Units with positive weights are displayed in black in the upper
#' portion of the plot, and units with negative weights are displayed in red in
#' the lower portion. Having many and large red points indicates a high degree
#' of extrapolation. All points are equally transparent, so darker regions
#' indicate multiple points with the same value. The vertical lines indicates
#' the weighted mean of the covariate in each group, and the X indicates the
#' mean of the covariate in the target sample as determined by the
#' `estimand` argument in the original call to `lmw()`. A large
#' discrepancy between the vertical lines and Xs indicates a lack of balance
#' between the treatment group and target sample. When `estimand = "CATE"`
#' in the original call to `lmw()`, any variables supplied to `variables`
#' that were not given a target value will not have the target mean displayed.
#'
#' When `type = "influence"`, `plot.lmw()` produces a plot of the
#' scaled sample influence curve (SIC) for each unit by index. It does so by
#' calling [influence.lmw()], which fits the outcome model to extract
#' residuals and compute the SIC as `SIC = (N-1) * w * r / (1 - h)`, where
#' `N` is the sample size, `w` are the units' implied regression
#' weights, `r` are the residuals, and `h` are the hat values. SIC
#' values are scaled to have a maximum of 1. Higher values indicate greater
#' relative influence.
#'
#' @param x an `lmw` object; the output of a call to [lmw()].
#' @param type the type of plot to display. Allowable options include
#' `"weights"`, `"extrapolation"`, and `"influence"`. See
#' Details. Abbreviations allowed.
#' @param \dots further arguments passed to specific types of plots.
#'
#' When `type = "weights"`, the following are accepted:
#' \describe{
#' \item{`rug`}{`logical`; whether to display a rug plot of the
#' weights. Default is `TRUE`.}
#' \item{`mean`}{whether to display a
#' red line indicating the mean of the weights. Default is `TRUE`.}
#' \item{`ess`}{whether to display the original and weighted effective
#' sample size in the top right corner. Default is `TRUE`.}
#' }
#' Other arguments are passed to [density()].
#'
#' When `type = "extrapolation"`, the following are accepted:
#' \describe{
#' \item{`variables`}{required; a right-sided formula or character vector
#' containing the names of the covariates for which extrapolation is to be
#' assessed.}
#' \item{`data`}{an optional data frame containing the
#' variables named in `variables`.}
#' }
#'
#' When `type = "influence"`, the
#' following are accepted:
#' \describe{
#' \item{`outcome`}{the name of the
#' outcome variable. Can be supplied as a string containing the name of the
#' outcome variable or as the outcome variable itself. If not supplied, the
#' outcome variable in the `formula` supplied to `lmw()`, if any,
#' will be used.}
#' \item{`data`}{an optional data frame containing the
#' outcome variable named in `outcome`.}
#' \item{`id.n`}{the number of
#' points to be labelled in the plot, starting with the most extreme.}
#' }
#'
#' @return A plot is displayed, and `x` is invisibly returned.
#'
#' @seealso [lmw()], [summary.lmw()],
#' [plot.summary.lmw()]
#'
#' @examples
#' data("lalonde")
#'
#' # URI regression for ATT
#' lmw.out1 <- lmw(~ treat + age + education + race + married +
#'                    nodegree + re74 + re75, data = lalonde,
#'                 estimand = "ATT", method = "URI",
#'                 treat = "treat")
#' lmw.out1
#'
#' # Distribution of weights
#' plot(lmw.out1, type = "weights")
#'
#' # Extrapolation/representativeness for age and married
#' plot(lmw.out1, type = "extrapolation",
#'      variables = ~age + married)
#'
#' # Extrapolation/representativeness for race
#' plot(lmw.out1, type = "extrapolation",
#'      variables = ~race)
#'
#' # Influence for re78 outcome
#' plot(lmw.out1, type = "influence", outcome = "re78")

#' @exportS3Method plot lmw
plot.lmw <- function(x, type = "weights", ...) {
  chk::chk_string(type)
  type <- tolower(type)
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

  grDevices::dev.hold()
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
  grDevices::dev.flush()
}

extrapolation_plot <- function(x, variables, data = NULL, ...) {
  .pardefault <- par(no.readonly = TRUE)
  on.exit(par(.pardefault))

  par(mar = c(2.75, 3, 1.75, 1),
      mgp = c(1.5, 0.5, 0))

  # set.seed(1234)

  t <- x$treat
  w <- x$weights

  tlevs <- if (is.null(x$contrast)) levels(t) else x$contrast

  data <- get_data(data, x)

  if (missing(variables)) {
    chk::err("the variables for which extrapolation is to be assessed must be named in the `variables` argument")
  }
  if (is.character(variables)) {
    if (is.null(data) || !is.data.frame(data)) {
      chk::err("if `variables` is specified as a string, a data frame argument must be supplied to `data`")
    }

    if (!all(variables %in% names(data))) {
      chk::err("all variables in `variables` must be in `data`")
    }

    covs <- covs_df_to_matrix(data[variables])
  }
  else if (inherits(variables, "formula")) {
    vars.in.formula <- all.vars(variables)
    if (!is.null(data) && is.data.frame(data)) {
      data <- cbind(data[names(data) %in% vars.in.formula],
                    x$covs[names(data) %in% setdiff(vars.in.formula, names(data))])
    }
    else data <- x$covs

    covs <- covs_df_to_matrix(model.frame(variables, data = data))
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

  cex.text <- if (length(v) <= 2) .8  else .8/.66

  #Transparency of points depends on number of points logarithmically
  alpha <- vapply(tlevs, function(i) min(1, .6/log10(sum(abs(w[t == i]) > 1e-8))),
                  numeric(1L))

  col <- character(length(t))
  for (i in tlevs) {
    col[t == i & w >= 0] <- grDevices::adjustcolor("black", alpha.f = alpha[i])
    col[t == i & w < 0] <- grDevices::adjustcolor("red", alpha.f = alpha[i])
  }

  grDevices::dev.hold()
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
  grDevices::dev.flush()

}

influence_plot <- function(x, outcome, data = NULL, id.n = 3, ...) {
  .pardefault <- par(no.readonly = TRUE)
  on.exit(par(.pardefault))

  # par(mar = c(5.1, 4.1, 1.75, 1),
  #     mgp = c(3, 1, 0))

  inf <- do.call("influence", list(x, substitute(outcome), data))

  SIC <- inf$sic

  ## Scaled SIC
  SIC_std <- SIC/max(SIC)

  labels.id <- seq_along(SIC)

  N <- length(SIC_std)
  indices <- c(1, pretty(seq_len(N))[-1])
  # indices <- indices[-length(indices)]
  # indices[length(indices)] <- N

  grDevices::dev.hold()
  plot(SIC_std,
       type = "h",
       lty = "solid",
       xlab = "Obs. number",
       ylab = "Scaled SIC",
       ylim = c(0, 1.075),
       # cex.lab = 0.8,
       # cex.axis = 0.8,
       col = grDevices::gray(0),
       xaxt = "n",
       yaxt = "n")

  show.r <- order(SIC_std, decreasing = TRUE)[1:id.n]
  y.id <- show.r
  y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
  labpos <- c(4, 2)[(y.id > mean(par("usr")[1:2])) + 1L]

  text(y.id, SIC_std[show.r], labels.id[show.r], cex = 0.75, xpd = TRUE,
       pos = labpos, offset = 0.25)

  axis(side = 1,
       at = indices)
  axis(side = 2,
       at = seq(0, 1, 0.25))
  mtext(grDevices::as.graphicsAnnot("Sample Influence Curve"), 3, 0.25, cex = 1)
  grDevices::dev.flush()
}
