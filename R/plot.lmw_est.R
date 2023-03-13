#' Plot diagnostics for an `lmw_est` object
#'
#' @description
#' Produces plots to diagnose the regression model fit to estimate the
#' treatment effect. These include an influence plot based on the sample
#' influence curve (SIC) and the regression diagnostics plots available for
#' `lm` objects in [plot.lm()].
#'
#' @details
#' When `type = "influence"`, `plot.lmw_est()` produces a plot of the
#' scaled sample influence curve (SIC) for each unit by index. It does so by
#' calling [influence.lmw_est()], which extract the model residuals
#' and computes the SIC as `SIC = (N-1) * w * r / (1 - h)`, where `N`
#' is the sample size, `w` are the units' implied regression weights,
#' `r` are the residuals, and `h` are the hat values. SIC values are
#' scaled to have a maximum of 1. Higher values indicate greater relative
#' influence.
#'
#' When `type = "lm"`, `plot.lmw_est()` produces several plots
#' displayed sequentially according to the arguments supplied to `plot()`.
#' These plots are produced by [plot.lm()] to diagnose the
#' distribution of residuals and other measures of leverage and influence.
#'
#' @param x an `lmw_est` object; the output of a call to
#' [lmw_est()].
#' @param type the type of plot to display. Allowable options include
#' `"influence"` and `"lm"`. See Details. Abbreviations allowed.
#' @param \dots When `type = "influence"`, the following are accepted:
#' \describe{
#' \item{`outcome`}{the name of the outcome variable. Can be
#' supplied as a string containing the name of the outcome variable or as the
#' outcome variable itself. If not supplied, the outcome variable in the
#' `formula` supplied to `lmw()`, if any, will be used.}
#' \item{`data`}{an optional data frame containing the outcome variable
#' named in `outcome`.}
#' \item{`id.n`}{the number of points to be
#' labelled in the plot, starting with the most extreme.}
#' }
#'
#' When `type =
#' "lm"`, any arguments passed to [plot.lm()] are accepted and passed
#' directly to `plot.lm`.
#'
#' @return A plot is displayed, and `x` is invisibly returned.
#'
#' @seealso [lmw_est()], [influence.lmw_est()],
#' [plot.lm()]
#'
#' @examples
#' data("lalonde")
#'
#' # URI regression for ATT
#' lmw.out1 <- lmw(~ treat + age + education + race + married +
#'                    nodegree + re74 + re75, data = lalonde,
#'                 estimand = "ATT", method = "URI",
#'                 treat = "treat")
#'
#' lmw.fit1 <- lmw_est(lmw.out1, outcome = "re78")
#' lmw.fit1
#'
#' # Influence using SIC
#' plot(lmw.fit1, type = "influence")
#'
#' # Usual regression diagnostics
#' plot(lmw.fit1, type = "lm", which = 1)

#' @exportS3Method plot lmw_est
plot.lmw_est <- function(x, type = "influence", ...) {
  type <- match_arg(type, c("influence", "lm"))

  if (type == "influence") {
    influence_plot(x, ...)
  }
  else if (type == "lm") {
    # Need to assign `"lm"` as last class and explicitly call
    # plot.lm() to get lm plots to work.
    cl <- class(x)
    class(x) <- c(cl, "lm")
    plot.lm(x, ...)
    class(x) <- cl
  }

  invisible(x)
}

# To be used internally
plot.lm <- utils::getS3method("plot", "lm")
