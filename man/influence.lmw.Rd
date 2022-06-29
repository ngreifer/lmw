\name{influence.lmw}
\alias{influence.lmw}
\alias{influence.lmw_est}

\title{
Regression Diagnostics for \code{lmw} and \code{lmw_est} objects
}
\description{
\code{influence()} produces influence measures for \code{lmw} objects that can be used as regression diagnostics to identify influential cases. These functions produce similar outputs to \code{\link{lm.influence}} but also include the sample influence curve (SIC) values, which combine information about the hat values, residuals, and implied regression weights.
}
\usage{
\method{influence}{lmw}(model, outcome, data = NULL, ...)

\method{influence}{lmw_est}(model, ...)
}

\arguments{
  \item{model}{
an \code{lmw} or \code{lmw_est} object; the output of a call to \code{\link{lmw}} or \code{\link{lmw_est}}.
}
  \item{outcome}{
the name of the outcome variable. Can be supplied as a string containing the name of the outcome variable or as the outcome variable itself. If not supplied, the outcome variable in the \code{formula} supplied to \code{lmw()}, if any, will be used.
}
  \item{data}{
an optional data frame containing the outcome variable named in \code{outcome}.
}
  \item{\dots}{
ignored.
}
}
\details{
\code{influence()} computes the hat values, (weighted) residuals, and sample influence curve (SIC) values for each unit, which can be used as regression diagnostics to assess influence. The weighted residuals are weighted by the sampling weights (if supplied), not the implied regression weights. The SIC values are computed as \code{SIC = (N-1) * w * r / (1 - h)}, where \code{N} is the sample size, \code{w} are the units' implied regression weights, \code{r} are the (weighted) residuals, and \code{h} are the hat values. SIC values are scaled to have a maximum of 1. Higher values indicate greater relative influence.
}
\value{
A list with the following components:
\item{hat}{a vector containing the diagonal of the hat matrix.}
\item{wt.res}{a vector of (weighted) residuals.}
\item{sic}{a vector containing the scaled SIC values.}
}

\note{
\code{influence.lmw()} uses non-standard evaluation to interpret its \code{outcome} argument. For programmers who wish to use \code{influence.lmw()} inside other functions, an effective way to pass the name of an arbitrary outcome (e.g., \code{y} passed as a string) is to use \code{\link{do.call}}, for example:
\preformatted{fun <- function(m, y, d) {
   do.call("influence", list(m, y, d))
}
}
When using \code{influence.lmw()} inside \code{\link{lapply}} or \code{purrr::map} to loop over outcomes, this syntax must be used as well.
}

\seealso{
\code{\link{plot.lmw}} for plotting the SIC values; \code{\link{lm.influence}} for influence measures for \code{lm} objects, which do not include SIC values; \code{\link{hatvalues}} for hat values for \code{lm} objects (note that \code{lmw_est} objects also have a \code{hatvalues()} method).
}
\examples{
data("lalonde")

#URI regression for ATT
lmw.out1 <- lmw(~ treat + age + educ + race + married +
                     nodegree + re74 + re75,
                data = lalonde, estimand = "ATT",
                method = "URI", treat = "treat")

#Influence for re78 outcome
infl <- influence(lmw.out1, outcome = "re78")
str(infl)

#Can also be used after lmw_est():
lmw.est1 <- lmw_est(lmw.out1, outcome = "re78")
all.equal(infl,
          influence(lmw.est1))
}
