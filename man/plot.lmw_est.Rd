\name{plot.lmw_est}
\alias{plot.lmw_est}
\title{
Plot Diagnnostics for an \code{lmw_est} Object
}
\description{
Produces plots to diagnose the regression model fit to estimate the treatment effect. These include an influence plot based on the sample influence curve (SIC) and the regression diagnostics plots available for \code{lm} objects in \code{\link{plot.lm}}.
}
\usage{
\methd{plot}{lmw_est}(x, type = "influence", ...)
}

\arguments{
  \item{x}{
an \code{lmw_est} object; the output of a call to \code{\link{lmw_est}}.
}
  \item{type}{
the type of plot to display. Allowable options include \code{"influence"} and \code{"lm"}. See Details. Abbreviations allowed.
}
  \item{\dots}{
When \code{type = "influence"}, the following are accepted:
\describe{
  \item{\code{outcome}}{the name of the outcome variable. Can be supplied as a string containing the name of the outcome variable or as the outcome variable itself. If not supplied, the outcome variable in the \code{formula} supplied to \code{lmw()}, if any, will be used.}
  \item{\code{data}}{an optional data frame containing the outcome variable named in \code{outcome}.}
  \item{\code{id.n}}{the number of points to be labelled in the plot, starting with the most extreme.}
}
When \code{type = "lm"}, any arguments passed to \code{\link{plot.lm}} are accepted and passed directly to \code{plot.lm}.
}
}
\details{
When \code{type = "influence"}, \code{plot.lmw_est()} produces a plot of the scaled sample influence curve (SIC) for each unit by index. It does so by calling \code{\link{influence.lmw_est}}, which extract the model residuals and computes the SIC as \code{SIC = (N-1) * w * r / (1 - h)}, where \code{N} is the sample size, \code{w} are the units' implied regression weights, \code{r} are the residuals, and \code{h} are the hat values. SIC values are scaled to have a maximum of 1. Higher values indicate greater relative influence.

When \code{type = "lm"}, \code{plot.lmw_est()} produces several plots displayed sequentially according to the arguments supplied to \code{plot()}. These plots are produced by \code{\link{plot.lm}} to diagnose the distribution of residuals and other measures of leverage and influence.
}
\value{
A plot is displayed, and \code{x} is invisibly returned.
}
\seealso{
\code{\link{lmw_est}}, \code{\link{influence.lmw_est}}, \code{\link{plot.lm}}
}
\examples{
data("lalonde")

#URI regression for ATT
lmw.out1 <- lmw(~ treat + age + educ + race + married +
                   nodegree + re74 + re75, data = lalonde,
                estimand = "ATT", method = "URI",
                treat = "treat")

lmw.fit1 <- lmw_est(lmw.out1, outcome = "re78")
lmw.fit1

#Influence using SIC
plot(lmw.fit1, type = "influence")

#Usual regression diagnostics
plot(lmw.fit1, type = "lm", which = 1)
}
