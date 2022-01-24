\name{plot.lmw}
\alias{plot.lmw}

\title{
Plots diagnosing the weights
}
\description{
Produces plots to diagnose properties of the weights, including their distribution, to what degree the distribution of covariates involves extrapolation in the weighted sample, and how much influence each unit has on the effect estimate.
}
\usage{
\method{plot}{lmw}(x, type = "weights", ...)
}
\arguments{
  \item{x}{
an \code{lmw} object; the output of a call to \code{\link{lmw}}.
}
  \item{type}{
the type of plot to display. Allowable options include \code{"weights"}, \code{"extrapolation"}, and \code{"influence"}. See Details. Abbreviations allowed.
}
  \item{\dots}{
further arguments passed to specific types of plots.

When \code{type = "weights"}, the following are accepted:
\describe{
  \item{\code{rug}}{\code{logical}; whether to display a rug plot of the weights. Default is \code{TRUE}.}
  \item{\code{mean}}{whether to display a red line indicating the mean of the weights. Default is \code{TRUE}.}
  \item{\code{ess}}{whether to display the original and weighted effective sample size in the top right corner. Default is \code{TRUE}.}
  }
Other arguments are passed to \code{\link{density}}.

When \code{type = "extrapolation"}, the following are accepted:
\describe{
  \item{\code{var}}{required; a right-sided formula or character vector containing the names of the covariates for which extrapolation is to be assessed.}
  \item{\code{data}}{an optional data frame containing the variables named in \code{var}.}
}
When \code{type = "influence"}, the following are accepted:
\describe{
  \item{\code{outcome}}{the name of the outcome variable. Can be supplied as a string containing the name of the outcome variable or as the outcome variable itself. If not supplied, the outcome variable in the \code{formula} supplied to \code{lmw()}, if any, will be used.}
  \item{\code{data}}{an optional data frame containing the outcome variable named in \code{outcome}.}
}
}
}
\details{
When \code{type = "weights"}, \code{plot.lmw()} produces a density plot of the weights within each treatment group. By construction, these weights will have a mean of 1. Some weights may be negative. The effective sample size (ESS) and original sample size (N) will be displayed in the upper right corner of the plot when \code{ess = TRUE}.

When \code{type = "extrapolation"}, \code{plot.lmw()} produces a plot of the distribution of weights and covariates for each treatment group. Each dot represents a unit, with values arranged on the x-axis according to their covariate value and the size of the dots corresponding to the magnitude of the weight. Units with positive weights are displayed in black in the upper portion of the plot, and units with negative weights are displayed in red in the lower portion. Having many and large red points indicates a high degree of extrapolation. All points are equally transparent, so darker regions indicate multiple points with the same value. The vertical lines indicates the weighted mean of the covariate in each group, and the X indicates the mean of the covariate in the target sample as determined by the \code{estimand} argument in the original call to \code{lmw()}. A large discrepancy between the vertical lines and Xs indicates a lack of balance between the treatment group and target sample.

When \code{type = "influence"}, \code{plot.lmw()} produces a plot of the scaled sample influence curve (SIC) for each unit by index. It does so by calling \code{\link{influence.lmw}}, which fits the outcome model to extract residuals and compute the SIC as \code{SIC = (N-1) * w * r / (1 - h)}, where \code{N} is the sample size, \code{w} are the units' implied regression weights, \code{r} are the residuals, and \code{h} are the hat values. SIC values are scaled to have a maximum of 1. Higher values indicate greater relative influence.
}
\value{
A plot is displayed, and \code{x} is invisibly returned.
}

\seealso{
\code{\link{lmw}}, \code{\link{summary.lmw}}, \code{\link{plot.summary.lmw}}
}
\examples{
data("lalonde")

#URI regression for ATT
lmw.out1 <- lmw(~ treat + age + educ + race + married +
                   nodegree + re74 + re75, data = lalonde,
                estimand = "ATT", method = "URI",
                treat = "treat")
lmw.out1

#Distribution of weights
plot(lmw.out1, type = "weights")

#Extrapolation/representativeness for age and married
plot(lmw.out1, type = "extrapolation",
     var = ~age + married)

#Extrapolation/representativeness for race
plot(lmw.out1, type = "extrapolation",
     var = ~race)

#Influence for re78 outcome
plot(lmw.out1, type = "influence", outcome = "re78")
}
