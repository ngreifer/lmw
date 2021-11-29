\name{plot.lmw}
\alias{plot.lmw}

\title{
Plots diagnosing the weights
}
\description{
Produces plots to diagnose properties of the weights, including their distribution and to what degree the distribution of covariates involves extrapolation in the weighted sample.
}
\usage{
\method{plot}{lmw}(x, type = "weights", ...)
}
\arguments{
  \item{x}{
an \code{lmw} object; the output of a call to \code{\link{lmw}}.
}
  \item{type}{
the type of plot to display. Allowable options include \code{"weights"} and \code{"extrapolation"}. See Details. Abbreviations allowed.
}
  \item{\dots}{
further arguments passed to specific types of plots. When \code{type = "weights"}, the following are accepted:
\describe{
  \item{\code{rug}}{\code{logical}; whether to display a rug plot of the weights. Default is \code{TRUE}.}
  \item{\code{mean}}{whether to display a red line indicating the mean of the weights. Default is \code{TRUE}.}
  }
When \code{type = "extrapolation"}, the follwing is required:
\describe{
  \item{\code{var}}{the name of the covariate in the original call to \code{lmw()} for which extrapolation is to be assessed.}
}
}
}
\details{
When \code{type = "weights"}, \code{plot.lmw()} produces a density plot of the weights within each treatment group. By construction, these weights will have a mean of 1. Some weights may be negative. The effective sample size (ESS) and original sample size (N) will be displayed in the upper left corner of the plot.

When \code{type = "extrapolation"}, \code{plot.lmw()} produces a plot of the distribution of weights and covariates for each treatment group. Each dot represents a unit, with values arrange on the x-axis according to their covariate value and the size of the dots corresponding to the magnitude of the weight. Units with positive weights are displayed in black in the upper portion of the plot, and units with negative weights are displayed in red in the lower portion. Having many and large red points indicates a high degree of extrapolation. All points are equally transparent, so darker regions indicate multiple points with the same value. The vertical lines indicates the weighted mean of the covariate in the treatment group, and the X indicates the mean of the covariate in the target sample as determined by the \code{estimand} argument in the original call to \code{lmw()}. A large discrepency between the vertical lines and Xs indicates a lack of balance between teh treatment group and target sample.
}
\value{
A plot is displayed, and \code{x} is invisibly returned.
}
\references{
%% ~put references to the literature/web site here ~
}

\seealso{
\code{\link{lmw}}, \code{\link{summary.lmw}}, \code{\link{plot.summary.lmw}}
}
\examples{

}
