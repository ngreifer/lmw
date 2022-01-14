\name{summary.lmw_est}
\alias{summary.lmw_est}

\title{
Extract effect estimates and standard errors from \code{lmw_est} fits
}
\description{
\code{summary()} computes the treatment effect and potential outcome mean estimates from the supplied \code{lmw_est} object. It functions similarly to \code{\link{summary.lm}} in producing estimate tables with the estimates, standard errors, t-statistics, and p-values. Other model statistics can be additionally requested.
}
\usage{
\method{summary}{lmw_est}(object, model = FALSE, ci = TRUE,
        alpha = .05, ...)
}
\arguments{
  \item{object}{
an \code{lmw_est} object; the output of a call to \code{lmw_est}.
}
  \item{model}{
\code{logical}; whether to produce a coefficient table for the model coefficients other than the treatment effect estimate. Note that these values should not be interpreted or reported so they are not produced by default.
}
  \item{ci}{
\code{logical}; whether to include confidence intervals in the output.
}
  \item{alpha}{
when \code{ci = TRUE}, the alpha value used to compute the critical test statistic for the confidence interval; equivalently, 1 minus the confidence level (e.g., for a 99\% confidence interval, \code{alpha = .01} should be specified). Default is .05 for a 95\% confidence interval.
}
  \item{\dots}{
ignored.
}
}
\details{
\code{summary.lmw_est()} produces a table of treatment effect estimates corresponding to all possible pairwise contrasts between the treatment levels. These treatment effects generalize to the population implied by the regression weights, which depends on the supplied estimand, whether sampling weights were provided, and which of the MRI or URI models was requested. The treatment effects are computed using linear contrasts of the outcome model coefficients.

When \code{type = "MRI"}, the potential outcome mean estimates are also reported. These correspond to the potential outcome means in the population implied by the regression weights. When \code{type = "URI"}, only the treatment effects are estimated; the model-implied outcome means do not correspond to the potential outcome means for the population implied by the regression weights. That is, while the treatment effect generalizes to the population defined by the regression weights, the estimated potential outcome means do not and so are not reported.

When \code{model = TRUE}, the model coefficients and their tests statistics are additionally produced. It is inappropriate to interpret or report these values as they have no causal interpretation.
}
\value{
A \code{summary.lmw_est} object with the following components:
\item{call}{The original call to \code{lmw_est()}}
\item{coefficients}{A matrix the treatment effect and potential outcome mean estimates, and their standard errors, t-statistics, and p-values. The row names of the matrix are \code{"E[Y0]"}, \code{"E[Y1]"}, and \code{"E[Y1-Y0]"}, corresponding to the potential outcome mean under control, the potential outcome mean under treatment, and the treatment effect, regardless of how they are printed. The columns are \code{"Estimate"}, \code{"Std. Error"}, \code{"t value"}, and \code{"Pr(>|t|)"}. When \code{ci = TRUE}, the confidence limits \code{"95\%" CI L} (lower) and \code{"95\%" CI U} (upper) will be included between \code{"Std. Error"} and \code{"t value"}.}
\item{model.coefficients}{When \code{model = TRUE}, the coefficient table of the model coefficients, which has the same columns as \code{coefficients.}}
\item{aliased}{When \code{model = TRUE}, a named logical vector showing if the original coefficients are aliased (i.e., \code{NA}).}
\item{sigma, df, r.squared, adj.r.squared}{The residual standard deviation, degrees of freedom components, R-squared, and adjusted R-squared. See \code{\link{summary.lm}}.}
\item{estimand}{The requested estimand}
When \code{type = "URI"}, the mean potential outcomes (\code{"E[Y0]"} and \code{"E[Y1]"}) and their corresponding estimates will be set to \code{NA} and they will not be printed in the output.
}

\seealso{
\code{\link{lmw_est}} for fitting the outcome regression model, \code{\link{summary.lm}} for more information on the output components
}
\examples{
# See examples at help("lmw_est")
}
