\name{summary.lmw}
\alias{summary.lmw}

\title{
Assess balance on an lmw object
}
\description{
Computes balance statsitics for an \code{lmw} object created by \code{\link{lmw}}. Balance involves not only the similarity between the treatment groups but also the similarity between each treatment group and the target population.
}
\usage{
\method{summary}{lmw}(object, un = TRUE, addlvariables = NULL,
        standardize = TRUE, data = NULL, ...)
}

\arguments{
  \item{object}{
an \code{lmw} object; the output of a call to \code{\link{lmw}}.
}
  \item{un}{
\code{logical}; whether to display balance statistics for the sample prior to weighting. If \code{s.weights} were supplied to \code{lmw}, the unadjusted sample will be weighted by the sampling weights.
}
  \item{addlvariables}{
additional variable for which balance statistics are to be computed along with the covariates in the \code{lmw} object. Can be entered in one of three ways: as a data frame of covariates with as many rows as there were units in the original \code{lmw} call, as a string containing the names of variables in \code{data}, or as a right-sided formula with the additional variables (and possibly their transformations) found in \code{data}, the environment, or the \code{lmw} object.
}
  \item{standardize}{
\code{logical}; whether to compute standardized (\code{TRUE}) or unstandardized (\code{FALSE}) mean differences. Default is \code{TRUE}.
}
  \item{data}{
a optional data frame containing variables named in \code{addlvariables} if specified as a string or formula.
}
  \item{\dots}{
ignored.
}
}
\details{
\code{summary.lmw} produces covariate balance statistics and effective samples size before and after weighting. For each covariate, the following statistics are computed: the standardized mean difference between the treated and control groups (SMD), the standardized mean difference between the treated group and target sample (TSMD Treated), the standardized mean difference between between the control group and target sample (TSMD Control), the Kolmogorov-Smirnov statistic between the treated and control groups (KS), the the Kolmogorov-Smirnov statistic between the treated group and target sample (KS Treated), and the the Kolmogorov-Smirnov statistic between the control group and target sample (KS Control). For binary covariates, the Kolmogorov-Smirnov statistic is equal to the unstandardized difference in means and is computed as such.

After weighting, the standardized mean difference between the treated and control groups of each covariate included in the original call to \code{lmw()} will be equal to zero. However, the difference between each treatment group and the target sample may not be equal to zero when \code{type = "URI"} in the call to \code{lmw()}, and covariates supplied to \code{addlvariables} not included in the call to \code{lmw()} may not be well balanced.

When \code{s.weights} are supplied to \code{lmw()}, the unadjusted balance statistics (if requested) will incorporate the sampling weights. When \code{base.weights} are supplied to \code{lmw()}, the unadjusted balance statistics will \emph{not} incorporate the sampling weights.

Standardized mean differences are computed as the difference between the (weighted) means divided by a standardization factor, which is the standard deviation of the covariate in the target sample. When \code{estimand = "ATT"} in the call to \code{lmw()}, the standardization factor is the standard deviation in the treated group; when \code{estimand = "ATC"}, the standardization factor is the standard deviation in the control group; when \code{estimand = "ATE"} or \code{"CATE"}, the standardization factor is the square root of the average of the variances of both treatment groups. When \code{s.weights} is supplied, the standardization factor is computed including the sampling weights; otherwise it is computed in the unweighted sample.

The effective sample size (ESS) is computed as \eqn{(\sum|w|)^2/\sum w^2}. With uniform weights, this is equal to the sample size.
}
\value{
A \code{summary.lmw} object, which contains the following components:
\item{call}{The original call to \code{lmw()}.}
\item{nn}{The (effective) sample sizes before and after weighting.}
\item{sum.un}{If \code{un = TRUE}, the balance statistics prior to weighting.}
\item{sum.weighted}{The balance statistics after weighting.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}

\seealso{
\code{\link{lmw}}
}
\examples{

}
