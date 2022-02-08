\name{lmw_iv}
\alias{lmw_iv}

\title{
Compute instrumental variable regression-implied weights
}
\description{
Computes the weights implied by an instrumental variable model that would estimate a weighted difference in outcome means equal to the treatment effect resulting from the supplied model fit with two-stage least squares.
}
\usage{
lmw_iv(formula, data, estimand = "ATE", method = "URI",
       treat = NULL, iv, base.weights = NULL,
       s.weights = NULL, obj = NULL, target = NULL,
       target.weights = NULL, contrast = NULL,
       focal = NULL)
}
\arguments{
  \item{formula}{
a one-sided formula with the treatment and covariates on the right-hand side corresponding to the second-stage (reduced form) outcome regression model to be fit. If an outcome variable is supplied on the left-hand side, it will be ignored. This model should not include an instrumental variable. See Details for how this formula is interpreted in light of other options.
}
  \item{data}{
a data frame containing the variables named in \code{formula}, \code{treat}, and \code{iv}.
}
  \item{estimand}{
the estimand of interest, which determines how covariates are centered. Should be one of \code{"ATE"} for the average treatment effect, \code{"ATT"} for the average treatment effect in the treated, \code{"ATC"} for the average treatment effect in the control, or \code{"CATE"} for the conditional average treatment effect. When \code{estimand = "CATE"}, an argument to \code{target} must be supplied. This argument only affects what \code{\link{summary.lmw}} considers to be the target population. Default is \code{"ATE"} unless \code{obj} is specified, in which case it takes its value from the supplied object.
}
  \item{method}{
the method used to estimate the weights; currently, only \code{"URI"} for uni-regression imputation weights, where a single model is fit to the whole dataset, is allowed.
}
  \item{treat}{
the name of the treatment variable in \code{data}. If unspecified, the first variable present in \code{formula} will be taken as the treatment variable with a message. Currently, only binary treatments are supported. See Details.
}
  \item{iv}{
a character vector or one-sided formula containing the names of the instrumental variables in \code{data}. These variables should not appear in \code{formula}. See Details.
}
  \item{base.weights}{
a vector of base weights. See Details. If omitted and \code{obj} is specified, the weights from the supplied object will be used.
}
  \item{s.weights}{
a vector of sampling weights. See Details. If omitted and \code{obj} is specified, the sampling weights from the supplied object will be used.
}
  \item{obj}{
a \code{matchit} or \code{weightit} object corresponding to the matched or weighted sample in which the implied outcome regression would take place. See Details.
}
  \item{target}{
a list or data frame containing the target values for each covariate included in \code{formula}. Ignored with a warning when \code{estimand} is not \code{"CATE"}.
}
  \item{target.weights}{
a vector of sampling weights to be applied to \code{target} when supplied as a data frame. Ignored with a warning when \code{estimand} is not \code{"CATE"}.
}
  \item{contrast}{
ignored.
}
  \item{focal}{
the level of the treatment variable to be considered "focal" (i.e., the "treated" level when \code{estimand = "ATT"} or the control level when \code{estimand = "ATC"}). Ignored when \code{estimand} is \code{"ATE"} or \code{"CATE"}.
}
}
\details{
\code{lmw_iv()} computes weights that, when used in the weighted difference in outcome means between the treatment groups, corresponds to the two-stage least squares (2SLS) estimate of the treatment effect. \code{formula} corresponds to the second-stage (reduced form) model, with the treatment replaced by its fitted values resulting from the first stage model. The first stage is fit by removing the treatment from the supplied \code{formula}, adding the instrumental variables named in \code{iv} as main effects, and using the treatment as the outcome. The treatment is assumed to be endogenous and the supplied instrumental variables assumed to be instruments conditional on the other covariates, which are assumed to to be exogenous.

Centering covariates (i.e., by supplying arguments to \code{estimand}, \code{target}, or \code{focal}) has no effect on the resulting weights or the treatment effect estimated. It will affect the intercept of the outcome model and the target population to which the weighted distribution of covariates is compared using \code{\link{summary.lmw}} and \code{\link{plot.lmw}}. See \code{\link{lmw}} for more information on covariate centering.

\subsection{Base weights and sampling weights}{

Base weights (\code{base.weights}) and sampling weights (\code{s.weights}) are similar in that they both involve combining weights with an outcome regression model. However, they differ in a few ways. Sampling weights are primarily used to adjust the target population; when the outcome model is fit, it is fit using weighted least squares, and when target balance is assessed, it is assessed using the sampling weighted population as the target population. Centering of covariates in the outcome model is done using the sampling weighted covariate means. Base weights are primarily used to offer a second level of balancing beyond the implied regression weights, i.e., to fit the 2SLS models in the base-weighted sample. Base weights do not change the target population, so when target balance is assessed, it is assessed using the unweighted population as the target population.

Some forms of weights both change the target population and provide an extra layer of balancing, like propensity score weights that target estimands other than the ATT, ATC, or ATE (e.g., overlap weights), or matching weights where the target population is defined by the matching (e.g., matching with a caliper, cardinality matching, or coarsened exact matching). Because these weights change the target population, they should be supplied to \code{s.weights} to ensure covariates are appropriately centered. When there are no treatment-by-covariate interactions and \code{method = "URI"}, whether weights are supplied to \code{base.weights} or \code{s.weights} will not matter for the estimation of the weights but will affect the target population in \link[=summary.lmw]{balance assessment}.

When both \code{base.weights} and \code{s.weights} are supplied, e.g., when the base weights are the result of a propensity score model fit with sampling weights, it is assumed the base weights do not incorporate the sampling weights; that is, it is assumed that to estimate a treatment effect \emph{without} regression adjustment, the base weights and the sampling weights would have to be multiplied together. This is true, for example, for the weights in a \code{matchit} or \code{weightit} object (see below) but not for weights in the output of \code{MatchIt::match.data()} unless called with \code{include.s.weights = FALSE} or weights resulting from \code{CBPS::CBPS()}.
}

\subsection{2SLS after using \pkg{MatchIt} or \pkg{WeightIt}}{
Instrumental variable regression weights can be computed in a matched or weighted sample by supplying a \code{matchit} or \code{weightit} object (from \pkg{MatchIt} or \pkg{WeightIt}, respectively) to the \code{obj} argument of \code{lmw()}. The estimand, base weights, and sampling weights (if any) will be taken from the supplied object and used in the calculation of the implied regression weights, unless these have been supplied separately to \code{lmw_iv()}. The \code{weights} component of the supplied object containing the matching or balancing weights will be passed to \code{base.weights} and the \code{s.weights} component will be passed to \code{s.weights}. Arguments supplied to \code{lmw_iv()} will take precedence over the corresponding components in the \code{obj} object.
}

\subsection{Multi-category treatments and instruments}{

Multi-category treatments are not currently supported for \code{lmw_iv()}.
}

}

\value{
An \code{lmw_iv} object, which inherits from \code{lmw} objects and contains the following components:
\item{treat}{the treatment variable, given as a factor.}
\item{weights}{the computed implied regression weights.}
\item{covs}{a data frame containing the covariates included the model formula.}
\item{estimand}{the requested estimand.}
\item{method}{the method used to estimate the weights (\code{"URI"} or \code{"MRI"}).}
\item{base.weights}{the weights supplied to \code{base.weights}.}
\item{s.weights}{the weights supplied to \code{s.weights}.}
\item{call}{the original call to \code{lmw_iv()}.}
\item{formula}{the model formula.}
\item{iv}{the instrumental variables, given as a one-sided formula.}
\item{target}{the supplied covariate target values when \code{estimand = "CATE"}, after some initial processing.}
\item{contrast}{the contrasted treatment groups.}
\item{focal}{the focal treatment levels when \code{estimand} is \code{"ATT"} or \code{"ATC"}.}

All functions that lack a specific \code{lmw_iv} method work with \code{lmw_iv} objects as they do for \code{lmw} objects, link \code{\link{summary.lmw}}, \code{\link{plot.lmw}}, etc.
}
\references{
Chattopadhyay, A., & Zubizarreta, J. R. (2021). On the implied weights of linear regression for causal inference. ArXiv:2104.06581 [Stat]. \url{https://arxiv.org/abs/2104.06581}
}

\seealso{
\code{\link{summary.lmw}} for summarizing balance and representativeness; \code{\link{plot.lmw}} for plotting features of the weights; \code{\link{lmw_est}} for estimating treatment effects from \code{lmw_iv} objects; \code{\link{influence.lmw}} for influence measures; \code{ivreg()} in the \pkg{ivreg} package for fitting 2SLS models.
}
\examples{
#stop("Examples needed for lmw_iv")
}
