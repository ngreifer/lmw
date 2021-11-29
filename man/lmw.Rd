\name{lmw}
\alias{lmw}

\title{
Compute regression-implied weights
}
\description{
Computes the weights implied by an outcome regression model that would estimate a weighted difference in outcome means equal the covariate-adjusted treatment effect resulting from the supplied regression model.
}
\usage{
lmw(formula, data, type = "URI", estimand = "ATE",
    treat = NULL, target = NULL, base.weights = NULL,
    s.weights = NULL)
}
\arguments{
  \item{formula}{
a one-sided formula with the treatment and covariates on the right-hand side corresponding to the outcome regression model to be fit. If an outcome variable is supplied on the left-hand side, it will be ignored. See Details for how this formula is interpreted in light of other options.
}
  \item{data}{
a data frame containing the variables named in \code{formula} and \code{treat}.
}
  \item{type}{
the type of weights estimated; either \code{"URI"} (the default) for uni-regression imputation weights, where a single model is fit to the whole dataset, or \code{"MRI"} for multi-regression imputation, where the model is fit seperately in the treatment groups. This affects the interpretation of \code{formula}. See Details.
}
  \item{estimand}{
the estimand of interest, which determines how covariates are centered. Should be one of \code{"ATE"} for the average treatment effect, \code{"ATT"} for the average treatment effect in the treated, \code{"ATC"} for the average treatment effect in the control, or \code{"CATE"} for the conditional average treatment effect. When \code{estimand = "CATE"}, an argument to \code{target} must be supplied. This argument also affects what \code{\link{summary.lmw}} considers to be the target population. Default is \code{"ATE"}.
}
  \item{treat}{
the name of the treatment variable in \code{data}. If unspecified, the first variable present in \code{formula} will be taken as the treatment variable with a message.
}
  \item{target}{
%%     ~~Describe \code{target} here~~
}
  \item{base.weights}{
a vector of base weights. See Details.
}
  \item{s.weights}{
a vector of sampling weights. See Details.
}
}
\details{
\code{formula} is interpreted differently depending on whether \code{type} is \code{"URI"} or \code{"MRI"}, When \code{type = "URI"}, the formula is taken literally as the right-hand side of the outcome model formula. The only difference is that the covariates will be centered based on the argument to \code{estimand} (see below). When \code{type = "MRI"}, all references to the treatment are removed (i.e., covariate interactions with treatment become covariate main effects if not already present), and the new formula is taken as the right-hand side of the model formula fit within each treatment group. This is equivalent to allowing all covariates to have both main effects and interactions with treatment after centering the covariates based on the argument to \code{estimand}. Allowing the treatment to interact with all covariates with \code{type = "URI"} is equivalent to specifying \code{type = "MRI"}, and the returned weights will be the same.

When any treatment-by-covariate inteactions are present in \code{formula} or when \code{type = "MRI"}, covariates are centered at specific values to ensure the resulting weights correspond to ther desired estimand as supplied to the \code{estimand} argument. For the ATE, the covariates are centered at their means in the full sample. For the ATT and ATC, the covariates are centered at their means in the treatment or control group, respectively. For the CATE, the covariates are centered at the values supplied to \code{target}.

Note that when covariate-by-covariate interactions are present, they will be centered after computing the interaction rather than the interaction being computed on the centered covariates unless \code{estimand = "CATE"}, in which case the covariates will be centered prior to involvement in interactions.

\subsection{Base weights and sampling weights}{

Base weights (\code{base.weights}) and sampling weights (\code{s.weights}) are similar in that they both imply an outcome model fit with weighted least squares. However, they differ in the following way: sampling weights change the target population of the estimand to the weighted population, whereas base weights retain the original target population. If the weights are the result of a first-stage estimation procedure that aims to target the same estimand as the implied regression weights, they should be supplied to \code{base.weights}. If the weights are sampling/survey weights or otherwise are meant to change the target population, they should be supplied to \code{s.weights}.

Some examples of types of weights that should be supplied to \code{base.weights} are propensity score weights for the ATE, ATT, or ATC and matching weights when a specific estimand is targeted (i.e., when no calipers are used). Some examples of types of weights that should be supplied to \code{s.weights} are survey/sampling weights, propensity score weights for estimands defined by the weights (e.g., overlap weights), and matching weights where the target population is defined by the matching (e.g., matching with a caliper, cardinality matching, or coarsened exact matching). When in doubt, it is usually best to use \code{s.weights}. When there are no treatment-by-covariate interactions and \code{type = "URI"}, whether weights are supplied to \code{base.weights} or \code{s.weights} will not matter for the estimation of the weights but will affect the target population in balance assessment and the estimation of the potential outcome means.
}
}
\value{
An \code{lmw} object, which contains the following components:
\item{treat}{the treatment variable.}
\item{weights}{the computed weights.}
\item{covs}{a data frame containing the covariates included the model formula.}
\item{estimand}{the request estimand.}
\item{type}{the type of weights (\code{"URI"} or \code{"MRI"}).}
\item{s.weights}{the weights supplied to \code{s.weights}.}
\item{formula}{the model formula.}
\item{call}{the original call to \code{lmw()}.}
}
\references{
%% ~put references to the literature/web site here ~
}

\seealso{
\code{\link{lm}}
}
\examples{

}
