\name{lmw}
\alias{lmw}

\title{
Compute regression-implied weights
}
\description{
Computes the weights implied by an outcome regression model that would estimate a weighted difference in outcome means equal to the covariate-adjusted treatment effect resulting from the supplied regression model.
}
\usage{
lmw(formula, data, method = "URI", estimand = "ATE",
    treat = NULL, target = NULL, base.weights = NULL,
    s.weights = NULL, dr.method = "WLS",
    obj = NULL, contrast = NULL, focal = NULL)
}
\arguments{
  \item{formula}{
a one-sided formula with the treatment and covariates on the right-hand side corresponding to the outcome regression model to be fit. If an outcome variable is supplied on the left-hand side, it will be ignored. See Details for how this formula is interpreted in light of other options.
}
  \item{data}{
a data frame containing the variables named in \code{formula} and \code{treat}.
}
  \item{method}{
the method used to estimate the weights; either \code{"URI"} (the default) for uni-regression imputation weights, where a single model is fit to the whole dataset, or \code{"MRI"} for multi-regression imputation, where the model is fit separately in the treatment groups. This affects the interpretation of \code{formula}. See Details.
}
  \item{estimand}{
the estimand of interest, which determines how covariates are centered. Should be one of \code{"ATE"} for the average treatment effect, \code{"ATT"} for the average treatment effect in the treated, \code{"ATC"} for the average treatment effect in the control, or \code{"CATE"} for the conditional average treatment effect. When \code{estimand = "CATE"}, an argument to \code{target} must be supplied. This argument also affects what \code{\link{summary.lmw}} considers to be the target population. Default is \code{"ATE"} unless \code{obj} is specified, in which case it takes its value from the supplied object.
}
  \item{treat}{
the name of the treatment variable in \code{data}. If unspecified, the first variable present in \code{formula} will be taken as the treatment variable with a message. See Details.
}
  \item{target}{
a list or data frame containing the target values for each covariate included in \code{formula}. Ignored with a warning when \code{estimand} is not \code{"CATE"}.
}
  \item{base.weights}{
a vector of base weights. See Details. If omitted and \code{obj} is specified, the weights from the supplied object will be used.
}
  \item{s.weights}{
a vector of sampling weights. See Details. If omitted and \code{obj} is specified, the sampling weights from the supplied object will be used.
}
  \item{dr.method}{
the method used to incorporate the \code{base.weights} into a doubly-robust estimator. Can be one of \code{"WLS"} for weighted least squares or \code{"AIPW"} for augmented inverse probability weighting. Only \code{"WLS"} is allowed with \code{method = "URI"}. Ignored when \code{base.weights} is \code{NULL}. Currently, only \code{"WLS"} is available.
}
  \item{obj}{
a \code{matchit} or \code{weightit} object corresponding to the matched or weighted sample in which the implied outcome regression would take place. See Details.
}
  \item{contrast}{
for multi-category treatments with \code{method = "URI"}, a vector containing the names or indices of the two treatment levels to be contrasted (since in this case the weights depend on the specific contrast). See \code{Details}.
}
  \item{focal}{
the level of the treatment variable to be considered "focal" (i.e., the "treated" level when \code{estimand = "ATT"} or the control level when \code{estimand = "ATC"}). Ignored when \code{estimand} is \code{"ATE"} or \code{"CATE"}. Otherwise, if unspecified, the second value of \code{contrast} will be considered focal when \code{estimand = "ATT"} and the first value of \code{contrast} will be considered focal when \code{estimand = "ATC"}. For binary treatments, this generally does not need to be supplied. See \code{Details}.
}
}
\details{
\code{formula} is interpreted differently depending on whether \code{method} is \code{"URI"} or \code{"MRI"}. When \code{method = "URI"}, the formula is taken literally as the right-hand side of the outcome model formula. The only difference is that the covariates will be centered based on the argument to \code{estimand} (see below). When \code{method = "MRI"}, all references to the treatment are removed (i.e., covariate interactions with treatment become covariate main effects if not already present), and the new formula is taken as the right-hand side of the model formula fit within each treatment group. This is equivalent to allowing all covariates to have both main effects and interactions with treatment after centering the covariates based on the argument to \code{estimand}. Allowing the treatment to interact with all covariates with \code{method = "URI"} is equivalent to specifying \code{method = "MRI"}, and, for binary treatments, the returned weights will be the same.

When any treatment-by-covariate interactions are present in \code{formula} or when \code{method = "MRI"}, covariates are centered at specific values to ensure the resulting weights correspond to the desired estimand as supplied to the \code{estimand} argument. For the ATE, the covariates are centered at their means in the full sample. For the ATT and ATC, the covariates are centered at their means in the treatment or control group (i.e., the \code{focal} group), respectively. For the CATE, the covariates are centered at the values supplied to \code{target}.

When \code{estimand = "CATE"} and \code{target} is specified, the named entries to \code{target} must correspond to all the named \emph{covariates} in \code{formula}; for example, if \code{formula = ~ X1 + log(X1) + X2 + X1:X2}, values in \code{target} must be given for \code{X1} and \code{X2}, but not \code{log(X1)} or \code{X1:X2}. For factor variables, the value in \code{target} must correspond to one of their levels. To choose a target value for a factor corresponding to a proportion (e.g., a target value of .5 for a variable like \code{sex} indicating a target population with a 50-50 sex split), the factor variable must be split into a numeric variable beforehand, e.g., using \code{\link{model.matrix}} or \code{cobalt::splitfactor()}. \code{target} values cannot be given to variables specified using \code{$}, \code{[[]]}, or \code{[]} (e.g., \code{data$X1}), so an error will be thrown if they are used in \code{formula}.

Note that when covariate-by-covariate interactions are present, they will be centered after computing the interaction rather than the interaction being computed on the centered covariates unless \code{estimand = "CATE"}, in which case the covariates will be centered at the values specified in \code{target} prior to involvement in interactions.

\subsection{Base weights and sampling weights}{

Base weights (\code{base.weights}) and sampling weights (\code{s.weights}) are similar in that they both involve combining weights with an outcome regression model. However, they differ in a few ways. Sampling weights are primarily used to adjust the target population; when the outcome model is fit, it is fit using weighted least squares, and when target balance is assessed, it is assessed using the sampling weighted population as the target population. Centering of covariates in the outcome model is done using the sampling weighted covariate means. Base weights are primarily used to offer a second level of balancing beyond the implied regression weights; they can be incorporated into the effect estimator either using weighted least squares or using the augmented inverse probability weighting estimator. Base weights do not change the target population, so when target balance is assessed, it is assessed using the unweighted population as the target population.

Some forms of weights both change the target population and provide an extra layer of balancing, like propensity score weights that target estimands other than the ATT, ATC, or ATE (e.g., overlap weights), or matching weights where the target population is defined by the matching (e.g., matching with a caliper, cardinality matching, or coarsened exact matching). Because these weights change the target population, they should be supplied to \code{s.weights} to ensure covariates are appropriately centered. When there are no treatment-by-covariate interactions and \code{method = "URI"}, whether weights are supplied to \code{base.weights} or \code{s.weights} will not matter for the estimation of the weights but will affect the target population in \link[=summary.lmw]{balance assessment}.

When both \code{base.weights} and \code{s.weights} are supplied, e.g., when the base weights are the result of a propensity score model fit with sampling weights, it is assumed the base weights do not incorporate the sampling weights; that is, it is assumed that to estimate a treatment effect \emph{without} regression adjustment, the base weights and the sampling weights would have to be multiplied together. This is true, for example, for the weights in a \code{matchit} or \code{weightit} object (see below) but not for weights in the output of \code{MatchIt::match.data()} unless called with \code{include.s.weights = FALSE} or weights resulting from \code{CBPS::CBPS()}.
}

\subsection{Regression after using \pkg{MatchIt} or \pkg{WeightIt}}{
Regression weights can be computed in a matched or weighted sample by supplying a \code{matchit} or \code{weightit} object (from \pkg{MatchIt} or \pkg{WeightIt}, respectively) to the \code{obj} argument of \code{lmw()}. The estimand, base weights, and sampling weights (if any) will be taken from the supplied object and used in the calculation of the implied regression weights, unless these have been supplied separately to \code{lmw()}. The \code{weights} component of the supplied object containing the matching or balancing weights will be passed to \code{base.weights} and the \code{s.weights} component will be passed to \code{s.weights}. Arguments supplied to \code{lmw()} will take precedence over the corresponding components in the \code{obj} object.
}

\subsection{Multi-category treatments}{
There are a few differences when the treatment has multiple (i.e., more than 2) categories. If \code{estimand} is \code{"ATT"} or \code{"ATC"}, an argument should be supplied to \code{focal} identifying which group is the treated or control (i.e., "focal") group, respectively.

The key difference, though, is when \code{method = "URI"}, because in this case the contrast between each pair of treatment groups has its own weights and its own implied target population. Because \code{lmw()} only produces one set of weights, an argument must be supplied to \code{contrast} identifying which groups are to be used as the contrast for computing the weights. In addition, to compute the treatment effect corresponding to the chosen contrast as a weighted difference in outcome means, the difference must be taken between the weighted mean of the non-reference group and the weighted mean of \emph{all other groups combined}, rather than simply the weighted mean of the reference group.

The implication of this is that contrast statistics computed in the weighted sample involve all units, even those not in the contrasted groups, whereas statistics computed in the unweighted sample only involve units in the contrasted groups. See \code{\link{summary.lmw}} for more information on assessing balance using the regression weights for multi-category treatments. Given these complications, it is generally best to use \code{method = "MRI"} with multi-category treatments.
}

}
\value{
An \code{lmw} object, which contains the following components:
\item{treat}{the treatment variable, given as a factor.}
\item{weights}{the computed implied regression weights.}
\item{covs}{a data frame containing the covariates included the model formula.}
\item{estimand}{the requested estimand.}
\item{method}{the method used to estimate the weights (\code{"URI"} or \code{"MRI"}).}
\item{base.weights}{the weights supplied to \code{base.weights}.}
\item{s.weights}{the weights supplied to \code{s.weights}.}
\item{dr.method}{when \code{base.weights} are supplied, the method for computing the doubly-robust weights.}
\item{call}{the original call to \code{lmw()}.}
\item{formula}{the model formula.}
\item{target}{the supplied covariate target values when \code{estimand = "CATE"}, after some initial processing.}
\item{contrast}{the contrasted treatment groups.}
\item{focal}{the focal treatment levels when \code{estimand} is \code{"ATT"} or \code{"ATC"}.}
}
\references{
Chattopadhyay, A., & Zubizarreta, J. R. (2021). On the implied weights of linear regression for causal inference. ArXiv:2104.06581 [Stat]. \url{https://arxiv.org/abs/2104.06581}
}

\seealso{
\code{\link{summary.lmw}} for summarizing balance and representativeness; \code{\link{plot.lmw}} for plotting features of the weights; \code{\link{lmw_est}} for estimating treatment effects from \code{lmw} objects; \code{\link{influence.lmw}} for influence measures; \code{\link{lm}} for fitting standard regression models.
}
\examples{
data("lalonde", package = "MatchIt")

#URI regression for ATT
lmw.out1 <- lmw(~ treat + age + educ + race + married +
                  nodegree + re74 + re75, data = lalonde,
                method = "URI", estimand = "ATT",
                treat = "treat")
lmw.out1
summary(lmw.out1)

#MRI regression for ATT
lmw.out2 <- lmw(~ treat + age + educ + race + married +
                  nodegree + re74 + re75, data = lalonde,
                method = "MRI", estimand = "ATT",
                treat = "treat")
lmw.out2
summary(lmw.out2)

#MRI regression for ATT after propensity score matching
m.out <- MatchIt::matchit(treat ~ age + educ + race +
                            married + nodegree + re74 + re75,
                          data = lalonde, method = "nearest",
                          estimand = "ATT")
lmw.out3 <- lmw(~ treat + age + educ + race + married +
                  nodegree + re74 + re75, data = lalonde,
                method = "MRI", treat = "treat", obj = m.out)
lmw.out3
summary(lmw.out3)

#MRI regression for CATE with given target
target <- list(age = 25, educ = 11, race = "black",
               married = 0, nodegree = 1, re74 = 0,
               re75 = 0)
lmw.out4 <- lmw(~ treat + age + educ + race + married +
                  nodegree + re74 + re75, data = lalonde,
                method = "MRI", estimand = "CATE",
                treat = "treat",
                target = target)
lmw.out4
summary(lmw.out4)
}
