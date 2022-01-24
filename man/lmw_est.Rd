\name{lmw_est}
\alias{lmw_est}
\alias{lmw_est.lmw}
\alias{lmw_est.lmw_iv}

\title{
Fit an outcome regression model to estimate a treatment effect
}
\description{
\code{lmw_est()} fits the outcome regression corresponding to the model used to compute the weights in the supplied \code{lmw} object and returns the model coefficients and their covariance matrix. Use \code{\link{summary.lmw_est}} to compute and view the treatment effect and potential outcome mean estimates and their standard errors.
}
\usage{
lmw_est(x, ...)

\method{lmw_est}{lmw}(x, outcome, data = NULL, robust = TRUE, cluster = NULL, ...)

\method{lmw_est}{lmw_iv}(x, outcome, data = NULL, robust = TRUE, cluster = NULL, ...)
}

\arguments{
  \item{x}{
an \code{lmw} or \code{lmw_iv} object; the output of a call to \code{\link{lmw}} or \code{\link{lmw_iv}}.
}
  \item{outcome}{
the name of the outcome variable. Can be supplied as a string containing the name of the outcome variable or as the outcome variable itself. If not supplied, the outcome variable in the \code{formula} supplied to \code{lmw()} or \code{lmw_iv()}, if any, will be used.
}
  \item{data}{
an optional data frame containing the outcome variable named in \code{outcome} and the cluster variable(s) when \code{cluster} is supplied as a \code{formula}.
}
  \item{robust}{
whether to compute the robust covariance matrix for the model coefficients. Allowable values include those allowed for the \code{type} argument of \code{\link[sandwich]{vcovHC}} or \code{\link[sandwich]{vcovCL}} when \code{cluster} is specified. Can also be specified as \code{TRUE} (the default), which means \code{"HC3"} or \code{"HC1"} when \code{cluster} is specified, or \code{FALSE}, which means \code{"const"} (i.e., the standard non-robust covariance). When \code{cluster} is specified, \code{robust} will be set to \code{TRUE} if \code{FALSE}.
}
  \item{cluster}{
the clustering variable(s) for computing a cluster-robust covariance matrix. See \code{\link[sandwich]{vcovCL}}. If supplied as a \code{formula}, the clustering variables must be present in the original dataset used to compute the weights or \code{data}.
}
  \item{\dots}{
other arguments passed to \code{\link[sandwich]{vcovHC}} or \code{\link[sandwich]{vcovCL}}.
}
}
\details{
\code{lmw_est()} uses \code{\link{lm.fit}} or \code{\link{lm.wfit}} to fit the outcome regression model (and first stage model for \code{lmw_iv} objects) and returns the output of these functions augmented with other components related to the estimation of the weights. Unlike with \code{lm.[w]fit()}, the covariance matrix of the parameter estimates is also included in the output.

For \code{lmw} objects, the model fit is that supplied to the \code{formula} input to \code{lmw()} except that it is fit in a dataset appropriately centered to ensure the estimand corresponds with the one requested. When \code{method = "MRI"} in the call to \code{lmw()}, the model is fit as an interaction between the treatment and all the (centered) terms in the model formula. The results will be similar to those from using \code{\link{lm}} on the model formula and supplied data except that the covariates are centered beforehand.

For \code{lmw_iv} objects, the first stage model is constructed by removing the treatment from the supplied model formula, adding the instrumental variable as a main effect, and using the treatment variable as the outcome. For the second stage (reduced form) model, the fitted values of the treatment from the first stage model are used in place of the treatment in the outcome model. The results are similar to those from using \code{ivreg::ivreg()}, and the coefficients estimates will be the same except for the intercept due to the centering of covariates.

Although some coefficients in the model may be interpretable as treatment effect estimates, \code{\link{summary.lmw_est}} should be used to view and extract the treatment effect and potential outcome mean estimates, standard errors, and other model statistics. The output of \code{lmw_est()} should rarely be used except to be supplied to \code{summary()}.
}
\value{
An \code{lmw_est} object with the following components:
\item{coefficients, residuals, fitted.values, effects, weights, rank, df.residual, qr}{for \code{lmw} objects, the output of the \code{\link{lm.fit}} or \code{\link{lm.wfit}} call used to fit the outcome model. For \code{lmw_iv} objects, the output of the \code{\link{lm.fit}} or \code{\link{lm.wfit}} call used to fit the the second stage model, with \code{residuals} corresponding to the residuals computed when using the true treatment variable in place of the fitted treatment values in the model.}
\item{model.matrix}{the model matrix (supplied to the \code{x} argument of \code{lm.fit}).}
\item{vcov}{the estimated covariance matrix of the parameter estimates as produced by \code{\link[sandwich]{vcovHC}} or \code{\link[sandwich]{vcovCL}}.}
\item{call}{the call to \code{lmw_est()}.}
\item{estimand}{the requested estimand.}

Methods exist for \code{lmw_est} objects for \code{\link{model.matrix}}, \code{\link{vcov}}, \code{\link{hatvalues}}, \code{\link[sandwich]{bread}}, and \code{\link[sandwich]{estfun}}, all of which are used internally to compute the parameter estimate covariance matrix. The first two simply extract the corresponding component from the \code{lmw_est} object and the last three imitate the corresponding methods for \code{lm} objects (or \code{ivreg} objects for \code{lmw_iv} inputs). Other regression-related functions, such as \code{\link{coef}}, \code{\link{residuals}}, and \code{\link{fitted}}, use the default methods and should work correctly with \code{lmw_est} objects.
}
\note{
\code{lmw_est()} uses non-standard evaluation to interpret its \code{outcome} argument. For programmers who wish to use \code{lmw_est()} inside other functions, an effective way to pass the name of an arbitrary outcome (e.g., \code{y} passed as a string) is to use \code{\link{do.call}}, for example:
\preformatted{fun <- function(m, y, d) {
   do.call("lmw_est", list(m, y, d))
}
}
When using \code{lmw_est()} inside \code{\link{lapply}} or \code{purrr::map} to loop over outcomes, this syntax must be used as well.
}
\seealso{
\code{\link{summary.lmw_est}} for viewing and extracting the treatment effect and potential outcome mean estimates, standard errors, and other model statistics; \code{\link{lmw}} or \code{\link{lmw_iv}} for estimating the weights that correspond to the model estimated by \code{lmw_est()}; \code{\link{lm}} and \code{\link{lm.fit}} for fitting the corresponding model; \code{ivreg()} in the \pkg{ivreg} package for fitting 2SLS models; \code{\link{influence.lmw_est}} for influence measures
}
\examples{
data("lalonde", package = "MatchIt")

# MRI regression for ATT
lmw.out1 <- lmw(~ treat + age + educ + race + married +
                  nodegree + re74 + re75, data = lalonde,
                  method = "MRI", estimand = "ATT",
                  treat = "treat")

lmw.fit1 <- lmw_est(lmw.out1, outcome = "re78")
lmw.fit1

summary(lmw.fit1)

# MRI regression for ATT after propensity score matching
m.out <- MatchIt::matchit(treat ~ age + educ + race +
                            married + nodegree + re74 + re75,
                          data = lalonde, method = "nearest",
                          estimand = "ATT")
lmw.out2 <- lmw(~ treat + age + educ + race + married +
                  nodegree + re74 + re75, data = lalonde,
                method = "MRI", treat = "treat", obj = m.out)

## Using a cluster-robust SE with subclass (pair membership)
## as the cluster variable
lmw.fit2 <- lmw_est(lmw.out2, outcome = "re78", cluster = ~subclass)
lmw.fit2

summary(lmw.fit2)
}
