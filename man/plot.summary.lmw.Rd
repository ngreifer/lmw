\name{plot.summary.lmw}
\alias{plot.summary.lmw}
\title{
Produce a Love plot of balance statistics
}
\description{
Produces Love plots (also known as dot plots) of balance statistics to summarize balance visually. The plots are generated using \code{\link{dotchart}} and \code{\link{points}}.
}
\usage{
\method{plot}{summary.lmw}(x,
     stats, abs = TRUE, var.order = "data",
     threshold = NULL, layout = "vertical", ...)
}
\arguments{
  \item{x}{
a \code{summary.lmw} object; the output of a call to \code{\link{summary.lmw}} with \code{standardize = TRUE}.
}
  \item{stats}{
a vector of the names of the columns in the \code{summary.lmw} output to plot; more than one is allowed. Abbreviations allowed.
}
  \item{abs}{
\code{logical}; whether the statistics should be plotted in absolute value (\code{TRUE}) or not (\code{FALSE}). Default is \code{TRUE}. This does not affect the display of KS statistics (which are always non-negative).
}
  \item{var.order}{
how the variables should be ordered. Allowable options include \code{"data"}, ordering the variables as they appear in the summary output, \code{"alphabetical"}, ordering the variables alphabetically, and, when \code{un = TRUE} in the call to \code{summar.lmw()}, \code{"unadjusted"}, ordering the variables by the first statistic in \code{stats} in the unadjusted sample. Default is \code{"data"}. Abbreviations allowed.
}
  \item{threshold}{
numeric values at which to place vertical lines indicating a balance threshold. These can make it easier to see for which variables balance has been achieved given a threshold. Multiple values can be supplied to add multiple lines. When \code{abs = FALSE}, the lines will be displayed on both sides of zero. The lines are drawn with \code{abline} with the linetype (\code{lty}) argument corresponding to the order of the entered variables (see options at \code{\link{par}}). Enter a value as \code{NA} to skip that value of \code{lty} (e.g., \code{c(NA, .05)} to have only a dashed vertical line at .05).
}
  \item{layout}{
how the multiple plots should be laid out. Allowable options include \code{"vertical"} (the default) and \code{"horizontal"}. Abbreviations allowed.
}
  \item{\dots}{
further arguments passed to \code{\link{dotplot}}.
}
}
\details{
Love plots will be produced for the requested statistics in the \code{summary.lmw} output. How these plots are arranged depends on the value supplied to \code{layout}, which uses \code{\link{layout}} to arrange the plots.
}
\value{
A plot is displayed, and \code{x} is invisibly returned.
}

\seealso{
\code{\link{summary.lmw}}
}
\examples{
data("lalonde")

#URI regression for ATT
lmw.out1 <- lmw(~ treat + age + educ + race + married +
                  nodegree + re74 + re75, data = lalonde,
                method = "URI", estimand = "ATT", treat = "treat")
lmw.out1
(s <- summary(lmw.out1))

plot(s)
plot(s, stats = "SMD", abs = FALSE)
}
