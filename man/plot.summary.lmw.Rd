\name{plot.summary.lmw}
\alias{plot.summary.lmw}
\title{
Produce a Love plot of balance statistics
}
\description{
Produces Love plots (also known as dot plots) of balance statistics to summarize balance visually. The plots are generated using \code{\link{dotchart}} and \code{\link{points}}.
}
\usage{
\method{plot}{summary.lmw}(x, abs = TRUE, var.order = "data",
     threshold = NULL, layout = "vertical", ...)
}
\arguments{
  \item{x}{
a \code{summary.lmw} object; the output of a call to \code{\link{summary.lmw}} with \code{standardize = TRUE}.
}
  \item{abs}{
\code{logical}; whether the statistics should be plotted in absolute value (\code{TRUE}) or not (\code{FALSE}). Default is \code{TRUE}.
}
  \item{var.order}{
how the variables should be ordered. Allowable options include \code{"data"}, ordering the variables as they appear in the summary output, and \code{"alphabetical"}, ordering the variables alphabetically. Default is \code{"data"}. Abbreviations allowed.
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
Love plots will be produced for the TSMD Treated and TSMD Control values in the \code{summary.lmw} output. How these plots are arranged depends on the value supplied to \code{layout}, which uses \code{\link{layout}} to arrange the plots.
}
\value{
A plot is displayed, and \code{x} is invisibly returned.
}
\references{
%% ~put references to the literature/web site here ~
}

\seealso{
\code{\link{summary.lmw}}
}
\examples{

}
