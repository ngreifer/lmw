#' Data from National Supported Work Demonstration and PSID, as analyzed by
#' Dehejia and Wahba (1999).
#'
#' This is a subsample of the data from the treated group in the National
#' Supported Work Demonstration (NSW) and the comparison sample from the
#' Population Survey of Income Dynamics (PSID). This data was previously
#' analyzed extensively by Lalonde (1986) and Dehejia and Wahba (1999). The
#' original dataset is available at
#' \url{https://users.nber.org/~rdehejia/nswdata2.html}.
#'
#' The data corresponds to the NSW treated sample and the PSID control sample
#' with 1974 earnings included. This specific dataset is different from the one in the
#' \pkg{MatchIt} and \pkg{cobalt} packages, which is a subset of this dataset.
#'
#'
#' @docType data
#' @format A data frame with 2675 observations (185 treated, 2490 control). There
#' are 9 variables measured for each individual. In addition, two constructed variables are included: \code{treat_multi}, which splits the original control group into two, and \code{Ins}, which is a constructed instrumental variable.
#' \itemize{
#' \item "treat" is the
#' treatment assignment (1=treated, 0=control).
#'  \item "age" is age in years.
#'  \item "education" is education in number of years of schooling.
#'  \item "race" is
#' the individual's race/ethnicity, (Black, Hispanic, or White).
#' \item "married" is an indicator for married (1=married, 0=not married).
#' \item "nodegree" is an indicator for whether the individual lacks a high school
#' degree (i.e., has fewer than 12 years of schooling; 1=no degree, 0=degree).
#' \item "re74" is income in 1974, in U.S. dollars.
#' \item "re75" is income in
#' 1975, in U.S. dollars.
#' \item "re78" is income in 1978, in U.S. dollars.
#' \item "treat_multi" is a constructed version of "treat" that splits the control group into to levels (1=treated, 2=control group A, 3=control group B).
#' \item "Ins" is a binary instrumental variable.
#'
#'  }
#' "treat" is the treatment variable, "re78" is the outcome, and the others are
#' pre-treatment covariates. Note that in the original data, "race" is instead
#' coded as two dummy variables, "black" and "hispan".
#' @references Lalonde, R. (1986). Evaluating the econometric evaluations of
#' training programs with experimental data. American Economic Review 76:
#' 604-620.
#'
#' Dehejia, R.H. and Wahba, S. (1999).  Causal Effects in Nonexperimental
#' Studies: Re-Evaluating the Evaluation of Training Programs.  Journal of the
#' American Statistical Association 94: 1053-1062.
#' @keywords datasets
"lalonde"
