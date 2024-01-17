
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lmw: Linear Model Weights

<!-- <img src="man/figures/logo.png" align="right" width="150"/> -->

## [![CRAN_Status_Badge](https://img.shields.io/cran/v/lmw)](https://cran.r-project.org/package=lmw) [![CRAN_Downloads_Badge](https://cranlogs.r-pkg.org/badges/lmw)](https://cran.r-project.org/package=lmw)

### Overview

`lmw` computes weights implied by a linear regression model used to
estimate an average treatment effect and provides diagnostics that
incorporate these weights as described in [Chattopadhyay & Zubizarreta
(2023)](https://doi.org/10.1093/biomet/asac058). The treatment effect
resulting from this model can be represented as a difference between the
weighted outcome means in the treatment groups, similar to inverse
probability weighting.

The weights implied by the regression model have several interesting
characteristics: they yield exact mean balance between the treatment
groups for every covariate included in the model, they have minimum
variance among all weights that do so, and they may be negative. Despite
yielding exact mean balance between treatment groups, they may not yield
balance between the treatment groups and the target population
corresponding to the desired estimand; in addition, the negative weights
they produce indicate extrapolation beyond the support of the
covariates. `lmw` provides tools to compute the implied weights and
perform diagnostics to assess balance, extrapolation, influence, and
distributional properties of the weights. In addition, `lmw` provides
tools to estimate average treatment effects from the specified models
that correspond to selected target populations.

Below is an example of the use of `lmw` to obtain and evaluate
regression weights for estimating the average treatment effect on the
treated (ATT) of a job training program on earnings using the Lalonde
dataset (see `help("lalonde", package = "lmw")` for details). Here, the
treatment variable is `treat`, the outcome is `re78` (1978 earnings) and
`age`, `education`, `race`, `married`, `nodegree`, `re74`, and `re75`
are pretreatment covariates.

``` r
#Load lmw and the data
library("lmw")
data("lalonde")

#Estimate the weights
lmw.out <- lmw(~ treat + age + education + race + married +
                 nodegree + re74 + re75,
               data = lalonde, treat = "treat",
               estimand = "ATT", method = "URI")
```

`lmw.out` contains the implied regression weights. We can see that the
weighted outcome mean difference between the treatment groups is equal
to the coefficient on the treatment in the corresponding outcome
regression model:

``` r
#The weighted difference in means
w <- lmw.out$weights
with(lalonde, 
     weighted.mean(re78[treat == 1], w[treat == 1]) - 
       weighted.mean(re78[treat == 0], w[treat == 0]))
```

    #> [1] 751.9464

``` r
#The coefficient in the corresponding outcome model
fit <- lm(re78 ~ treat + age + education + race + married +
            nodegree + re74 + re75,
          data = lalonde)
coef(fit)["treat"]
```

    #>    treat 
    #> 751.9464

How well do the weights do at balancing the treatment groups to the
target population (in the case, the treated sample)?

``` r
(s <- summary(lmw.out))
```

    #> 
    #> Call:
    #> lmw(formula = ~treat + age + education + race + married + nodegree + 
    #>     re74 + re75, data = lalonde, estimand = "ATT", method = "URI", 
    #>     treat = "treat")
    #> 
    #> Summary of Balance for Unweighted Data:
    #>                 SMD TSMD Control TSMD Treated    KS TKS Control TKS Treated
    #> age          -1.263        1.263            0 0.377       0.377           0
    #> education    -0.881        0.881            0 0.403       0.403           0
    #> raceblack     1.630       -1.630            0 0.593       0.593           0
    #> racehispanic  0.114       -0.114            0 0.027       0.027           0
    #> racewhite    -2.091        2.091            0 0.620       0.620           0
    #> married      -1.729        1.729            0 0.677       0.677           0
    #> nodegree      0.886       -0.886            0 0.403       0.403           0
    #> re74         -3.547        3.547            0 0.729       0.729           0
    #> re75         -5.446        5.446            0 0.774       0.774           0
    #> 
    #> Summary of Balance for Weighted Data:
    #>              SMD TSMD Control TSMD Treated    KS TKS Control TKS Treated
    #> age            0        0.060        0.060 0.127       0.147       0.026
    #> education      0        0.024        0.024 0.081       0.083       0.023
    #> raceblack      0       -0.045       -0.045 0.000       0.017       0.017
    #> racehispanic   0        0.008        0.008 0.000       0.002       0.002
    #> racewhite      0        0.049        0.049 0.000       0.015       0.015
    #> married        0        0.135        0.135 0.000       0.053       0.053
    #> nodegree       0       -0.051       -0.051 0.000       0.023       0.023
    #> re74           0        0.044        0.044 0.578       0.590       0.015
    #> re75           0        0.049        0.049 0.566       0.577       0.015
    #> 
    #> Effective Sample Sizes:
    #>          Control Treated
    #> All       2490.    185. 
    #> Weighted   367.3   180.6

The `SMD` column contains the standardized mean differences for each
covariate between the treatment groups; after weighting, the values are
all equal to zero because of the properties of the implied regression
weights. However, the difference between each treatment group and target
sample remain, as displayed in the `TSMD Treated` and `TSMD Control`
columns, which contain the standardized mean differences between each
treatment group and the target sample (which in this case is the treated
group because the ATT was requested).

We can summarize this balance table in a plot:

``` r
plot(s)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

We can examine the distribution of weights to see to what degree
negative weights are present:

``` r
plot(lmw.out, type = "weights")
```

<img src="man/figures/README-unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

A red line indicates the average of the weights in that group (the
weights are scaled so this is always equal to 1). Negative weights are
present in the control group. Weights are also highly variable in the
control group, leading to a decreased effective sample size (ESS).

We can further examine extrapolation for specific covariates:

``` r
plot(lmw.out, type = "extrapolation",
     variables = ~married + re75)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

The $\times$ indicates the mean of the covariate in the target
population (the treated group), and the vertical line indicates the mean
of the covariate weighted by the implied regression weights. For these
covariates, the implied regression weights yield a sample fairly
representative of the target population.

We can examine how influential individual points are using their sample
influence curves, which are a function of their residuals, leverages,
and implied weights:

``` r
plot(lmw.out, type = "influence", outcome = re78)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

Finally, we can fit the outcome model and extract the average treatment
effect estimates:

``` r
lmw.fit <- lmw_est(lmw.out, outcome = re78)

summary(lmw.fit)
```

    #> 
    #> Effect estimates:
    #>              Estimate Std. Error 95% CI L 95% CI U t value Pr(>|t|)
    #> E[Y₁-Y₀|A=1]    751.9      788.9   -795.0   2298.9   0.953    0.341
    #> 
    #> Residual standard error: 10070 on 2665 degrees of freedom

In addition, `lmw` can be used with multi-category treatments, two-stage
least squares estimation of instrumental variable models, and
doubly-robust estimators by interfacing with the `MatchIt` and
`WeightIt` packages to implement these diagnostics for regression in a
matched or weighted sample.

### References

Chattopadhyay, A., & Zubizarreta, J. R. (2023). On the implied weights
of linear regression for causal inference. *Biometrika*, 110(3),
615–629. <https://doi.org/10.1093/biomet/asac058>
