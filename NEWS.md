`lmw` News and Updates
======

# `lmw` 0.0.2

* Changed the `var` argument in `plot(., type = "extrapolation")` to `variables` (`var` still works as an abbreviation).

* Added a unit testing framework.

* Weights now have an average of 1 in each treatment group. This only affects the display of the weights, making it easier to assess extreme weights, but does not affect balance or any other property of the weights.

* When using cluster-robust standard errors in `lmw_est()`, the degrees of freedom for the corresponding t-testing as now computed as one less than the number of clusters to be consistent with other software.

* `"IPWRA"` is now an acceptable alias for `"WLS"` when supplied to the `dr.method` of `lmw()` with base weights. IPWRA is the name Stata gives to the same estimator.

* When a `weightit` or `matchit` object is supplied to `obj` in `lmw()` or `lmw_iv()`, the `focal` component of the supplied object is automatically supplied and does not need to be specified. A warning will be thrown if `focal` is supplied and differs from the `focal` component of the supplied object. 

* Errors are prettier.

# `lmw` 0.0.1

* First version!
