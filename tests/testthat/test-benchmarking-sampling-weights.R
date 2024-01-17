test_that("URI, binary treatment, sampling weights", {
  skip_if_not_installed("estimatr")

  data("lalonde"); lalonde$sw <- runif(nrow(lalonde))

  expect_no_condition(
    l <- lmw(re78 ~ treat + age + education + race + married + re74 + re75,
             data = lalonde, method = "URI", treat = "treat",
             s.weights = sw)
  )

  ### Weighted difference in means
  est <- weighted_mean_diff(lalonde$re78, lalonde$treat, l$weights)

  ### lmw_est() output
  expect_no_condition(
    e <- lmw_est(l)
  )

  expect_no_condition(
    s <- summary(e)
  )

  expect_equal_est(s$coefficients[1, "Estimate"],
                   est)

  ### estimatr::lm_robust() output
  f <- estimatr::lm_robust(re78 ~ treat + age + education + race + married + re74 + re75,
                           data = lalonde, se_type = "HC3", weights = sw)

  expect_equal_est(s$coefficients[1, c("Estimate", "Std. Error", "95% CI L", "95% CI U", "t value", "Pr(>|t|)")],
                   summary(f)$coefficients["treat", c("Estimate", "Std. Error", "CI Lower", "CI Upper", "t value", "Pr(>|t|)")])
})

test_that("MRI, binary treatment, ATE, sampling weights", {
  skip_if_not_installed("marginaleffects")
  skip_if_not_installed("estimatr")

  data("lalonde"); lalonde$sw <- runif(nrow(lalonde))

  expect_no_condition(
    l <- lmw(re78 ~ treat + age + education + race + married + re74 + re75,
             data = lalonde, method = "MRI", treat = "treat",
             estimand = "ATE", s.weights = sw)
  )

  ### Weighted difference in means
  est <- weighted_mean_diff(lalonde$re78, lalonde$treat, l$weights)

  ### lmw_est() output
  expect_no_condition(
    e <- lmw_est(l)
  )

  expect_no_condition(
    s <- summary(e)
  )

  expect_equal_est(s$coefficients[1, "Estimate"],
                   est)

  ### avg_comparisons() output
  f <- lm(re78 ~ treat * (age + education + race + married + re74 + re75),
          data = lalonde, weights = sw)

  ac <- marginaleffects::avg_comparisons(f, variables = "treat",
                                         vcov = "HC3", wts = "sw")

  expect_equal_est(s$coefficients[1, c("Estimate", "Std. Error", "t value")],
                   unlist(ac[1, c("estimate", "std.error", "statistic")]))

  ### estimatr output
  fl <- estimatr::lm_lin(re78 ~ treat,
                         ~ age + education + race + married + re74 + re75,
                         data = lalonde, se_type = "HC3", weights = sw)

  expect_equal_est(s$coefficients[1, c("Estimate", "Std. Error", "95% CI L", "95% CI U", "t value", "Pr(>|t|)")],
                   summary(fl)$coefficients["treat", c("Estimate", "Std. Error", "CI Lower", "CI Upper", "t value", "Pr(>|t|)")])
})

test_that("MRI, binary treatment, ATT, sampling weights", {
  skip_if_not_installed("marginaleffects")

  data("lalonde"); lalonde$sw <- runif(nrow(lalonde))

  expect_no_condition(
    l <- lmw(re78 ~ treat + age + education + race + married + re74 + re75,
             data = lalonde, method = "MRI", treat = "treat",
             estimand = "ATT",
             s.weights = sw)
  )

  ### Weighted difference in means
  est <- weighted_mean_diff(lalonde$re78, lalonde$treat, l$weights)

  ### lmw_est() output
  expect_no_condition(
    e <- lmw_est(l)
  )

  expect_no_condition(
    s <- summary(e)
  )

  expect_equal_est(s$coefficients[1, "Estimate"],
                   est)

  ### avg_comparisons() output
  f <- lm(re78 ~ treat * (age + education + race + married + re74 + re75),
          data = lalonde, weights = sw)

  ac <- marginaleffects::avg_comparisons(f, variables = "treat",
                                         vcov = "HC3",
                                         newdata = subset(lalonde, treat == 1),
                                         wts = "sw")

  expect_equal_est(s$coefficients[1, c("Estimate", "Std. Error", "t value")],
                   unlist(ac[1, c("estimate", "std.error", "statistic")]))
})

test_that("URI, binary treatment, fixed effects, sampling weights", {
  skip_if_not_installed("estimatr")

  data("lalonde"); lalonde$sw <- runif(nrow(lalonde))

  expect_no_condition(
    l <- lmw(re78 ~ treat + age + education + married + re74 + re75,
             data = lalonde, method = "URI", treat = "treat",
             fixef = ~race,
             s.weights = sw)
  )

  ### Weighted difference in means
  est <- weighted_mean_diff(lalonde$re78, lalonde$treat, l$weights)

  ### lmw_est() output
  expect_no_condition(
    e <- lmw_est(l)
  )

  expect_no_condition(
    s <- summary(e)
  )

  expect_equal_est(s$coefficients[1, "Estimate"],
                   est)

  ### estimatr::lm_robust() output
  f <- estimatr::lm_robust(re78 ~ treat + age + education + married + re74 + re75,
                           data = lalonde, fixed_effects = ~race,
                           se_type = "HC3", weights = sw)

  expect_equal_est(s$coefficients[1, c("Estimate", "Std. Error", "95% CI L", "95% CI U", "t value", "Pr(>|t|)")],
                   summary(f)$coefficients["treat", c("Estimate", "Std. Error", "CI Lower", "CI Upper", "t value", "Pr(>|t|)")])
})

test_that("MRI, binary treatment, ATE, fixed effects, sampling weights", {
  skip_if_not_installed("marginaleffects")

  data("lalonde"); lalonde$sw <- runif(nrow(lalonde))

  expect_no_condition(
    l <- lmw(re78 ~ treat + age + education + married + re74 + re75,
             data = lalonde, method = "MRI", treat = "treat",
             fixef = ~race,
             s.weights = sw)
  )

  ### Weighted difference in means
  est <- weighted_mean_diff(lalonde$re78, lalonde$treat, l$weights)

  ### lmw_est() output
  expect_no_condition(
    e <- lmw_est(l)
  )

  expect_no_condition(
    s <- summary(e)
  )

  expect_equal_est(s$coefficients[1, "Estimate"],
                   est)

  ### URI + covariate version
  l2 <- lmw(re78 ~ treat * (age + education + married + re74 + re75) + race,
            data = lalonde, method = "URI", treat = "treat",
            s.weights = sw)

  est2 <- weighted_mean_diff(lalonde$re78, lalonde$treat, l2$weights)

  expect_equal_est(est, est2)

  ### avg_comparisons() output
  f <- lm(re78 ~ treat * (age + education + married + re74 + re75) + race,
          data = lalonde, weights = sw)

  ac <- marginaleffects::avg_comparisons(f, variables = "treat",
                                         vcov = "HC3", wts = "sw")

  expect_equal_est(s$coefficients[1, c("Estimate", "Std. Error", "t value")],
                   unlist(ac[1, c("estimate", "std.error", "statistic")]))
})

test_that("URI, binary treatment, 2SLS, sampling weights", {
  skip_if_not_installed("estimatr")

  data("lalonde"); lalonde$sw <- runif(nrow(lalonde))

  expect_no_condition(
    l <- lmw_iv(re78 ~ treat + age + education + race + married + re74 + re75,
                data = lalonde, method = "URI", treat = "treat",
                iv = ~Ins,
                s.weights = sw)
  )

  ### Weighted difference in means
  est <- weighted_mean_diff(lalonde$re78, lalonde$treat, l$weights)

  ### lmw_est() output
  expect_no_condition(
    e <- lmw_est(l, robust = "HC1")
  )

  expect_no_condition(
    s <- summary(e)
  )

  expect_equal_est(s$coefficients[1, "Estimate"],
                   est)

  ### estimatr::lm_robust() output
  f <- estimatr::iv_robust(re78 ~ treat + age + education + race + married + re74 + re75 |
                             Ins + age + education + race + married + re74 + re75,
                           data = lalonde, se_type = "HC1", weights = sw)

  expect_equal_est(s$coefficients[1, c("Estimate", "Std. Error", "95% CI L", "95% CI U", "t value", "Pr(>|t|)")],
                   summary(f)$coefficients["treat", c("Estimate", "Std. Error", "CI Lower", "CI Upper", "t value", "Pr(>|t|)")])
})

test_that("URI, binary treatment, 2SLS, fixed effects, CR SEs, sampling weights", {
  skip_if_not_installed("estimatr")

  data("lalonde"); lalonde$sw <- runif(nrow(lalonde))

  expect_no_condition(
    l <- lmw_iv(re78 ~ treat + age + education + married + re74 + re75,
                data = lalonde, method = "URI", treat = "treat",
                iv = ~Ins, fixef = ~race,
                s.weights = sw)
  )

  ### Weighted difference in means
  est <- weighted_mean_diff(lalonde$re78, lalonde$treat, l$weights)

  ### lmw_est() output
  expect_no_condition(
    e <- lmw_est(l, robust = "HC1", cluster = ~race)
  )

  expect_no_condition(
    s <- summary(e)
  )

  expect_equal_est(s$coefficients[1, "Estimate"],
                   est)

  ### estimatr::lm_robust() output
  f <- estimatr::iv_robust(re78 ~ treat + age + education + married + re74 + re75 |
                             Ins + age + education + married + re74 + re75,
                           data = lalonde, fixed_effects = ~race,
                           se_type = "stata", cluster = race,
                           weights = sw)

  expect_equal_est(s$coefficients[1, c("Estimate", "Std. Error", "95% CI L", "95% CI U", "t value", "Pr(>|t|)")],
                   summary(f)$coefficients["treat", c("Estimate", "Std. Error", "CI Lower", "CI Upper", "t value", "Pr(>|t|)")])
})

test_that("MRI, binary treatment, 2SLS, sampling weights", {
  skip_if_not_installed("estimatr")

  data("lalonde"); lalonde$sw <- runif(nrow(lalonde))

  expect_no_condition(
    l <- lmw_iv(re78 ~ treat + age + education + married + re74 + re75,
                data = lalonde, method = "MRI", treat = "treat",
                iv = ~Ins,
                s.weights = sw)
  )

  ### Weighted difference in means
  est <- weighted_mean_diff(lalonde$re78, lalonde$treat, l$weights)

  ### lmw_est() output
  expect_no_condition(
    e <- lmw_est(l, robust = "HC1")
  )

  expect_no_condition(
    s <- summary(e)
  )

  expect_equal_est(s$coefficients[1, "Estimate"],
                   est)

  ### estimatr::lm_robust() output
  f <- estimatr::iv_robust(re78 ~ treat * (age + (education) + (married) + (re74) + (re75)) |
                             Ins * ((age) + (education) + (married) + (re74) + (re75)),
                           data = transform(lalonde,
                                            age = age - weighted.mean(age, sw),
                                            education = education - weighted.mean(education, sw),
                                            married = married - weighted.mean(married, sw),
                                            re74 = re74 - weighted.mean(re74, sw),
                                            re75 = re75 - weighted.mean(re75, sw)),
                           se_type = "HC1", weights = sw)

  expect_equal_est(s$coefficients[1, c("Estimate", "Std. Error", "95% CI L", "95% CI U", "t value", "Pr(>|t|)")],
                   summary(f)$coefficients["treat", c("Estimate", "Std. Error", "CI Lower", "CI Upper", "t value", "Pr(>|t|)")])
})

test_that("URI, multi-category treatment, sampling weights", {
  skip_if_not_installed("estimatr")

  data("lalonde"); lalonde$sw <- runif(nrow(lalonde))

  expect_no_condition(
    l <- lmw(re78 ~ treat_multi + age + education + race + married + re74 + re75,
             data = lalonde, method = "URI", treat = "treat_multi",
             contrast = c("3", "1"),
             s.weights = sw)
  )

  ### Weighted difference in means
  est <- weighted_mean_diff(lalonde$re78, lalonde$treat_multi == "3", l$weights)

  ### lmw_est() output
  expect_no_condition(
    e <- lmw_est(l)
  )

  expect_no_condition(
    s <- summary(e)
  )

  expect_equal_est(s$coefficients["E[Y3-Y1]", "Estimate"],
                   est)

  ### estimatr::lm_robust() output
  f <- estimatr::lm_robust(re78 ~ treat_multi + age + education + race + married + re74 + re75,
                           data = lalonde, se_type = "HC3", weights = sw)

  expect_equal_est(s$coefficients["E[Y2-Y1]", c("Estimate", "Std. Error", "95% CI L", "95% CI U", "t value", "Pr(>|t|)")],
                   summary(f)$coefficients["treat_multi2", c("Estimate", "Std. Error", "CI Lower", "CI Upper", "t value", "Pr(>|t|)")])

  expect_equal_est(s$coefficients["E[Y3-Y1]", c("Estimate", "Std. Error", "95% CI L", "95% CI U", "t value", "Pr(>|t|)")],
                   summary(f)$coefficients["treat_multi3", c("Estimate", "Std. Error", "CI Lower", "CI Upper", "t value", "Pr(>|t|)")])
})

test_that("MRI, multi-category treatment, ATE, sampling weights", {
  skip_if_not_installed("estimatr")

  data("lalonde"); lalonde$sw <- runif(nrow(lalonde))

  expect_no_condition(
    l <- lmw(re78 ~ treat_multi + age + education + race + married + re74 + re75,
             data = lalonde, method = "MRI", treat = "treat_multi",
             s.weights = sw)
  )

  ### Weighted difference in means
  est <- weighted_mean_diff(lalonde$re78, lalonde$treat_multi == "3", l$weights,
                            subset = lalonde$treat_multi %in% c("2", "3"))

  ### lmw_est() output
  expect_no_condition(
    e <- lmw_est(l)
  )

  expect_no_condition(
    s <- summary(e)
  )

  expect_equal_est(s$coefficients["E[Y3-Y2]", "Estimate"],
                   est)

  ### estimatr::lm_robust() output
  f <- estimatr::lm_lin(re78 ~ treat_multi,
                        ~ age + education + race + married + re74 + re75,
                        data = lalonde, se_type = "HC3", weights = sw)

  expect_equal_est(s$coefficients["E[Y2-Y1]", c("Estimate", "Std. Error", "95% CI L", "95% CI U", "t value", "Pr(>|t|)")],
                   summary(f)$coefficients["treat_multi2", c("Estimate", "Std. Error", "CI Lower", "CI Upper", "t value", "Pr(>|t|)")])

  expect_equal_est(s$coefficients["E[Y3-Y1]", c("Estimate", "Std. Error", "95% CI L", "95% CI U", "t value", "Pr(>|t|)")],
                   summary(f)$coefficients["treat_multi3", c("Estimate", "Std. Error", "CI Lower", "CI Upper", "t value", "Pr(>|t|)")])
})

test_that("MRI, multi-category treatment, ATT, sampling weights", {
  skip_if_not_installed("marginaleffects")

  data("lalonde"); lalonde$sw <- runif(nrow(lalonde))

  expect_no_condition(
    l <- lmw(re78 ~ treat_multi + age + education + race + married + re74 + re75,
             data = lalonde, method = "MRI", treat = "treat_multi",
             estimand = "ATT", focal = "1",
             s.weights = sw)
  )

  ### Weighted difference in means
  est <- weighted_mean_diff(lalonde$re78, lalonde$treat_multi == "1", l$weights,
                            subset = lalonde$treat_multi %in% c("1", "3"))

  ### lmw_est() output
  expect_no_condition(
    e <- lmw_est(l)
  )

  expect_no_condition(
    s <- summary(e)
  )

  expect_equal_est(s$coefficients["E[Y1-Y3]", "Estimate"],
                   est)

  ### avg_comparisons() output
  f <- lm(re78 ~ treat_multi * (age + education + race + married + re74 + re75),
          data = lalonde, weights = sw)

  ac <- marginaleffects::avg_predictions(f, variables = "treat_multi",
                                         vcov = "HC3",
                                         newdata = subset(lalonde, treat_multi == "1"),
                                         wts = "sw",
                                         hypothesis = "pairwise")

  expect_equal_est(s$coefficients[c("E[Y1-Y2]", "E[Y1-Y3]", "E[Y2-Y3]"), c("Estimate", "Std. Error", "t value")],
                   as.matrix(ac[1:3, c("estimate", "std.error", "statistic")]))
})

test_that("URI, multi-category treatment, fixed effects, sampling weights", {
  skip_if_not_installed("estimatr")

  data("lalonde"); lalonde$sw <- runif(nrow(lalonde))

  expect_no_condition(
    l <- lmw(re78 ~ treat_multi + age + education + married + re74 + re75,
             data = lalonde, method = "URI", treat = "treat_multi",
             contrast = c("3", "1"), fixef = ~race,
             s.weights = sw)
  )

  ### Weighted difference in means
  est <- weighted_mean_diff(lalonde$re78, lalonde$treat_multi == "3", l$weights)

  ### lmw_est() output
  expect_no_condition(
    e <- lmw_est(l)
  )

  expect_no_condition(
    s <- summary(e)
  )

  expect_equal_est(s$coefficients["E[Y3-Y1]", "Estimate"],
                   est)

  ### estimatr::lm_robust() output
  f <- estimatr::lm_robust(re78 ~ treat_multi + age + education + married + re74 + re75,
                           data = lalonde, se_type = "HC3", fixed_effects = ~race,
                           weights = sw)

  expect_equal_est(s$coefficients["E[Y2-Y1]", c("Estimate", "Std. Error", "95% CI L", "95% CI U", "t value", "Pr(>|t|)")],
                   summary(f)$coefficients["treat_multi2", c("Estimate", "Std. Error", "CI Lower", "CI Upper", "t value", "Pr(>|t|)")])

  expect_equal_est(s$coefficients["E[Y3-Y1]", c("Estimate", "Std. Error", "95% CI L", "95% CI U", "t value", "Pr(>|t|)")],
                   summary(f)$coefficients["treat_multi3", c("Estimate", "Std. Error", "CI Lower", "CI Upper", "t value", "Pr(>|t|)")])
})
