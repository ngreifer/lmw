skip_if_not_installed("MatchIt")

test_that("URI, binary treatment, MatchIt", {
  skip_if_not_installed("estimatr")

  data("lalonde")

  M <- MatchIt::matchit(treat ~ age + education + race + married + re74 + re75,
                        data = lalonde)
  md <- MatchIt::match.data(M, data = lalonde)

  expect_no_condition(
    l <- lmw(re78 ~ treat + age + education + race + married + re74 + re75,
             data = lalonde, method = "URI", treat = "treat",
             base.weights = M$weights)
  )

  ### Weighted difference in means
  est <- weighted_mean_diff(lalonde$re78, lalonde$treat, l$weights)

  ### lmw_est() output
  expect_no_condition(
    e <- lmw_est(l, cluster = ~M$subclass)
  )

  expect_no_condition(
    s <- summary(e)
  )

  expect_equal_est(s$coefficients[1, "Estimate"],
                   est)

  ### estimatr::lm_robust() output
  f <- estimatr::lm_robust(re78 ~ treat + age + education + race + married + re74 + re75,
                           data = md, cluster = subclass, weights = weights,
                           se_type = "stata")

  expect_equal_est(s$coefficients[1, c("Estimate", "Std. Error", "95% CI L", "95% CI U", "t value", "Pr(>|t|)")],
                   summary(f)$coefficients["treat", c("Estimate", "Std. Error", "CI Lower", "CI Upper", "t value", "Pr(>|t|)")])
})

test_that("MRI, binary treatment, ATT, MatchIt", {
  skip_if_not_installed("marginaleffects")

  data("lalonde")

  M <- MatchIt::matchit(treat ~ age + education + race + married + re74 + re75,
                        data = lalonde)
  md <- MatchIt::match.data(M, data = lalonde)

  expect_no_condition(
    l <- lmw(re78 ~ treat + age + education + race + married + re74 + re75,
             data = lalonde, method = "MRI", treat = "treat",
             estimand = "ATT",
             base.weights = M$weights)
  )

  ### Weighted difference in means
  est <- weighted_mean_diff(lalonde$re78, lalonde$treat, l$weights)

  ### lmw_est() output
  expect_no_condition(
    e <- lmw_est(l, cluster = ~M$subclass)
  )

  expect_no_condition(
    s <- summary(e)
  )

  expect_equal_est(s$coefficients[1, "Estimate"],
                   est)

  ### avg_comparisons() output
  f <- lm(re78 ~ treat * (age + education + race + married + re74 + re75),
          data = md, weights = weights)

  ac <- marginaleffects::avg_comparisons(f, variables = "treat",
                                         vcov = ~subclass,
                                         newdata = subset(md, treat == 1))

  expect_equal_est(s$coefficients[1, c("Estimate", "Std. Error", "t value")],
                   unlist(ac[1, c("estimate", "std.error", "statistic")]))
})

test_that("URI, binary treatment, fixed effects, MatchIt", {
  skip_if_not_installed("estimatr")

  data("lalonde")

  M <- MatchIt::matchit(treat ~ age + education + race + married + re74 + re75,
                        data = lalonde)
  md <- MatchIt::match.data(M, data = lalonde)

  expect_no_condition(
    l <- lmw(re78 ~ treat + age + education + married + re74 + re75,
             data = lalonde, method = "URI", treat = "treat",
             fixef = ~race,
             base.weights = M$weights)
  )

  ### Weighted difference in means
  est <- weighted_mean_diff(lalonde$re78, lalonde$treat, l$weights)

  ### lmw_est() output
  expect_no_condition(
    e <- lmw_est(l, cluster = ~M$subclass)
  )

  expect_no_condition(
    s <- summary(e)
  )

  expect_equal_est(s$coefficients[1, "Estimate"],
                   est)

  ### estimatr::lm_robust() output
  f <- estimatr::lm_robust(re78 ~ treat + age + education + married + re74 + re75,
                           data = md, fixed_effects = ~race, weights = weights,
                           se_type = "stata", cluster = subclass)

  expect_equal_est(s$coefficients[1, c("Estimate", "Std. Error", "95% CI L", "95% CI U", "t value", "Pr(>|t|)")],
                   summary(f)$coefficients["treat", c("Estimate", "Std. Error", "CI Lower", "CI Upper", "t value", "Pr(>|t|)")])
})

test_that("MRI, binary treatment, ATT, fixed effects, MatchIt", {
  skip_if_not_installed("marginaleffects")

  data("lalonde")

  M <- MatchIt::matchit(treat ~ age + education + race + married + re74 + re75,
                        data = lalonde)
  md <- MatchIt::match.data(M, data = lalonde)

  expect_no_condition(
    l <- lmw(re78 ~ treat + age + education + married + re74 + re75,
             data = lalonde, method = "MRI", treat = "treat",
             estimand = "ATT",
             fixef = ~race,
             base.weights = M$weights)
  )

  ### Weighted difference in means
  est <- weighted_mean_diff(lalonde$re78, lalonde$treat, l$weights)

  ### lmw_est() output
  expect_no_condition(
    e <- lmw_est(l, cluster = M$subclass)
  )

  expect_no_condition(
    s <- summary(e)
  )

  expect_equal_est(s$coefficients[1, "Estimate"],
                   est)

  ### URI + covariate version
  l2 <- lmw(re78 ~ treat * (age + education + married + re74 + re75) + race,
            data = lalonde, method = "URI", treat = "treat",
            estimand = "ATT",
            base.weights = M$weights)

  est2 <- weighted_mean_diff(lalonde$re78, lalonde$treat, l2$weights)

  expect_equal_est(est, est2)

  ### avg_comparisons() output
  f <- lm(re78 ~ treat * (age + education + married + re74 + re75) + race,
          data = md, weights = weights)

  ac <- marginaleffects::avg_comparisons(f, variables = "treat",
                                         vcov = ~subclass,
                                         newdata = subset(md, treat == 1))

  expect_equal_est(s$coefficients[1, c("Estimate", "Std. Error", "t value")],
                   unlist(ac[1, c("estimate", "std.error", "statistic")]))
})

test_that("URI, binary treatment, 2SLS, MatchIt", {
  skip_if_not_installed("estimatr")

  data("lalonde")

  M <- MatchIt::matchit(treat ~ age + education + race + married + re74 + re75,
                        data = lalonde)
  md <- MatchIt::match.data(M, data = lalonde)

  expect_no_condition(
    l <- lmw_iv(re78 ~ treat + age + education + race + married + re74 + re75,
                data = lalonde, method = "URI", treat = "treat",
                iv = ~Ins, base.weights = M$weights)
  )

  ### Weighted difference in means
  est <- weighted_mean_diff(lalonde$re78, lalonde$treat, l$weights)

  ### lmw_est() output
  expect_no_condition(
    e <- lmw_est(l, cluster = M$subclass)
  )

  expect_no_condition(
    s <- summary(e)
  )

  expect_equal_est(s$coefficients[1, "Estimate"],
                   est)

  ### estimatr::lm_robust() output
  f <- estimatr::iv_robust(re78 ~ treat + age + education + race + married + re74 + re75 |
                             Ins + age + education + race + married + re74 + re75,
                           data = md, se_type = "stata", cluster = subclass,
                           weights = weights)

  expect_equal_est(s$coefficients[1, c("Estimate", "Std. Error", "95% CI L", "95% CI U", "t value", "Pr(>|t|)")],
                   summary(f)$coefficients["treat", c("Estimate", "Std. Error", "CI Lower", "CI Upper", "t value", "Pr(>|t|)")])
})

test_that("URI, binary treatment, 2SLS, fixed effects, CR SEs, MatchIt", {
  skip_if_not_installed("ivreg")
  skip_if_not_installed("lmtest")

  data("lalonde")

  M <- MatchIt::matchit(treat ~ age + education + race + married + re74 + re75,
                        data = lalonde)
  md <- MatchIt::match.data(M, data = lalonde)

  expect_no_condition(
    l <- lmw_iv(re78 ~ treat + age + education + married + re74 + re75,
                data = lalonde, method = "URI", treat = "treat",
                iv = ~Ins, fixef = ~race,
                base.weights = M$weights)
  )

  ### Weighted difference in means
  est <- weighted_mean_diff(lalonde$re78, lalonde$treat, l$weights)

  ### lmw_est() output
  expect_no_condition(
    e <- lmw_est(l, cluster = ~race + M$subclass)
  )

  expect_no_condition(
    s <- summary(e)
  )

  expect_equal_est(s$coefficients[1, "Estimate"],
                   est)

  ### ivreg::ivreg() output; need for 2-way CR SEs
  f <- ivreg::ivreg(re78 ~ treat + age + education + married + race + re74 + re75 |
                             Ins + age + education + married + race + re74 + re75,
                           data = md, weights = weights)

  ss <- lmtest::coeftest(f, vcov. = sandwich::vcovCL, cluster = ~subclass + race,
                         type = "HC1", df = 2)

  expect_equal_est(s$coefficients[1, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")],
                   ss["treat", c("Estimate", "Std. Error", "t value", "Pr(>|t|)")])
})
