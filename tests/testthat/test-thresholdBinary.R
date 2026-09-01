test_that("default threshold and controls work with one covariate", {
  calls <- list()
  cv_control <- list(V = 2, stratifyCV = TRUE)

  local_mocked_bindings(
    tmleThreshold.auto = function(threshold, W, A, Y, Delta, weights,
                                  sl_library, method, cvControl, ...) {
      calls[[length(calls) + 1L]] <<- list(
        threshold = threshold,
        W = W,
        A = A,
        Y = Y,
        Delta = Delta,
        weights = weights,
        sl_library = sl_library,
        method = method,
        cvControl = cvControl
      )
      list(
        threshold = threshold,
        estimate = 0.2,
        se = 0.05,
        ci_lo = 0.102,
        ci_hi = 0.298
      )
    },
    .package = "npthreshold"
  )

  data <- data.frame(
    adjustment = c(-1, 0, 1, 2),
    response = c(0, 1, 0, 1),
    biomarker_value = c(NA, 4, 2, 3),
    analysis_weight = c(0, 1, 1, 1)
  )

  result <- thresholdBinary(
    data = data,
    covariates = "adjustment",
    outcome = "response",
    marker = "biomarker_value",
    weights = "analysis_weight",
    method = "custom.method",
    cvControl = cv_control
  )

  expect_equal(result$threshold, 2)
  expect_equal(result$n_in_bin, 3)
  expect_equal(result$n_events_in_bin, 2)
  expect_equal(ncol(calls[[1]]$W), 1)
  expect_equal(calls[[1]]$A, c(0, 4, 2, 3))
  expect_identical(calls[[1]]$method, "custom.method")
  expect_identical(calls[[1]]$cvControl, cv_control)
})

test_that("missing outcomes are allowed only when Delta is zero", {
  captured_y <- NULL

  local_mocked_bindings(
    tmleThreshold.auto = function(threshold, W, A, Y, ...) {
      captured_y <<- Y
      list(
        threshold = threshold,
        estimate = 0.2,
        se = 0.05,
        ci_lo = 0.102,
        ci_hi = 0.298
      )
    },
    .package = "npthreshold"
  )

  data <- data.frame(
    adjustment = 1:4,
    response = c(1, NA, 0, 1),
    biomarker_value = 1:4,
    observed = c(1, 0, 1, 1)
  )

  result <- thresholdBinary(
    data = data,
    covariates = "adjustment",
    outcome = "response",
    marker = "biomarker_value",
    Delta = "observed",
    threshold_list = 1
  )

  expect_equal(captured_y, c(1, 0, 0, 1))
  expect_equal(result$n_events_in_bin, 2)

  data$observed[2] <- 1
  expect_error(
    thresholdBinary(
      data,
      covariates = "adjustment",
      outcome = "response",
      marker = "biomarker_value",
      Delta = "observed",
      threshold_list = 1
    ),
    "may be missing only where `Delta` is 0",
    fixed = TRUE
  )
})

test_that("invalid inputs fail before fitting", {
  data <- data.frame(
    adjustment = 1:4,
    response = c(0, 1, 0, 1),
    biomarker_value = 1:4,
    analysis_weight = c(1, 1, 1, -1)
  )

  expect_error(
    thresholdBinary(
      data,
      covariates = "missing_column",
      outcome = "response",
      marker = "biomarker_value",
      threshold_list = 1
    ),
    "missing from `data`",
    fixed = TRUE
  )

  expect_error(
    thresholdBinary(
      data,
      covariates = "adjustment",
      outcome = "response",
      marker = "biomarker_value",
      weights = "analysis_weight",
      threshold_list = 1
    ),
    "finite, nonnegative",
    fixed = TRUE
  )

  data$analysis_weight <- 1
  expect_error(
    thresholdBinary(
      data,
      covariates = "adjustment",
      outcome = "response",
      marker = "biomarker_value",
      threshold_list = 10
    ),
    "No positive-weight observed outcomes",
    fixed = TRUE
  )
})

test_that("Super Learner fits receive weights and full cvControl", {
  calls <- list()
  cv_control <- list(V = 2, stratifyCV = TRUE)

  local_mocked_bindings(
    SuperLearner = function(Y, X, newX = NULL, family, SL.library,
                            method, cvControl, obsWeights, ...) {
      calls[[length(calls) + 1L]] <<- list(
        Y = Y,
        X = X,
        newX = newX,
        cvControl = cvControl,
        obsWeights = obsWeights
      )
      prediction_length <- if (is.null(newX)) nrow(X) else nrow(newX)
      list(SL.predict = rep(0.5, prediction_length))
    },
    .package = "SuperLearner"
  )
  local_mocked_bindings(
    tmle.inefficient = function(...) {
      list(psi = 0.2, EIF = rep(0, 4))
    },
    .package = "npthreshold"
  )

  npthreshold:::tmleThreshold.auto(
    threshold = 1,
    W = matrix(1:4, ncol = 1),
    A = c(0, 1, 2, 3),
    Y = c(0, 1, 0, 1),
    Delta = c(1, 0, 1, 1),
    weights = c(1, 2, 0, 4),
    sl_library = "SL.mean",
    method = "custom.method",
    cvControl = cv_control
  )

  expect_length(calls, 3)
  expect_equal(calls[[1]]$obsWeights, 4)
  expect_equal(calls[[2]]$obsWeights, c(1, 2, 0, 4))
  expect_equal(calls[[3]]$obsWeights, c(1, 2, 0, 4))
  expect_true(all(vapply(calls, function(x) {
    identical(x$cvControl, cv_control)
  }, logical(1))))
})

test_that("adaptive truncation uses the binary threshold indicator", {
  truncation_inputs <- list()

  local_mocked_bindings(
    truncate_pscore_adaptive = function(A, pi, ...) {
      truncation_inputs[[length(truncation_inputs) + 1L]] <<- A
      pi
    },
    .package = "npthreshold"
  )

  suppressWarnings(
    npthreshold:::tmle.inefficient(
      threshold = 0,
      A = c(-2, -1, 1, 2),
      Delta = rep(1, 4),
      Y = c(0, 1, 0, 1),
      gv = rep(0.5, 4),
      Qv = rep(0.5, 4),
      Gv = rep(1, 4),
      weights = rep(1, 4)
    )
  )

  expect_equal(truncation_inputs[[1]], c(0, 0, 1, 1))
  expect_equal(truncation_inputs[[2]], rep(1, 4))
})

test_that("fit errors identify the threshold", {
  local_mocked_bindings(
    tmleThreshold.auto = function(...) stop("learner failed"),
    .package = "npthreshold"
  )

  data <- data.frame(
    adjustment = 1:4,
    response = c(0, 1, 0, 1),
    biomarker_value = 1:4
  )

  expect_error(
    thresholdBinary(
      data,
      covariates = "adjustment",
      outcome = "response",
      marker = "biomarker_value",
      threshold_list = 1
    ),
    "thresholdBinary fit failed at threshold 1: learner failed",
    fixed = TRUE
  )
})

test_that("a basic Super Learner fit completes", {
  set.seed(1)
  n <- 60
  data <- data.frame(
    adjustment = seq(-1, 1, length.out = n),
    response = rep(c(0, 1, 0), length.out = n),
    biomarker_value = seq(-2, 2, length.out = n)
  )

  result <- thresholdBinary(
    data = data,
    covariates = "adjustment",
    outcome = "response",
    marker = "biomarker_value",
    threshold_list = c(-0.5, 0.5),
    sl_library = "SL.mean",
    method = "method.NNLS",
    cvControl = list(V = 2)
  )

  expect_equal(nrow(result), 2)
  expect_named(
    result,
    c(
      "threshold", "estimate", "se", "ci_lo", "ci_hi",
      "n_in_bin", "n_events_in_bin", "estimate_monotone",
      "ci_lo_monotone", "ci_hi_monotone"
    )
  )
  expect_true(all(is.finite(result$estimate)))
  expect_true(all(is.finite(result$se)))
})
