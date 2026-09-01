mock_survival_result <- function(target_failure_time, ...) {
  data.table::data.table(
    treatment = 1,
    event_type = 1,
    times = list(target_failure_time),
    estimates = list(0.2),
    se = list(0.05),
    EIF = list(0),
    CI = list(matrix(c(0.102, 0.298), nrow = 1))
  )
}

test_that("events after tf are assigned to a later time bin", {
  captured <- NULL

  local_mocked_bindings(
    survtmle3_discrete = function(failure_time, target_failure_time,
                                  cross_fit, nfolds, ...) {
      captured <<- list(
        failure_time = failure_time,
        target_failure_time = target_failure_time,
        cross_fit = cross_fit,
        nfolds = nfolds
      )
      mock_survival_result(target_failure_time)
    },
    .package = "npthreshold"
  )

  data <- data.frame(
    adjustment = 1:4,
    time = c(12, 15, 17, 18),
    event = c(1, 1, 1, 0),
    biomarker = 1:4
  )

  thresholdSurv(
    data = data,
    covariates = "adjustment",
    failure_time = "time",
    event_type = "event",
    marker = "biomarker",
    tf = 15,
    threshold_list = 1,
    nbins_time = 2,
    learner.treatment = NULL,
    learner.event_type = NULL,
    learner.failure_time = NULL,
    learner.censoring_time = NULL,
    cross_fit = TRUE,
    nfolds = 3
  )

  expect_equal(captured$failure_time, c(1, 1, 2, 2))
  expect_equal(captured$target_failure_time, 1)
  expect_true(captured$cross_fit)
  expect_equal(captured$nfolds, 3)
})

test_that("events after an earliest-time tf are assigned to a later bin", {
  captured_failure_time <- NULL

  local_mocked_bindings(
    survtmle3_discrete = function(failure_time, target_failure_time, ...) {
      captured_failure_time <<- failure_time
      mock_survival_result(target_failure_time)
    },
    .package = "npthreshold"
  )

  thresholdSurv(
    data = data.frame(
      adjustment = 1:3,
      time = 1:3,
      event = c(1, 1, 1),
      biomarker = 1:3
    ),
    covariates = "adjustment",
    failure_time = "time",
    event_type = "event",
    marker = "biomarker",
    tf = 1,
    threshold_list = 1,
    nbins_time = 1,
    learner.treatment = NULL,
    learner.event_type = NULL,
    learner.failure_time = NULL,
    learner.censoring_time = NULL
  )

  expect_equal(captured_failure_time, c(1, 2, 2))
})

test_that("zero-weight rows are excluded from summaries", {
  captured_weights <- NULL

  local_mocked_bindings(
    survtmle3_discrete = function(weights, target_failure_time, ...) {
      captured_weights <<- weights
      mock_survival_result(target_failure_time)
    },
    .package = "npthreshold"
  )

  data <- data.frame(
    adjustment = 1:3,
    time = 1:3,
    event = c(0, 1, 1),
    biomarker = c(NA, 2, 3),
    analysis_weight = c(0, 1, 1)
  )

  result <- thresholdSurv(
    data = data,
    covariates = "adjustment",
    failure_time = "time",
    event_type = "event",
    marker = "biomarker",
    tf = 2,
    weights = "analysis_weight",
    learner.treatment = NULL,
    learner.event_type = NULL,
    learner.failure_time = NULL,
    learner.censoring_time = NULL
  )

  expect_equal(captured_weights, c(0, 1, 1))
  expect_equal(result$threshold, 2)
  expect_equal(result$rcdf, 1)
  expect_equal(result$n_in_bin, 2)
  expect_equal(result$n_events_in_bin, 1)
  expect_equal(result$person_time_in_bin, 4)
  expect_true(is.na(data$biomarker[1]))
})

test_that("left-point replacement uses the smallest threshold and keeps its CI", {
  local_mocked_bindings(
    survtmle3_discrete = mock_survival_result,
    .package = "npthreshold"
  )

  data <- data.frame(
    adjustment = 1:4,
    time = 1:4,
    event = c(0, 1, 0, 1),
    biomarker = 1:4
  )

  result <- thresholdSurv(
    data = data,
    covariates = "adjustment",
    failure_time = "time",
    event_type = "event",
    marker = "biomarker",
    tf = 3,
    threshold_list = c(2, 1),
    learner.treatment = NULL,
    learner.event_type = NULL,
    learner.failure_time = NULL,
    learner.censoring_time = NULL,
    modify_left_point = TRUE,
    left_estimate = 0.8,
    left_se = 0.02
  )

  left <- result[result$threshold == 1, ]
  right <- result[result$threshold == 2, ]

  expect_equal(left$estimate, 0.8)
  expect_equal(left$se, 0.02)
  expect_equal(left$ci_lo, 0.8 - 1.96 * 0.02)
  expect_equal(left$ci_hi, 0.8 + 1.96 * 0.02)
  expect_equal(right$estimate, 0.2)
  expect_equal(right$ci_lo, 0.102)
  expect_equal(right$ci_hi, 0.298)
})

test_that("tf outside observed follow-up fails before fitting", {
  data <- data.frame(
    adjustment = 1:3,
    time = 1:3,
    event = c(0, 1, 1),
    biomarker = 1:3
  )

  expect_error(
    thresholdSurv(
      data = data,
      covariates = "adjustment",
      failure_time = "time",
      event_type = "event",
      marker = "biomarker",
      tf = 4,
      threshold_list = 1,
      learner.treatment = NULL,
      learner.event_type = NULL,
      learner.failure_time = NULL,
      learner.censoring_time = NULL
    ),
    "`tf` must fall within the observed failure-time range.",
    fixed = TRUE
  )
})

test_that("survival fitting works without censoring", {
  n <- 80
  covariates <- matrix(
    seq(-1, 1, length.out = n),
    ncol = 1,
    dimnames = list(NULL, "W1")
  )
  learner <- sl3::Lrnr_glm$new()

  result <- survtmle3_discrete(
    failure_time = rep(1:4, length.out = n),
    event_type = rep(1, n),
    treatment = rep(c(0, 1), length.out = n),
    covariates = covariates,
    learner.treatment = learner,
    learner.failure_time = learner,
    learner.censoring_time = learner,
    learner.event_type = learner,
    target_failure_time = 3,
    target_treatment = 1,
    target_event_type = 1,
    cross_fit = FALSE,
    calibrate = FALSE,
    verbose = FALSE
  )

  expect_equal(nrow(result), 1)
  expect_true(is.finite(unlist(result$estimates)))
})

test_that("competing-risk fitting excludes censored outcomes", {
  n <- 120
  covariates <- matrix(
    seq(-1, 1, length.out = n),
    ncol = 1,
    dimnames = list(NULL, "W1")
  )
  learner <- sl3::Lrnr_glm$new()

  result <- survtmle3_discrete(
    failure_time = rep(1:4, length.out = n),
    event_type = rep(c(0, 1, 2, 1, 2), length.out = n),
    treatment = rep(c(0, 1), length.out = n),
    covariates = covariates,
    learner.treatment = learner,
    learner.failure_time = learner,
    learner.censoring_time = learner,
    learner.event_type = learner,
    target_failure_time = 3,
    target_treatment = 1,
    target_event_type = 1,
    cross_fit = FALSE,
    calibrate = FALSE,
    verbose = FALSE
  )

  expect_equal(nrow(result), 1)
  expect_true(is.finite(unlist(result$estimates)))
})

test_that("non-converged targeting is reported as an error", {
  n <- 80
  covariates <- matrix(
    seq(-1, 1, length.out = n),
    ncol = 1,
    dimnames = list(NULL, "W1")
  )
  learner <- sl3::Lrnr_glm$new()

  expect_error(
    survtmle3_discrete(
      failure_time = rep(1:4, length.out = n),
      event_type = rep(c(1, 0, 1, 1), length.out = n),
      treatment = rep(c(0, 1), length.out = n),
      covariates = covariates,
      learner.treatment = learner,
      learner.failure_time = learner,
      learner.censoring_time = learner,
      learner.event_type = learner,
      target_failure_time = 3,
      target_treatment = 1,
      target_event_type = 1,
      cross_fit = FALSE,
      calibrate = FALSE,
      max_iter = 1,
      tol = 0,
      verbose = FALSE
    ),
    "TMLE targeting did not converge after 1 iterations",
    fixed = TRUE
  )
})
