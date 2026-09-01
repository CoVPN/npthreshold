#' Non-parametric estimation of the covariate-adjusted threshold-response function
#'
#' This function estimates the covariate-adjusted probability of a binary outcome above
#' a range of thresholds for a continuous marker.
#'
#' @param data A \code{data.frame} containing the dataset to be analyzed.
#' @param covariates A \code{character} vector specifying the names of columns to be used as adjustment covariates.
#' @param outcome A \code{character} string specifying the name of the binary outcome column.
#' @param marker A \code{character} string specifying the name of the column representing the continuous marker.
#' @param Delta A \code{character} string specifying the name of the binary
#'   column indicating whether the outcome is observed. If \code{NULL}, all
#'   outcomes are assumed to be observed. When outcomes are missing, this
#'   binary-treatment estimator assumes outcome observation is independent of
#'   the outcome conditional on the covariates and the threshold indicator.
#' @param weights A \code{character} string specifying the name of a
#'   nonnegative observation-weight column. The weights are used in the Super
#'   Learner fits, targeting step, averaging, and influence function, and are
#'   not automatically normalized. Observations with zero weight are excluded
#'   from counts and estimation. If \code{NULL}, all observations receive
#'   weight one.
#' @param threshold_list A \code{numeric} vector of threshold values for estimation.
#' If \code{NULL}, estimation is performed using the minimum marker value among
#' observations with positive weight. A marker value is considered above a
#' threshold when it is greater than or equal to that threshold.
#' @param sl_library A \code{character} vector specifying the Super Learner libraries to be used for estimation.
#' The default is \code{c("SL.mean", "SL.glm")}.
#' @param method A \code{character} string specifying the method used by Super Learner. The default is \code{"method.CC_nloglik"}.
#' @param cvControl A list specifying the cross-validation control parameters for Super Learner. The default is \code{list(V = 5)}.
#'
#' @return A \code{data.frame} containing the following columns for each threshold in \code{threshold_list}:
#' \itemize{
#'   \item \code{threshold}: The threshold value.
#'   \item \code{estimate}: The estimated value at the threshold.
#'   \item \code{se}: The standard error of the estimate.
#'   \item \code{ci_lo}: The lower bound of the pointwise 95% confidence interval.
#'   \item \code{ci_hi}: The upper bound of the pointwise 95% confidence interval.
#'   \item \code{estimate_monotone}: The monotone-adjusted estimate.
#'   \item \code{ci_lo_monotone}: The lower bound of the monotone-adjusted pointwise 95% confidence interval.
#'   \item \code{ci_hi_monotone}: The upper bound of the monotone-adjusted pointwise 95% confidence interval.
#'   \item \code{n_in_bin}: The number of positive-weight observations at or above the threshold.
#'   \item \code{n_events_in_bin}: The number of observed events among positive-weight observations at or above the threshold.
#' }
#' Thresholds with fewer than 10 observed events receive very little weight in
#' the monotone projection.
#'
#' @import data.table isotone SuperLearner
#' @export
thresholdBinary <- function(data,
                            covariates,
                            outcome,
                            marker,
                            Delta = NULL,
                            weights = NULL,
                            threshold_list = NULL,
                            sl_library = c("SL.mean", "SL.glm"),
                            method = "method.CC_nloglik",
                            cvControl = list(V = 5)){

  data <- data.frame(data)

  if (!is.character(covariates) ||
      length(covariates) < 1L ||
      anyNA(covariates)) {
    stop("`covariates` must contain at least one column name.", call. = FALSE)
  }

  named_arguments <- list(outcome = outcome, marker = marker)
  if (!is.null(Delta)) named_arguments$Delta <- Delta
  if (!is.null(weights)) named_arguments$weights <- weights

  invalid_names <- vapply(
    named_arguments,
    function(x) !is.character(x) || length(x) != 1L || is.na(x),
    logical(1)
  )
  if (any(invalid_names)) {
    stop(
      sprintf(
        "The following arguments must each name one column: %s.",
        paste(names(named_arguments)[invalid_names], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  required_columns <- unique(c(covariates, unlist(named_arguments)))
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "The following columns are missing from `data`: %s.",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  n <- nrow(data)
  if (n < 2L) {
    stop("`data` must contain at least two observations.", call. = FALSE)
  }

  outcome_values <- data[[outcome]]
  marker_values <- data[[marker]]
  Delta_values <- if (is.null(Delta)) rep(1, n) else data[[Delta]]
  weights_values <- if (is.null(weights)) rep(1, n) else data[[weights]]

  if (!is.numeric(weights_values) ||
      any(!is.finite(weights_values)) ||
      any(weights_values < 0) ||
      sum(weights_values) <= 0) {
    stop("`weights` must be finite, nonnegative, and have a positive sum.",
         call. = FALSE)
  }

  is_binary <- function(x) {
    (is.numeric(x) || is.logical(x)) &&
      all(x[!is.na(x)] %in% c(0, 1))
  }

  if (!is_binary(Delta_values) || anyNA(Delta_values)) {
    stop("`Delta` must contain only 0 and 1, with no missing values.",
         call. = FALSE)
  }
  Delta_values <- as.numeric(Delta_values)

  if (!is_binary(outcome_values)) {
    stop("`outcome` must contain only 0, 1, or missing values.",
         call. = FALSE)
  }
  if (any(is.na(outcome_values) & Delta_values == 1)) {
    stop("`outcome` may be missing only where `Delta` is 0.", call. = FALSE)
  }
  outcome_values <- as.numeric(outcome_values)
  outcome_values[is.na(outcome_values)] <- 0

  if (!is.numeric(marker_values)) {
    stop("`marker` must be numeric.", call. = FALSE)
  }
  if (any(!is.finite(marker_values[weights_values > 0]))) {
    stop("`marker` must be finite for observations with positive weight.",
         call. = FALSE)
  }
  marker_values[weights_values == 0] <- 0

  covariate_data <- data[, covariates, drop = FALSE]
  if (anyNA(covariate_data)) {
    stop("Adjustment covariates must not contain missing values.",
         call. = FALSE)
  }
  numeric_covariates <- vapply(covariate_data, is.numeric, logical(1))
  if (any(vapply(covariate_data[numeric_covariates],
                 function(x) any(!is.finite(x)), logical(1)))) {
    stop("Numeric adjustment covariates must be finite.", call. = FALSE)
  }
  covariates_matrix <- model.matrix(~ . - 1, data = covariate_data)

  if (is.null(threshold_list)) {
    threshold_list <- min(marker_values[weights_values > 0])
  }
  if (!is.numeric(threshold_list) ||
      length(threshold_list) < 1L ||
      any(!is.finite(threshold_list))) {
    stop("`threshold_list` must contain at least one finite numeric value.",
         call. = FALSE)
  }

  if (!is.character(method) || length(method) != 1L || is.na(method)) {
    stop("`method` must be a single character string.", call. = FALSE)
  }
  if (!is.list(cvControl)) {
    stop("`cvControl` must be a list.", call. = FALSE)
  }

  observed_in_bin <- vapply(threshold_list, function(threshold) {
    sum(weights_values > 0 &
          Delta_values == 1 &
          marker_values >= threshold)
  }, numeric(1))
  if (any(observed_in_bin == 0)) {
    stop(
      sprintf(
        "No positive-weight observed outcomes are at or above threshold(s): %s.",
        paste(threshold_list[observed_in_bin == 0], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  results <- data.frame()

  for (threshold in threshold_list) {
    res_temp <- tryCatch(
      tmleThreshold.auto(threshold = threshold,
                         W = covariates_matrix,
                         A = marker_values,
                         Y = outcome_values,
                         Delta = Delta_values,
                         weights = weights_values,
                         sl_library = sl_library,
                         method = method,
                         cvControl = cvControl),
      error = function(e) {
        stop(
          sprintf(
            "thresholdBinary fit failed at threshold %s: %s",
            format(threshold),
            conditionMessage(e)
          ),
          call. = FALSE
        )
      }
    )

    results <- rbind(results, res_temp)
  }

  n_in_bin <- sapply(threshold_list, function(thresh) {
    sum(weights_values > 0 & marker_values >= thresh)
  })

  n_events_in_bin <- sapply(threshold_list, function(thresh) {
    sum(weights_values > 0 &
          marker_values >= thresh &
          Delta_values == 1 &
          outcome_values == 1)
  })

  results$n_in_bin <- n_in_bin
  results$n_events_in_bin <- n_events_in_bin

  n_event_cutoff <- 10

  weights_for_iso <- sqrt(n_in_bin)
  weights_for_iso[n_events_in_bin < n_event_cutoff] <- 0.01
  weights_for_iso <- weights_for_iso / sum(weights_for_iso)

  results$estimate_monotone <- -isotone::gpava(results$threshold,
                                               -results$estimate,
                                               weights = weights_for_iso)$x

  results$ci_lo_monotone <- results$estimate_monotone - (1.96 * results$se)
  results$ci_hi_monotone <- results$estimate_monotone + (1.96 * results$se)

  results

}


tmleThreshold.auto <- function(threshold, W, A, Y, Delta = rep(1, length(A)),
                               weights = rep(1, length(A)),
                               sl_library = c("SL.mean", "SL.glm"),
                               method = "method.CC_nloglik",
                               cvControl = list(V = 5),
                               run_mono = FALSE) {
  if (is.null(sl_library)) {
    sl_library <- c("SL.mean", "SL.glm", "SL.gam", "SL.ranger", "SL.earth")
  }

  n <- length(A)
  W <- as.matrix(W)
  colnames(W) <- paste0("W", 1:ncol(W))

  Av <- as.numeric(A >= threshold)
  positive_weight <- weights > 0
  q_index <- Delta == 1 & Av == 1 & positive_weight

  lrnr_Qv <- SuperLearner::SuperLearner(
    Y = Y[q_index],
    X = as.data.frame(W[q_index, , drop = FALSE]),
    newX = as.data.frame(W),
    family = binomial(),
    SL.library = sl_library,
    method = method,
    cvControl = cvControl,
    obsWeights = weights[q_index]
  )

  if (all(Av[positive_weight] == 1)) {
    lrnr_gv <- SuperLearner::SuperLearner(Y = Av, X = as.data.frame(W), family = binomial(), SL.library = "SL.mean",
                            method = method, cvControl = cvControl,
                            obsWeights = weights)
  } else {
    lrnr_gv <- SuperLearner::SuperLearner(Y = Av, X = as.data.frame(W), family = binomial(), SL.library = sl_library,
                            method = method, cvControl = cvControl,
                            obsWeights = weights)
  }

  if (all(Delta[positive_weight] == 1)) {
    lrnr_Gv <- SuperLearner::SuperLearner(Y = Delta, X = as.data.frame(cbind(W, Av)), family = binomial(), SL.library = "SL.mean",
                            method = method, cvControl = cvControl,
                            obsWeights = weights)
  } else {
    lrnr_Gv <- SuperLearner::SuperLearner(Y = Delta, X = as.data.frame(cbind(W, Av)), family = binomial(), SL.library = sl_library,
                            method = method, cvControl = cvControl,
                            obsWeights = weights)
  }

  Qv <- lrnr_Qv$SL.predict
  gv <- lrnr_gv$SL.predict
  Gv <- lrnr_Gv$SL.predict

  out_v <- suppressWarnings(tmle.inefficient(threshold, A, Delta, Y,
                                             gv, Qv, Gv,
                                             weights))
  psi <- out_v$psi
  IF <- out_v$EIF
  se <- sd(IF) / sqrt(n)
  CI_left <- psi - 1.96 * se
  CI_right <- psi + 1.96 * se

  list(threshold = threshold,
       estimate = psi,
       se = se,
       ci_lo = CI_left,
       ci_hi = CI_right)
}


#The binary-treatment-based TMLE (binTMLE):
tmle.inefficient <- function(threshold, A, Delta, Y,
                             gv, Qv, Gv,
                             weights = rep(1, length(A))) {
  n <- length(A)
  Av <- as.numeric(A >= threshold)
  gv <- truncate_pscore_adaptive(Av, gv)
  Gv <- truncate_pscore_adaptive(Delta, Gv)
  Qv <- pmax(Qv, 1e-8)
  Qv <- pmin(Qv, 1-1e-8)
  eps <- coef(glm.fit(as.matrix(rep(1,n)), Y,
                      offset = qlogis(Qv), weights = weights * Delta*Av/(gv*Gv),
                      family = binomial(), start = 0))
  Qv_star <- plogis(qlogis(Qv) + eps )
  psi <- weighted.mean(Qv_star, weights)
  EIF <- Delta*Av/(gv*Gv) * (Y - Qv_star) + Qv_star - psi
  EIF <- weights * EIF
  return(list(psi = psi, EIF = EIF))
}


truncate_pscore_adaptive <- function(A, pi, min_trunc_level = 1e-8) {
  risk_function <- function(cutoff, level) {
    pi <- pmax(pi, cutoff)

    alpha <- A/pi #Riesz-representor
    alpha1 <- 1/pi
    mean(alpha^2 - 2*(alpha1))
  }
  cutoff <- optim(1e-3, fn = risk_function, method = "Brent", lower = min_trunc_level, upper = 0.5, level = 1)$par

  pi <- pmax(pi, cutoff)
  return(pi)
}
