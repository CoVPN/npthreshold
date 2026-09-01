#' Non-parametric estimation of the covariate-adjusted threshold-response function
#' adjusting for right-censoring
#'
#' This method estimates \eqn{E[E[Y|A \ge v, W]]} for a range of thresholds \eqn{v}, where \eqn{Y}
#' is a binary outcome of interest that is subject to right-censoring, \eqn{A} is a
#' continuous biomarker of interest, and \eqn{W} are baseline variables. To
#' account for missingness in the marker \eqn{A}, the user can
#' provide inverse probability observation weights in the function call.
#'
#' @param data A \code{data.frame} containing the dataset to be analyzed.
#' @param covariates A \code{character} vector specifying the names of columns to be used as adjustment covariates.
#' @param failure_time A \code{character} string indicating the name of the column representing the failure time variable.
#' @param event_type A \code{character} string specifying the name of the column representing the event type variable.
#' @param marker A \code{character} string for the name of the column representing the marker variable.
#' @param weights A \code{character} string specifying the name of a
#'   nonnegative observation-weight column. If \code{NULL}, all observations
#'   receive weight one. Observations with zero weight are excluded from
#'   descriptive counts and may have a missing marker value.
#' @param threshold_list A \code{numeric} vector of threshold values for
#'   estimation. If \code{NULL}, estimation uses the minimum marker value among
#'   observations with positive weight.
#' @param tf A \code{numeric} value indicating the reference time point for the analysis.
#' @param nbins_time A \code{numeric} value specifying the number of time bins to be used.
#' @param learner.censoring_time A required \code{binomial} \code{\link[sl3]{Lrnr_base}} learner object used for fitting conditional hazard model for censoring.
#' @param learner.event_type A required \code{binomial} \code{\link[sl3]{Lrnr_base}} learner object used for fitting the conditional probability distribution for failure time of the failure event type.
#' @param learner.treatment A required \code{binomial} \code{\link[sl3]{Lrnr_base}} learner object used for fitting the propensity score model for the treatment mechanism.
#' @param learner.failure_time A required \code{binomial} \code{\link[sl3]{Lrnr_base}} learner object used for fitting the conditional hazard model for failure.
#' @param verbose A logical value indicating whether to display progress and diagnostic messages during the computation. Default is FALSE.
#' @param placebo_risk Placebo risk (to calculate VE columns). If NULL (default), will not
#' include VE results.
#' @param placebo_se Placebo risk standard error (to calculate VE standard errors)
#' If NULL (default), will not include VE results.
#' @param modify_left_point If \code{TRUE}, replace the estimate and standard
#'   error at the smallest threshold with \code{left_estimate} and
#'   \code{left_se}.
#' @param left_estimate Replacement estimate used when
#'   \code{modify_left_point = TRUE}.
#' @param left_se Replacement standard error used when
#'   \code{modify_left_point = TRUE}.
#' @param cross_fit Whether to use cross-fitted nuisance predictions. The
#'   default is \code{FALSE} to preserve existing behavior; \code{TRUE} is
#'   recommended for flexible machine-learning learners.
#' @param nfolds Number of folds used when \code{cross_fit = TRUE}.
#'
#' @return A data.frame containing the following columns for each threshold in `threshold_list`:
#' - `estimate`: The estimated value.
#' - `se`: The standard error of the estimate.
#' - `ci_lo`: The lower bound of the pointwise 95% confidence interval.
#' - `ci_hi`: The upper bound of the pointwise 95% confidence interval.
#' - `estimate_monotone`: The monotone-adjusted estimate.
#' - `ci_lo_monotone`: The lower bound of the monotone-adjusted confidence interval.
#' - `ci_hi_monotone`: The upper bound of the monotone-adjusted confidence interval.
#' - `n_in_bin`: The number of positive-weight observations at or above the threshold.
#' - `n_events_in_bin`: The number of type-1 events by `tf` among positive-weight
#'   observations at or above the threshold.
#' - `person_time_in_bin`: The total observed person-time through `tf` among
#'   positive-weight observations at or above the threshold.
#'
#' @import sl3 data.table origami
#' @examples
#' @export
thresholdSurv <- function(data,
                          covariates,
                          failure_time,
                          event_type,
                          marker,
                          tf,
                          weights = NULL,
                          threshold_list = NULL,
                          nbins_time = 20,
                          verbose = FALSE,
                          learner.treatment,
                          learner.event_type,
                          learner.failure_time,
                          learner.censoring_time,
                          placebo_risk = NULL,
                          placebo_se = NULL,
                          modify_left_point = FALSE,
                          left_estimate = NA,
                          left_se = NA,
                          cross_fit = FALSE,
                          nfolds = 10
) {

  data <- data.table::copy(data.table::as.data.table(data))
  tf <- as.numeric(tf)

  if (length(tf) != 1L || !is.finite(tf)) {
    stop("`tf` must be one finite number.", call. = FALSE)
  }

  if (xor(is.null(placebo_risk), is.null(placebo_se))) {
    stop(
      "`placebo_risk` and `placebo_se` must either both be supplied or both be NULL.",
      call. = FALSE
    )
  }

  if (!is.null(placebo_risk)) {
    if (length(placebo_risk) != 1L ||
        !is.finite(placebo_risk) ||
        placebo_risk <= 0) {
      stop("`placebo_risk` must be one positive finite number.", call. = FALSE)
    }

    if (length(placebo_se) != 1L ||
        !is.finite(placebo_se) ||
        placebo_se < 0) {
      stop("`placebo_se` must be one non-negative finite number.", call. = FALSE)
    }
  }

  weight_values <- if (is.null(weights)) {
    rep(1, nrow(data))
  } else {
    data[[weights]]
  }

  if (!is.numeric(weight_values) ||
      any(!is.finite(weight_values)) ||
      any(weight_values < 0) ||
      sum(weight_values) <= 0) {
    stop("`weights` must be finite, nonnegative, and have a positive sum.",
         call. = FALSE)
  }

  # subset to relevant variables
  data_select <- data[, unique(c(covariates, failure_time,
                                 event_type, marker)), with = FALSE]

  marker_values <- data_select[[marker]]
  positive_weight <- weight_values > 0

  if (!is.numeric(marker_values) ||
      any(!is.finite(marker_values[positive_weight]))) {
    stop("`marker` must be finite for observations with positive weight.",
         call. = FALSE)
  }
  marker_values[!positive_weight] <- 0
  data.table::set(data_select, , marker, marker_values)

  # Retain follow-up in its original units for descriptive summaries. The
  # estimator below uses a discretized time index instead.
  failure_time_original <- data_select[[failure_time]]

  if (!is.numeric(failure_time_original) ||
      any(!is.finite(failure_time_original))) {
    stop("`failure_time` must contain only finite numeric values.",
         call. = FALSE)
  }
  if (tf < min(failure_time_original) || tf > max(failure_time_original)) {
    stop("`tf` must fall within the observed failure-time range.",
         call. = FALSE)
  }

  if (!is.logical(modify_left_point) || length(modify_left_point) != 1L ||
      is.na(modify_left_point)) {
    stop("`modify_left_point` must be TRUE or FALSE.", call. = FALSE)
  }
  if (modify_left_point &&
      (length(left_estimate) != 1L || !is.finite(left_estimate) ||
       length(left_se) != 1L || !is.finite(left_se) || left_se < 0)) {
    stop(
      "`left_estimate` must be finite and `left_se` must be finite and nonnegative.",
      call. = FALSE
    )
  }

  # discretize
  time_grid <- unique(quantile(data_select[[failure_time]],
                               seq(0, 1,length = nbins_time+1),
                               type = 1))

  # add time of interest to grid
  time_grid <- sort(union(time_grid, tf))
  failure_time_discrete <- findInterval(data_select[[failure_time]],
                                        time_grid,
                                        left.open = TRUE,
                                        all.inside = TRUE)
  tf_discrete <- findInterval(tf,
                              time_grid,
                              left.open = TRUE,
                              all.inside = TRUE)
  after_tf_in_target_bin <- failure_time_original > tf &
    failure_time_discrete <= tf_discrete
  failure_time_discrete[after_tf_in_target_bin] <- tf_discrete + 1L
  data_select[[failure_time]] <- failure_time_discrete

  # if threshold list is null, then estimate on whole dataset
  if (is.null(threshold_list)){
    threshold_list <- min(marker_values[positive_weight])
  }

  out_list <- list()

  covariates_matrix <- model.matrix(~ . - 1, data = data_select[, covariates, with = FALSE])


  for(threshold in threshold_list ) {
    if (verbose) print(paste0("THRESHOLD: ", threshold))

    data_select[["trt_temp"]] <- as.numeric(data_select[[marker]] >= threshold)
    if(sum(data_select[["trt_temp"]][positive_weight]) == 0){
      stop(sprintf("There were zero observations above the threshold %s.",
                   threshold),
           call. = FALSE)
    }

    survout <- tryCatch(
      survtmle3_discrete(data_select[[failure_time]],
                         data_select[[event_type]],
                         data_select[["trt_temp"]],
                         covariates_matrix,
                         weights = weight_values,
                         learner.treatment = learner.treatment,
                         learner.failure_time = learner.failure_time,
                         learner.censoring_time = learner.censoring_time,
                         learner.event_type = learner.event_type,
                         target_failure_time = tf_discrete,
                         target_treatment = c(1),
                         target_event_type = 1,
                         failure_time.stratify_by_time = FALSE,
                         censoring_time.stratify_by_time = FALSE,
                         cross_fit = cross_fit,
                         cross_validate = FALSE,
                         calibrate = FALSE,
                         nfolds = nfolds,
                         verbose = verbose,
                         max_iter = 100,
                         tol = 1e-7),
      error = function(e) {
        stop(sprintf("Estimation failed at threshold %s: %s",
                     threshold,
                     conditionMessage(e)),
             call. = FALSE)
      }
    )

    survout$threshold <- threshold
    out_list[[as.character(threshold)]] <- survout
  }

  res <- rbindlist(out_list)

  # make results table
  res$estimate <- unlist(res$estimates )
  res$se <- unlist(res$se )
  res$ci_lo <- unlist(res$CI)[c(TRUE, FALSE)]
  res$ci_hi <- unlist(res$CI)[c(FALSE, TRUE)]
  res$times <- tf

  n_in_bin <- sapply(threshold_list, function(thresh) {
    sum(positive_weight & marker_values >= thresh)
  })

  n_events_in_bin <- sapply(threshold_list, function(thresh) {
    in_bin <- positive_weight & marker_values >= thresh
    sum(in_bin &
          data_select[[event_type]] == 1 &
          failure_time_original <= tf)
  })

  person_time_in_bin <- sapply(threshold_list, function(thresh) {
    in_bin <- positive_weight & marker_values >= thresh
    sum(pmin(failure_time_original[in_bin], tf))
  })

  if (modify_left_point){
    left_index <- which.min(res$threshold)
    res[left_index, "estimate"] <- left_estimate
    res[left_index, "se"] <- left_se
    res[left_index, "ci_lo"] <- left_estimate - (1.96 * left_se)
    res[left_index, "ci_hi"] <- left_estimate + (1.96 * left_se)
  }

  RCDF <- function(a) {
    sum(weight_values * (marker_values >= a)) / sum(weight_values)
  }

  RCDF <- Vectorize(RCDF)

  res$rcdf <- RCDF(res$threshold)

  res$n_in_bin <- n_in_bin
  res$n_events_in_bin <- n_events_in_bin
  res$person_time_in_bin <- person_time_in_bin

  n_event_cutoff <- 10

  weights_for_iso <- sqrt(n_in_bin)
  weights_for_iso[n_events_in_bin < n_event_cutoff] <- 0.01
  weights_for_iso <- weights_for_iso / sum(weights_for_iso)

  res$estimate_monotone <- -isotone::gpava(res$threshold,
                                               -res$estimate,
                                               weights = weights_for_iso)$x

  res$ci_lo_monotone <- res$estimate_monotone - (1.96 * res$se)
  res$ci_hi_monotone <- res$estimate_monotone + (1.96 * res$se)

  res <- res[,c("threshold", "rcdf", "n_in_bin", "n_events_in_bin",
                "person_time_in_bin", "estimate", "se",
                "ci_lo", "ci_hi", "estimate_monotone", "ci_lo_monotone", "ci_hi_monotone")]

  if (!is.null(placebo_risk) & !is.null(placebo_se)){
    res <- res %>%
      dplyr::mutate(ve_monotone = 1 - (estimate_monotone / placebo_risk))

    res$ve_se <- apply(res, 1, calc_delta_method,
                       placebo_risk = placebo_risk,
                       placebo_se = placebo_se)

    res <- res %>%
      dplyr::mutate(ve_monotone_ci_lo = ve_monotone - (1.96 * ve_se),
                    ve_monotone_ci_hi = ve_monotone + (1.96 * ve_se))
  }

  res

}

calc_delta_method <- function(row, placebo_risk, placebo_se) {
  # Extract values from the row
  estimate_monotone <- row['estimate_monotone']
  estimate_se  <- row['se']

  # Variance-covariance matrix (assuming no covariance)
  cov_matrix <- matrix(c(estimate_se^2, 0, 0, placebo_se^2), nrow = 2)

  # Calculate the standard error using the delta method
  se_ratio <- msm::deltamethod(~ x1 / x2,
                               mean = c(estimate_monotone, placebo_risk),
                               cov = cov_matrix)

  return(se_ratio)
}
