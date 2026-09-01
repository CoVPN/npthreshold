
#' Causal Survival Analysis using TMLE for discrete time-to-event variables
#'
#' This function performs causal survival analysis using targeted maximum likelihood estimation (TMLE) approach.
#' It estimates the treatment effect on survival outcomes while accounting for confounding variables and potential biases.
#' This method assumes the time-to-event variables are discrete with not too many levels.
#'
#' @param failure_time A numeric vector representing the time-to-event outcomes (e.g., time until death, failure, or censoring) for each observation.
#' The computational time and complexity of `survtmle_discrete` can significantly increase with the number of unique failure times.
#' For continuous or multi-valued discrete times, it is advisable to discretize the time-to-event variables into a grid of 20-30 time points.
#' Specifically, the dimension of the hazard estimation tasks generally scales linearly with the length of unique event times.
#' @param event_type A numeric vector of integers (0, 1, 2, ...) indicating the type of event for each observation. 0 is assumed to be censoring, and other values represent different types of failure events.
#' @param treatment A binary (0 or 1) or discrete (integer valued) treatment indicator vector indicating the treatment assignment for each observation.
#' @param covariates A named matrix specifying the names of potential counfounding variables in the dataset that should be used as adjustment covariates in the causal analysis.
#' @param weights A numeric vector of weights for each observation. Default is equal weights for all observations.
#' @param id A vector of unique identifiers for each observation. Default is a sequence from 1 to the length of 'treatment'.
#' @param learner.failure_time A required \code{binomial} \code{\link[sl3]{Lrnr_base}} learner object used for fitting the conditional hazard model for failure.
#' This learner is used as \code{binomial_learner} for \code{\link[sl3]{Lrnr_pooled_hazards}}.
#' @param learner.censoring_time A required \code{binomial} \code{\link[sl3]{Lrnr_base}} learner object used for fitting conditional hazard model for censoring.
#' This learner is used as \code{binomial_learner} for \code{\link[sl3]{Lrnr_pooled_hazards}}.
#' @param learner.event_type A required \code{binomial} \code{\link[sl3]{Lrnr_base}} learner object used for fitting the conditional probability distribution for failure time of the failure event type.
#' For nonbinary categorical event type, this is used as \code{binomial_learner} in \code{\link[sl3]{Lrnr_independent_binomial}}.
#' @param learner.treatment A required \code{binomial} \code{\link[sl3]{Lrnr_base}} learner object used for fitting the propensity score model for the treatment mechanism.
#' For nonbinary categorical treatment, this is used as \code{binomial_learner} in \code{\link[sl3]{Lrnr_independent_binomial}}.
#' @param cross_fit A logical value indicating whether to use cross-fitting (TRUE) or not (FALSE) when fitting the TMLE models.
#' See \code{\link[causalutils]{make_cross_fitted}}.
#' @param cross_validate A logical value indicating whether to wrap the supplied learners in cross-validation selectors.
#' @param failure_time.stratify_by_time A logical value indicating whether to stratify the analysis by time when fitting the conditional hazard model for failure time.
#' @param censoring_time.stratify_by_time A logical value indicating whether to stratify the analysis by time when fitting the conditional hazard model for censoring time.
#' @param calibrate A logical value indicating whether to perform calibration of treatment effect estimates using cross-validation.
#' @param autoML A logical value indicating whether learner arguments explicitly supplied as \code{NULL} should be replaced by \code{causalutils::get_autoML()}. The default is FALSE.
#' @param nfolds The number of folds to be used in cross-fitting. Default is 10.
#' @param max_iter The maximum number of iterations for the TMLE algorithm. Default is 200.
#' @param tol The convergence tolerance for the TMLE algorithm. Default is calculated based on the number of observations.
#' @param target_failure_time The target time points for which the treatment effect estimates will be computed. Default is the median of 'failure_time'.
#' @param target_treatment The target treatment values for which the treatment effect estimates will be computed. Default includes all unique values of 'treatment'.
#' @param target_event_type The target event types for which the treatment effect estimates will be computed. Default includes all unique non-zero values of 'event_type'.
#' @param verbose A logical value indicating whether to display progress and diagnostic messages during the computation. Default is TRUE.
#'
#' @return A list containing the estimated treatment effects and other relevant information from the TMLE analysis.
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(1234)
#' n <- 100
#' treatment <- rbinom(n, 1, 0.5)
#' covariates <- data.frame(W1 = rbinom(n, 1, 0.5), W2 = runif(n))
#' time <- sample(1:5, n, replace = TRUE)
#' status <- rbinom(n, 1, 0.5)
#'
#' result <- survtmle3_discrete(
#'   failure_time = time,
#'   event_type = status,
#'   treatment = treatment,
#'   covariates = covariates,
#'   target_failure_time = 4,
#'   target_treatment = 1,
#'   learner.treatment = sl3::Lrnr_glm$new(),
#'   learner.failure_time = sl3::Lrnr_glm$new(),
#'   learner.censoring_time = sl3::Lrnr_glm$new(),
#'   learner.event_type = sl3::Lrnr_glm$new(),
#'   cross_fit = FALSE
#' )
#' }
#' @seealso
#' \code{\link{survtmle}}: The original survtmle package in R.
#'
#' @keywords survival analysis, causal inference, TMLE, treatment effect
#' @import sl3 data.table origami
#' @export
survtmle3_discrete <- function(failure_time, event_type, treatment, covariates,
                               weights = rep(1, length(treatment)),
                               id = 1:length(treatment),
                               learner.treatment,
                               learner.failure_time,
                               learner.censoring_time,
                               learner.event_type,
                               target_failure_time = median(failure_time),
                               target_treatment = sort(unique(treatment)),
                               target_event_type = setdiff(unique(event_type), 0),
                               failure_time.stratify_by_time = FALSE,
                               censoring_time.stratify_by_time = FALSE,
                               cross_fit = TRUE,
                               cross_validate = FALSE,
                               calibrate = TRUE,
                               autoML = FALSE,
                               nfolds = 10,
                               max_iter = 200,
                               tol = 0.1/sqrt(sum(weights))/log(length(treatment)),
                               verbose = TRUE
) {
    treatment <- as.numeric(as.vector(treatment))
    covariates <- as.matrix(covariates)
    failure_time <- as.numeric(as.vector(failure_time))
    event_type <- as.numeric(as.vector(event_type))
    fold_type <- ifelse(cross_fit , "validation", "full")

  if(autoML) {
    default.learner <- causalutils::get_autoML()
    cross_fit <- TRUE
  } else {
    default.learner <- NULL
  }

  if(is.null(learner.treatment)) learner.treatment <- default.learner
  if(is.null(learner.failure_time)) learner.failure_time <- default.learner
  if(is.null(learner.censoring_time)) learner.censoring_time <- default.learner
  if(is.null(learner.event_type)) learner.event_type <- default.learner




  n <- length(treatment)
  weights <- n*weights/sum(weights)
  target_event_type <- sort(target_event_type)
  target_treatment <- sort(target_treatment)
  target_failure_time <- sort(target_failure_time)
  # estimate total failure hazard
  # estimate distribution of event_type
  # estimate censoring hazard


  # make origami fold object
  if(is.numeric(nfolds)) {
    folds <- origami::folds_vfold(n, nfolds)
  }
  # convert to categorical learner
  learner.treatment <- Lrnr_independent_binomial$new(binomial_learner = learner.treatment)
  learner.event_type <- Lrnr_independent_binomial$new(binomial_learner = learner.event_type)

  if(cross_validate) {
    # in case, learner is a stack, use cross validation selector too.
    learner.failure_time <- Lrnr_cv$new(learner.failure_time)
    learner.censoring_time <- Lrnr_cv$new(learner.censoring_time)
    # add multinomial loss
    learner.event_type <- make_learner(Pipeline, Lrnr_cv$new(learner.event_type), Lrnr_cv_selector$new(loss_loglik_multinomial))
    learner.treatment <- make_learner(Pipeline, Lrnr_cv$new(learner.treatment), Lrnr_cv_selector$new(loss_loglik_multinomial))
  }
  # NOTE: Work around for bug when causalutils::Lrnr_stratified has binomial_learner being cross-validated
  if(failure_time.stratify_by_time) learner.failure_time <- causalutils::Lrnr_stratified$new(learner.failure_time, variable_stratify = "t")
  if(censoring_time.stratify_by_time) learner.censoring_time <- causalutils::Lrnr_stratified$new(learner.censoring_time, variable_stratify = "t")
  # Add cv-selectors
  if(cross_fit) learner.failure_time <- make_learner(Pipeline, learner.failure_time, Lrnr_cv_selector$new(loss_squared_error))
  if(cross_fit) learner.censoring_time <- make_learner(Pipeline, learner.censoring_time, Lrnr_cv_selector$new(loss_squared_error))

  # process covariate matrix
  if(!is.matrix(covariates)) covariates <- as.matrix(covariates)
  if(length(colnames(covariates)) == 0 || is.null(colnames(covariates))) {
    colnames(covariates) <- paste0("W", 1:ncol(covariates))
  }
  covariate_names <- colnames(covariates)
  treatment_levels <- sort(unique(treatment))
  event_type_levels <- sort(setdiff(unique(event_type), 0))
  time_grid <- 1:max(failure_time)



  # process survival data and convert to long format
  data <- data.table(failure_time, event_type, treatment, covariates, id, weights)
  data_long <- transform_data_to_long(data, node_list = list(failure_time = "failure_time", event_type = "event_type", id = "id"))
  # store failure types
  event_types <- setdiff(sort(unique(data$event_type)), 0) # censoring is treated as a failure type for computation
  # get pooled folds for cv
  pooled_folds <- origami::make_folds(
    cluster_ids = data_long$id,
    V = nfolds
  )
  ##################################################
  ### estimate nuisance functions using learners ###
  ##################################################

  #####
  # estimate propensity score:
  ####
  # TODO: nonbinary discrete treatments
  if(verbose) print("Estimating propensity score for treatment mechanism... ")
  data_treatment <- data
  levels_trt <- sort(unique(data$treatment))
  if(length(levels_trt) == 1){
    treatment.hat <- as.matrix(rep(1, nrow(data_treatment)))
  } else {
  data_treatment$treatment <- factor(data_treatment$treatment, levels = levels_trt)
  #treatment_task <-  sl3_Task$new(data_treatment, covariates = covariate_names, outcome = "treatment", id = "id", weights = "weights")

  treatment_task <-  (sl3_Task$new(data_treatment, covariates = covariate_names, outcome = "treatment", id = "id", weights = "weights", folds = folds, outcome_type = variable_type(type = "categorical", levels = levels_trt)))
  # Even if binary, Lrnr_independent_binomial predicts a matrix.
  #learner.treatment_cat <-  Lrnr_independent_binomial$new(binomial_learner = learner.treatment)
  learner.treatment_trained <- learner.treatment$train(treatment_task)

  #print(as.data.table(learner.treatment_cat_trained$predict_fold(treatment_task)))
  treatment.hat <- sl3::unpack_predictions(learner.treatment_trained$predict_fold(treatment_task, fold_type))
}




  ######
  # estimate total failure hazard
  ######
  if(verbose) print("Estimating conditional hazard function for failure mechanism...")
  # get indices for rows corresponding to individuals still at risk at time t.
  in_risk_set <- which(data_long$in_risk_set == 1)
  # make time-pooled task (use user-supplied weights and id)
  pooled_hazard_task <- suppressWarnings(sl3_Task$new(data_long, covariates = c(covariate_names, "treatment", "t"), outcome = "dN", weights = "weights", id = "id", folds = pooled_folds, time = "t"))
  # train on those at risk of censoring or failure

  learner.failure_trained <- learner.failure_time$train(pooled_hazard_task[in_risk_set])
  # predict on all observations
  total.hazard.hats <- do.call(cbind, lapply(treatment_levels, function(level) {
    tmp <- data_long
    tmp$treatment <- level
    pooled_hazard_task <- suppressWarnings(sl3_Task$new(tmp, covariates = c(covariate_names, "treatment", "t"), outcome = "dN", weights = "weights", id = "id", folds = pooled_folds, time = "t"))
    #print("Here")
    #total.hazard.hat <-  learner.failure_trained$predict(pooled_hazard_task)
    #print("Here1")
    total.hazard.hat <-  learner.failure_trained$predict_fold(pooled_hazard_task, fold_type)

    return(total.hazard.hat)
  }))
  colnames(total.hazard.hats) <- as.character(treatment_levels)
  # store in list indexed by status

  ######
  # estimate censoring hazard
  ######
  if(verbose) print("Estimating conditional hazard function for censoring mechanism...")
  if(any(event_type == 0)) {
    # get indices for rows corresponding to individuals still at risk at time t.
    in_risk_set <- which(data_long$in_risk_set == 1)
    # make time-pooled task (use user-supplied weights and id)
    pooled_hazard_task <- suppressWarnings(sl3_Task$new(data_long, covariates = c(covariate_names, "treatment", "t"), outcome = "dC", weights = "weights", id = "id", folds = pooled_folds))
    # train on those at risk of censoring or failure
    learner.censored_trained <- learner.censoring_time$train(pooled_hazard_task[in_risk_set])
    # predict on all observations
    censor.hazard.hats <- do.call(cbind, lapply(treatment_levels, function(level) {
      tmp <- data_long
      tmp$treatment <- level
      pooled_hazard_task <- suppressWarnings(sl3_Task$new(tmp, covariates = c(covariate_names, "treatment", "t"), outcome = "dC", weights = "weights", id = "id", folds = pooled_folds))
      censor.hazard.hat <- as.vector(unlist(learner.censored_trained$predict_fold(pooled_hazard_task, fold_type)))
      return(censor.hazard.hat)
    }))
    colnames(censor.hazard.hats) <- as.character(treatment_levels)
  } else {
    censor.hazard.hats <- do.call(cbind, lapply(treatment_levels, function(level) {
      return(rep(0, nrow(data_long)))
    }))
  }


  ######
  # estimate event_type distribution
  ######
  if(verbose) print("Estimating conditional distribution for failure type mechanism...")
  if(length(event_types) > 1) {

    outcome_type <- variable_type(type = "categorical", levels = event_type_levels)
    not_censored <- which(data$event_type != 0)
    # set up learning tasks
    data_event_type <- data[, c(covariate_names, "treatment", "failure_time", "event_type", "weights", "id"), with = FALSE]
    sub_folds <- subset_folds(folds, not_censored)
    task_event_type_train <- suppressWarnings(sl3_Task$new(data_event_type[not_censored], covariates = c(covariate_names, "treatment", "failure_time"), outcome = "event_type", outcome_type = outcome_type, weights = "weights", id = "id", folds = sub_folds))
    # train and predict learner

    #learner.event_type.categorical <- Lrnr_independent_binomial$new(binomial_learner = learner.event_type)
    learner.event_type_trained <- learner.event_type$train(task_event_type_train)
    if (verbose) print( learner.event_type_trained$fit_object$formula)
    # matrix of event_type-specific estimates with columns following order of event_type_levels
    # at each time, get the probabilities.
    event_type.distr.hats <- lapply(treatment_levels, function(level) {
      event_type.distr.hats <- do.call(rbind, lapply(time_grid, function(t) {
        tmp <- data_event_type
        tmp$treatment <- level
        tmp$failure_time <- t
        task_event_type_predict <- suppressWarnings(sl3_Task$new(tmp, covariates = c(covariate_names, "treatment", "failure_time"), outcome = "event_type", outcome_type = outcome_type, weights = "weights", id = "id", folds = folds))

        event_type.distr.hats <- learner.event_type_trained$predict_fold(task_event_type_predict, fold_type)
        #print(head(data.table(event_type.distr.hats)))
        return(event_type.distr.hats)
      }))
      #print("here")
      #print(head(data.table(event_type.distr.hats)))
      # predictions are packed
      event_type.distr.hats <-  do.call(cbind, lapply(seq_along(event_type_levels), function(j) {
        # extract the jth index of each row
        as.vector(unlist(apply(event_type.distr.hats, 1, FUN = function(row){
          as.vector(unlist(row, use.names = FALSE)[[j]])
        })))
      }))

      colnames(event_type.distr.hats) <- as.character(event_type_levels)
      return(event_type.distr.hats)
    })
    names(event_type.distr.hats) <- as.character(treatment_levels)
  } else {
    event_type.distr.hats <- lapply(treatment_levels, function(level) {
      as.matrix(rep(1, nrow(data_long)))
    })
    names(event_type.distr.hats) <- as.character(treatment_levels)
  }
  # reformat





  # create list of nuisance estimates and parameters for tmle



  ##################################################
  ### targeted maximum likelihood estimation ###
  ##################################################

  # calibrate or adaptively truncate propensity scores away from zero



  treatment.hat <- apply(treatment.hat, 2, pmax, 0.005)
  treatment.hat <- apply(treatment.hat, 2, pmin, 1 - 0.005)
  colnames(treatment.hat) <- as.character(treatment_levels)


  max_eps <- 0.1
  tmle.spec <- list(max_eps = max_eps,  tol = tol, data_long = data_long, treatment_levels = treatment_levels, event_type_levels= event_type_levels, treatment.hat = treatment.hat, censor.hazard.hats = censor.hazard.hats, total.hazard.hats = total.hazard.hats,  event_type.distr.hats = event_type.distr.hats, weights = weights)


  #### initial


  get_estimates <- function(tmle.spec, a0, j0, target_failure_time) {
    total.hazard.hats <- tmle.spec$total.hazard.hats
    event_type.distr.hat_a0_j0 <- tmle.spec$event_type.distr.hat_a0_j0

    S_failure_left <- long_hazard_to_survival_mats(total.hazard.hats, time_grid = time_grid, left_cont = TRUE)
    names(S_failure_left) <- as.character(treatment_levels)
    S_failure_left_a0 <- as.vector(S_failure_left[[as.character(a0)]])

    event_type.distr.hat_a0_j0 <- as.vector(event_type.distr.hats[[as.character(a0)]][, match(j0, event_type_levels)])
    total.hazard.hats_a0 <- total.hazard.hats[, match(a0, treatment_levels)]

    F_failure_a0_j0 <- S_failure_left_a0 * total.hazard.hats_a0 * event_type.distr.hat_a0_j0
    F_failure_a0_j0 <- matrix(F_failure_a0_j0, ncol = length(time_grid))

    F_failure_a0_j0 <- t(apply(F_failure_a0_j0, 1, cumsum))
    F_failure_a0_j0_t0 <- F_failure_a0_j0[, match(target_failure_time, time_grid), drop = FALSE]

    ests <- apply(as.matrix(F_failure_a0_j0_t0), 2, weighted.mean, weights)
    return(ests)
  }


  ####


  output <- as.data.table(expand.grid(target_treatment, target_event_type))
  names(output) <- c("treatment", "event_type")
  output$times <- list()
  output$estimates <- list()
  output$se <- list()
  output$EIF <- list()
  for(j0 in target_event_type) {
    for(a0 in target_treatment) {
      tmle.spec$target_event_type <- j0
      tmle.spec$target_failure_time <- target_failure_time
      tmle.spec$target_treatment <- a0
      if (verbose) print( "ESTIMATE")
      if (verbose) print( get_estimates(tmle.spec, a0, j0, target_failure_time))


      # target distribution of event_type
      if(length(event_types) > 1) {
        out_event_type <- suppressWarnings(sl3:::call_with_args(.target_event_type_distr, tmle.spec, silent = TRUE))
        tmle.spec$event_type.distr.hats <- out_event_type$event_type.distr.hats
        score_event_type <- out_event_type$score
        if (verbose) print( paste0("type score solved: ", score_event_type))
      } else {
        out_event_type <- list(EIF = 0)
      }
      if (verbose) print( "ESTIMATE after strain targeting")
      if (verbose) print( get_estimates(tmle.spec, a0, j0, target_failure_time))
      # target the failure hazard estimator
      converged_flag <- FALSE
      for(iter in 1:max_iter) {
        if (verbose) print( iter)
        out_hazard <- suppressWarnings(sl3:::call_with_args(.target_failure_hazard, tmle.spec, silent = TRUE))
        tmle.spec$total.hazard.hats <- out_hazard$total.hazard.hats
        if (verbose) print( "ESTIMATE after hazard targeting once")
        if (verbose) print( get_estimates(tmle.spec, a0, j0, target_failure_time))
        score_hazard <- out_hazard$score
        if (verbose) print( paste0("HAZARD SCORE: ", abs(score_hazard)))
        if (verbose) print( paste0("TOL: ", abs(tol)))
        if(abs(score_hazard) <= tol) {
          if (verbose) print( paste0("failiure score solved: ", score_hazard))
         # print(iter)
        #  print(abs(score_hazard))
          converged_flag <- TRUE
          break
        }
      }
      if (!converged_flag) {
        stop(
          sprintf(
            "TMLE targeting did not converge after %d iterations; score = %g, tolerance = %g.",
            max_iter,
            score_hazard,
            tol
          ),
          call. = FALSE
        )
      }
      total.hazard.hats <- tmle.spec$total.hazard.hats
      event_type.distr.hat_a0_j0 <- tmle.spec$event_type.distr.hat_a0_j0

      S_failure_left <- long_hazard_to_survival_mats(total.hazard.hats, time_grid = time_grid, left_cont = TRUE)
      names(S_failure_left) <- as.character(treatment_levels)
      S_failure_left_a0 <- as.vector(S_failure_left[[as.character(a0)]])

      event_type.distr.hat_a0_j0 <- as.vector(event_type.distr.hats[[as.character(a0)]][, match(j0, event_type_levels)])
      total.hazard.hats_a0 <- total.hazard.hats[, match(a0, treatment_levels)]

      F_failure_a0_j0 <- S_failure_left_a0 * total.hazard.hats_a0 * event_type.distr.hat_a0_j0
      F_failure_a0_j0 <- matrix(F_failure_a0_j0, ncol = length(time_grid))

      F_failure_a0_j0 <- t(apply(F_failure_a0_j0, 1, cumsum))
      F_failure_a0_j0_t0 <- F_failure_a0_j0[, match(target_failure_time, time_grid), drop = FALSE]

      ests <- apply(as.matrix(F_failure_a0_j0_t0), 2, weighted.mean, weights)

      # Compute covariate component of EIF
      EIF_W <- do.call(cbind,lapply(1:length(target_failure_time), function(index){
        comp <- F_failure_a0_j0_t0[, index] - ests[index]
        #rep(weights*comp, length(time_grid))
        return(weights*comp)
      }))
      # Hazard component of EIF
      EIF_hazard <- apply(out_hazard$EIF, 2, function(eif){
        rowSums(matrix(eif, ncol = length(time_grid) , byrow = FALSE))
      })
      # Event type component of EIF
      EIF_event_type <- apply(as.matrix(out_event_type$EIF), 2, function(eif){
        rowSums(matrix(eif, ncol = length(time_grid) , byrow = FALSE))
      })


      EIF <- EIF_hazard + EIF_event_type + EIF_W
      row_index <- output$event_type==j0 & output$treatment==a0

      #print(c(sd(EIF_hazard), sd(as.vector(EIF_event_type)), sd( EIF_W)))


      ses <- apply(EIF, 2, sd)/sqrt(sum(weights))
      #print(ses)
      #stop("hi")
      names(ests) <- paste0("time_", target_failure_time)
      names(ses) <- paste0("time_", target_failure_time)

      output[row_index, "EIF"] <- list(EIF)
      output[row_index, "estimates"] <- list(ests)
      output[row_index, "se"] <- list(ses)
      CI <- cbind(ests - qnorm(0.975) * ses, ests + qnorm(0.975) * ses)
      colnames(CI) <- c("CI_left", "CI_right")
      rownames(CI) <- paste0("time_", target_failure_time)
      output[row_index, "CI"] <- list(CI)
      output[row_index, "times"] <- list(target_failure_time)
    }
  }




  tmle.spec <- list(max_eps = max_eps,   data_long = data_long, treatment_levels = treatment_levels, event_type_levels= event_type_levels, treatment.hat = treatment.hat, censor.hazard.hats = censor.hazard.hats, total.hazard.hats = total.hazard.hats,  event_type.distr.hats = event_type.distr.hats, weights = weights)




  #CI <- data.table(cbind(output$estimate - qnorm(0.975) * output$se, output$estimate + qnorm(0.975) * output$se), conf = "95%")
  #names(CI) <- c("CI_left", "CI_right", "conf")

  return(output)

}
