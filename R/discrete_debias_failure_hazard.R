
# internal use
.target_failure_hazard <- function(max_eps = 0.05, target_failure_time, target_treatment, target_event_type, data_long, treatment_levels, event_type_levels, treatment.hat, censor.hazard.hats, total.hazard.hats, event_type.distr.hats, weights, tol,
                                   verbose = F) {
  if(any(is.na(unlist(treatment.hat)))){
    stop("NA values found in propensity score.")
  }
  if(any(is.na(unlist(censor.hazard.hats)))){
    stop("NA values found in censoring hazard.")
  }
  if(any(is.na(unlist(total.hazard.hats)))){
    stop("NA values found in censoring hazard.")
  }
  if(any(is.na(unlist(event_type.distr.hats)))){
    stop("NA values found in event type distribution.")
  }
  if(any(is.na(unlist(weights)))){
    stop("NA values found in weights.")
  }
  n <- sum(weights)

  time_grid <- 1:max(data_long$failure_time)

  if(!(target_failure_time %in% time_grid)) {
    stop(paste0("target_failure_time, t0", target_event_type, ", not found in time_grid = range(", range(time_grid),")", collapse = ", "))
  }
  if(!(target_treatment %in% treatment_levels)) {
    stop(paste0("target_treatment, a0=", target_treatment, ", not found in treatment_levels = c(", treatment_levels,")", collapse = ", "))
  }
  if(!(target_event_type %in% event_type_levels)) {
    stop(paste0("target_event_type, j0=", target_event_type, ", not found in event_type_levels = c(", event_type_levels,")", collapse = ", "))
  }
  # get propensity score
  treatment.a0.hat <- as.vector(treatment.hat[, match(target_treatment, treatment_levels)])
  # a list of matrix survival estimates indexed by treatment level
  # for censoring, P(C >= t| ... )
  S_censor <- long_hazard_to_survival_mats(censor.hazard.hats, time_grid = time_grid, left_cont = TRUE)
  names(S_censor) <- as.character(treatment_levels)
  # for total failure, P(T > t| ...)
  S_failure <- long_hazard_to_survival_mats(total.hazard.hats, time_grid = time_grid)
  names(S_failure) <- as.character(treatment_levels)
  S_failure_left <- long_hazard_to_survival_mats(total.hazard.hats, time_grid = time_grid, left_cont = TRUE)
  names(S_failure_left) <- as.character(treatment_levels)



  # survivals at target treatment
  S_failure_a0 <- as.vector(S_failure[[as.character(target_treatment)]])
  S_failure_left_a0 <- as.vector(S_failure_left[[as.character(target_treatment)]])
  S_censor_a0 <- as.vector(S_censor[[as.character(target_treatment)]])

  # event_type distribution at target event_type
  event_type.distr.hat_a0_j0 <- as.vector(event_type.distr.hats[[as.character(target_treatment)]][, match(target_event_type, event_type_levels)])


  # hazard at target treatment
  total.hazard.hats_a0 <- total.hazard.hats[, match(target_treatment, treatment_levels)]
  # cause-specific cumulative incidence
  F_failure_a0_j0 <- S_failure_left_a0 * total.hazard.hats_a0 * event_type.distr.hat_a0_j0

  F_failure_a0_j0 <- matrix(F_failure_a0_j0, ncol = length(time_grid))
  F_failure_a0_j0 <- t(apply(F_failure_a0_j0, 1, cumsum))


  F_failure_a0_j0_t0 <- F_failure_a0_j0[, match(target_failure_time, time_grid), drop = FALSE]
  F_failure_a0_j0 <- as.vector(F_failure_a0_j0)


  # indicator of intervention treatment
  treatment_ind <- 1*(data_long$treatment == target_treatment)


  if (verbose) print(paste0("est step", weighted.mean(F_failure_a0_j0_t0, weights)))
  #print(quantile(event_type.distr.hat_a0_j0))
  #stop("Wait")

  # indicator of not being censored prior to time t
  censor_ind <- data_long$in_risk_set

  # bound survival probability away from 0 for inverse weighting stability
  S_failure_a0 <- pmax(S_failure_a0, 25/sqrt(n)/log(n))
  S_censor_a0 <- pmax(S_censor_a0, 25/sqrt(n)/log(n))
  # make clever covariates
  treatment.a0.hat <- pmax(treatment.a0.hat, 1e-6)

  H1.hat <- as.matrix((treatment_ind/treatment.a0.hat) * (censor_ind/S_censor_a0))




  if(any(is.na(unlist(H1.hat)))) {
    print(data.table(H1.hat))
    stop("NA values found in clever weights H1.hat.")
  }
  # should be matrix with ncols number of tgt times
  H3.hat <- as.matrix(event_type.distr.hat_a0_j0 - apply(F_failure_a0_j0_t0, 2, function(v){v-F_failure_a0_j0}) / S_failure_a0)
  # zero out times after target time.
  H3.hat <- do.call(cbind, lapply(seq_along(target_failure_time), function(index){
    t0 <- target_failure_time[[index]]
    zero_out <- 1*(data_long$t <= t0)
    return(zero_out *H3.hat[,index])
  }))


  if(any(is.na(unlist(H3.hat)))) {
    print(paste0("target times: ", target_failure_time))
    print(data.table(S_failure_a0, F_failure_a0_j0, F_failure_a0_j0_t0, event_type.distr.hat_a0_j0, H3.hat))
    stop("NA values found in clever covariate H3.hat.")
  }


  dN <- data_long$dN
  in_risk_set <- data_long$in_risk_set




  IF <- as.matrix(in_risk_set * as.vector(H1.hat) * weights * H3.hat * as.vector(dN - total.hazard.hats_a0))
  IF_long <- IF
  # sum over time
  if(any(is.na(unlist(IF)))) {
    stop("NA values found in influence function.")
  }
  IF <- as.matrix(apply(IF, 2, function(v){rowSums(matrix(v, ncol = length(time_grid)))}))
  scores <- apply(IF, 2, mean)
  if (verbose) print("SCORES")
  if (verbose) print(scores)
  direction <- scores/sqrt(mean(scores^2))

  if(any(is.na(direction))) {
    print("Weighted mean etimate")
    try({print(weighted.mean(as.vector(F_failure_a0_j0), weights))})
    print(quantile(H3.hat))
    print(quantile(H1.hat))
    print(quantile(in_risk_set))
    print(quantile(as.vector(dN - total.hazard.hats_a0)))
    warning("Direction of least favorable submodel has NA values. Setting equal to zero.")
    direction[is.na(direction)] <- 0
  }



  #total.hazard.hats_a0 <- pmin(pmax(total.hazard.hats_a0, 1e-16), 1 - 1e-16)
  risk_function <- function(eps) {
    # logistic fluctuation submodel
    haz_eps <- plogis(qlogis(total.hazard.hats_a0) + eps * (H3.hat%*%direction))
    return(weighted.mean(-1 * ifelse(dN == 1, log(haz_eps), log(1- haz_eps)), in_risk_set  * H1.hat  * weights ))
    # weighted log-likelihood loss
    #mean(in_risk_set * H1.hat * weights * ifelse(dN == 1, -log(haz_eps), -log(1-haz_eps)))
  }
  scale <- 1 #sqrt(weighted.mean(H3.hat%*%direction^2,   in_risk_set  * H1.hat  * weights ))

  wts <- in_risk_set  * H1.hat  * weights



  #### Add up to t0
  max_eps <- max_eps/scale


  # only update if tolerance not reached.
  if(sqrt(mean(scores^2)) > tol) {
    eps.hat <- optim(0, risk_function, method = "Brent", lower = -max_eps, upper = max_eps)$par
  } else {
    eps.hat <- 0
  }

  total.hazard.hats_a0 <- plogis(qlogis(total.hazard.hats_a0) + eps.hat * as.vector(H3.hat%*%direction))
  total.hazard.hats[, match(target_treatment, treatment_levels)] <- total.hazard.hats_a0

  if(any(is.na(unlist(total.hazard.hats_a0)))) {
    stop("NAs found in targeted/fluctuated hazard.")
  }


  return(list(total.hazard.hats = total.hazard.hats, score = sqrt(mean(scores^2)), EIF = IF_long))
}
