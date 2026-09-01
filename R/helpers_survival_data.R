# internal use
long_hazard_to_survival_mats <- function(hazard.hats, time_grid, left_cont = FALSE, as_vector = TRUE) {
  apply(hazard.hats, 2, function(hazard.hat){
    # convert hazard to (n by t) matrix with columns corresponding to times.
    hazard.hat.mat <- matrix(hazard.hat, ncol = length(time_grid), byrow = FALSE)
    out <- t(apply(1 - hazard.hat.mat, 1, cumprod))
    if(left_cont) {
      # left continuous version
      # add columns of 1's and drop last column
      out <- cbind(1, out)
      out <- out[, - ncol(out) , drop = FALSE]
    }
    if(as_vector) as_vector <- as.vector(out)
    return(out)
  }, simplify = FALSE)
}

# internal use
transform_data_to_long = function(data, node_list) {
  as.data.table(data)
  failure_time_name <- node_list$failure_time
  event_type_name <- node_list$event_type
  failure_time_data <- data[[failure_time_name]]
  event_type_data <- data[[event_type_name]]


  if (is.null(node_list$id) & !("id" %in% names(data))) {
    id <- 1:nrow(data)
    data <- cbind(id = id, data)
    node_list$id <- "id"
  }

  t_grid <- 1:max(failure_time_data)
  type_grid <- sort(unique(event_type_data))
  all_times <- lapply(t_grid, function(t_current) {
      df_time <- copy(data)
      df_time$t <- t_current
      # set censoring hazard indicator (event_type == 0)
      df_time$dC <- as.numeric(t_current == failure_time_data & event_type_data == 0)
      # set hazard failure indicator (event_type == type_current)
      df_time$dN <- as.numeric(t_current == failure_time_data & event_type_data != 0)
      # ties not possible since only one type attributed per individual
      # set(df_time, df_time$dC == 1 & df_time$dNj == 1, "dC", 0)
      # set indicator that individual is still at risk (in_risk_set)
      df_time$in_risk_set <- as.numeric(t_current <= failure_time_data)
      return(df_time)
    })




  df_long <- rbindlist(all_times)

  return(df_long)
}
