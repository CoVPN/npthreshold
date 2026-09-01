#' Plot threshold-response estimates
#'
#' This function creates a plot of estimates and confidence intervals for
#' \eqn{E[E[Y|A \ge v, W]]} across various thresholds.
#'
#' @param res A \code{data.frame} obtained from \code{thresholdSurv()} or
#'   \code{thresholdBinary()}.
#' @param type Either `"raw"`, `"monotone"`, or `"ve"`.
#' @param ylabel Label for the y-axis.
#' @param xlabel Label for the x-axis.
#' @param cutoffs A numeric vector \code{c(lower, upper)} used to clip the lower
#'   and upper confidence bounds, respectively. The point estimates are not
#'   clipped.
#' @param exp10 If \code{TRUE}, label x-axis values as powers of ten.
#' @param plot_density If \code{TRUE}, add the weighted reverse cumulative
#'   distribution function (reverse CDF) of the marker.
#' @param plot_endpoints If \code{TRUE}, add points for individual disease
#'   endpoints.
#' @param data Original dataset. Required when \code{plot_density = TRUE},
#'   \code{plot_endpoints = TRUE}, or an annotation is supplied, and when
#'   \code{exp10 = TRUE} without \code{xlimit}.
#' @param weights Name of the weight variable in \code{data}. Required when
#'   \code{plot_density = TRUE}.
#' @param marker Name of the marker variable in \code{data}. Required for
#'   density or endpoint plotting and annotations, and when \code{exp10 = TRUE}
#'   without \code{xlimit}.
#' @param annotate String for text to annotate in the corner of the graph
#' @param event Name of the event variable in \code{data}. Required when
#'   \code{plot_endpoints = TRUE}.
#' @param time_var Name of the time-to-event variable in \code{data}. Required
#'   when \code{plot_endpoints = TRUE}.
#' @param tf Reference time point used to select endpoints occurring on or
#'   before that time. Required when \code{plot_endpoints = TRUE}.
#' @param ylimit Numeric vector \code{c(min, max)} for y-axis limits, optionally
#'   followed by a positive tick step: \code{c(min, max, tick_step)}.
#' @param xlimit Numeric vector \code{c(min, max)} for x-axis limits, optionally
#'   followed by a positive tick step: \code{c(min, max, tick_step)}. When the
#'   step is omitted, tick locations are chosen automatically.
#' @param yright_endpoint Either `"min"` or `"max"`; controls whether endpoint
#'   interpolation beyond the observed threshold range uses the minimum or
#'   maximum estimate.
#' @return A \code{ggplot2} object displaying the estimates and confidence
#' intervals for \eqn{E[E[Y|A \ge v, W]]} across the specified thresholds.
#' @examples
#' estimates <- data.frame(
#'   threshold = 1:3,
#'   estimate = c(0.20, 0.25, 0.31),
#'   ci_lo = c(0.15, 0.19, 0.24),
#'   ci_hi = c(0.25, 0.31, 0.38)
#' )
#' graphthresh(estimates, xlimit = c(1, 3))
#' graphthresh(estimates, xlimit = c(1, 3, 0.5))
#' @import ggplot2
#' @export
graphthresh <- function(res,
                        type = "raw",
                        ylabel = "Estimate",
                        xlabel = "Thresholds",
                        cutoffs = NULL,
                        exp10 = FALSE,
                        plot_density = FALSE,
                        plot_endpoints = FALSE,
                        data = NA,
                        weights = NA,
                        marker = NA,
                        annotate = NA,
                        event = NA,
                        time_var = NA,
                        tf = NA,
                        ylimit = NULL,
                        xlimit = NULL,
                        yright_endpoint = "min"){

  res <- data.frame(res)

  parse_axis_limit <- function(axis_limit, argument) {
    if (is.null(axis_limit)) {
      return(list(limits = NULL, breaks = NULL))
    }

    if (!is.numeric(axis_limit) ||
        !length(axis_limit) %in% c(2L, 3L) ||
        any(!is.finite(axis_limit))) {
      stop(
        sprintf(
          "`%s` must be c(min, max) or c(min, max, tick_step).",
          argument
        ),
        call. = FALSE
      )
    }

    if (axis_limit[1] >= axis_limit[2]) {
      stop(sprintf("`%s` must have min < max.", argument), call. = FALSE)
    }

    if (length(axis_limit) == 3L && axis_limit[3] <= 0) {
      stop(
        sprintf("`%s` tick_step must be positive.", argument),
        call. = FALSE
      )
    }

    list(
      limits = axis_limit[1:2],
      breaks = if (length(axis_limit) == 3L) {
        seq(axis_limit[1], axis_limit[2], by = axis_limit[3])
      } else {
        NULL
      }
    )
  }

  y_axis <- parse_axis_limit(ylimit, "ylimit")
  y_limits <- y_axis$limits
  y_breaks <- y_axis$breaks

  x_axis <- parse_axis_limit(xlimit, "xlimit")
  x_limits <- x_axis$limits
  x_breaks <- x_axis$breaks

  if (!is.character(yright_endpoint) ||
      length(yright_endpoint) != 1L ||
      !yright_endpoint %in% c("min", "max")) {
    stop("`yright_endpoint` must be either 'min' or 'max'.", call. = FALSE)
  }

  make_y_scale <- function(sec_axis_arg = waiver()) {
    scale_args <- list(sec.axis = sec_axis_arg)

    if (!is.null(y_limits)) {
      scale_args$limits <- y_limits
    }

    if (!is.null(y_breaks)) {
      scale_args$breaks <- y_breaks
    } else {
      scale_args$n.breaks <- 10
    }

    do.call(scale_y_continuous, scale_args)
  }

  if (type == "monotone"){
    y_var <- "estimate_monotone"
    ci_lo_var <- "ci_lo_monotone"
    ci_hi_var <- "ci_hi_monotone"
  }
  else if (type == "raw"){
    y_var <- "estimate"
    ci_lo_var <- "ci_lo"
    ci_hi_var <- "ci_hi"
  }
  else if (type == "ve"){
    y_var <- "ve_monotone"
    ci_lo_var <- "ve_monotone_ci_lo"
    ci_hi_var <- "ve_monotone_ci_hi"
  }
  else {
    stop("`type` must be one of: 'raw', 'monotone', or 've'.")
  }

  if (!is.null(cutoffs)){
    if (!is.numeric(cutoffs) ||
        length(cutoffs) != 2L ||
        any(!is.finite(cutoffs)) ||
        cutoffs[1] > cutoffs[2]) {
      stop(
        "`cutoffs` must be two finite values: c(lower, upper), with lower <= upper.",
        call. = FALSE
      )
    }

    res[,ci_lo_var] <- pmax(res[, ci_lo_var], cutoffs[1])
    res[,ci_hi_var] <- pmin(res[, ci_hi_var], cutoffs[2])
  }

  shift_coef <- 0

  if (type == "ve"){

    if (is.null(ylimit)){
      scale_coef <- max(res[[ci_hi_var]], na.rm = TRUE)
    }
    else{
      scale_coef <- ylimit[2] - ylimit[1]
      shift_coef <- ylimit[1]
    }

  }
  else{
    if (is.null(ylimit)){
      scale_coef <- max(res[[ci_hi_var]], na.rm = TRUE)
    }
    else{
      scale_coef <- ylimit[2]
    }
  }

  if (yright_endpoint == "min"){
    yright <- min(res[[y_var]], na.rm = TRUE)
  }
  else{
    yright <- max(res[[y_var]], na.rm = TRUE)
  }

  ggthresh <- ggplot(res, aes(x = threshold, y = !!rlang::sym(y_var))) +
    geom_point(size = 0.7) +
    geom_line() +
    geom_ribbon(aes(ymin = !!rlang::sym(ci_lo_var),
                    ymax = !!rlang::sym(ci_hi_var)),
                alpha = 0.3, color = NA) +
    labs(x = "Thresholds", y = "Estimates (CI)") +
    theme_minimal() +
    xlab(xlabel) +
    ylab(ylabel)  +
    theme(plot.title = element_text(hjust = 0.5))

  if (!is.null(x_limits) && !exp10){
    x_scale_args <- list(limits = x_limits)
    if (!is.null(x_breaks)) {
      x_scale_args$breaks <- x_breaks
    }
    ggthresh <- ggthresh + do.call(scale_x_continuous, x_scale_args)
  }

  if (plot_density){

    RCDF <- function(a) {
      (sum(data[[weights]] * (data[[marker]] >= a)) /
         sum(data[[weights]]) * scale_coef) + shift_coef
    }

    RCDF <- Vectorize(RCDF)

    col <- c(col2rgb("olivedrab3"))
    col <- rgb(col[1], col[2], col[3], alpha = 255, maxColorValue = 255)

    ggthresh <- ggthresh +
      stat_function(fun = RCDF, color = col, geom = "line") +
      make_y_scale(
        sec_axis(~ (. - shift_coef) / scale_coef,
                 name = "Reverse CDF",
                 breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))
      ) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 0, hjust = 1),
            axis.text.y = element_text(angle = 0, hjust = 1))
  }
  else{
    ggthresh <- ggthresh + make_y_scale()
  }

  if (plot_endpoints){
    data_event <- data %>%
      dplyr::filter(.data[[event]] == 1 & .data[[time_var]] <= tf)

    data_event$y_inter <- approx(res$threshold, res[[y_var]],
                                 xout = data_event[,marker],
                                 yright = yright)$y
    ggthresh <- ggthresh +
      geom_point(aes(x = .data[[marker]], y = .data$y_inter),
                 data = data_event,
                 color = "blue") +
      geom_vline(xintercept = max(data_event[,marker]),
                 linetype = "dotted",
                 color = "red")

  }

  if (exp10){

    if (is.null(x_limits)) {
      min_value <- min(data[[marker]], na.rm = TRUE)
      max_value <- max(data[[marker]], na.rm = TRUE)
      x_limits <- c(min_value, max_value)

      if (ceiling(min_value) == floor(max_value) ||
          max_value - min_value < 1) {
        x_breaks <- seq(ceiling(min_value * 10) / 10,
                        floor(max_value * 10) / 10, by = 0.1)
      }
      else {
        x_breaks <- seq(ceiling(min_value), floor(max_value))
      }
    }
    else if (is.null(x_breaks)) {
      x_breaks <- pretty(x_limits)
      x_breaks <- x_breaks[x_breaks >= x_limits[1] & x_breaks <= x_limits[2]]
    }

    ggthresh <- ggthresh +
      scale_x_continuous(
        limits = x_limits,
        breaks = x_breaks,
        labels = parse(text = paste0("10^", x_breaks))
      )
  }

  if (!is.na(annotate)){
    x_annotate_loc <- max(data[,marker], na.rm = TRUE) * 0.85

    ggthresh <- ggthresh  +
      ggplot2::annotate("text",
                        x = x_annotate_loc,
                        y = max(res[,ci_hi_var]) * 0.95,
                        label = annotate)
  }

  ggthresh
}
