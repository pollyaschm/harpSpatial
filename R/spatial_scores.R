# various spatial verification scores
#' Calculate spatial scores
#' @param score The score to calculate
#' @param obfield A matrix containing the observation field
#' @param fcfield A matrix containing the forecast field. Must have the same dimension as obfield.
#' @param ... Other options that may depend on the score (like scale, threshold, ...)
#' @export
spatial_scores <- function(score = NULL, obfield = NULL, fcfield = NULL, ...) {
  # TODO: add score options, plot_func and plot_opt
  score_list <- list(
                     "bias"   = list(fields = c("bias"), "func" = "score_sp_aggregated"),
                     "mse"    = list(fields = c("mse"), "func" = "score_sp_aggregated"),
#                     "gridded" = list(fields = c("bias", "mse"), "func" = "score_sp_gridded"),
                     "SAL"     = list(fields = c("S", "A", "L"), "func" = "SAL"),
                     "FSS"     = list(fields = c("threshold", "scale", "value"), "func" = "score_fss"),
                     "FSS_p"     = list(fields = c("percentile", "scale", "value"), "func" = "score_fss")
                     )
  # if called without "score", return a list of all scores
  if (is.null(score)) return(names(score_list))
  else if (!is.element(score, names(score_list))) error("Unknown score ", score, ".\n")

  # Derive table structure
  table_structure <- spatial_score_table(score_list[[score]]$fields)
  # if called without "obfield" and "fcfield", just return the table structure for the given score
  if (is.null(obfield) && is.null(fcfield)) {
    return(score_list[[score]])
  }

  # FIXME: we may be calling with options that are not recognised/used by the score
  do.call(score_list[[score]]$func, list(obfield = obfield, fcfield = fcfield, ...))
}




#' Run "fuzzy" spatial verification for 1 case
#'
#' @param obfield Observation grid.
#' @param fcfield Forecast field
#' @param thresholds A vector of thresholds
#' @param window_sizes A vector of (odd!) window sizes
#' @return A tibble with columns for threshold, window_size and various scores.
#' @export
verify_fuzzy <- function(obfield, fcfield, thresholds, window_sizes, ...) {
  # you  might as well calculate a fixed set: they're 'cheap'
  if (is.character(scores)) scores <- list(scores)
  if (any(window_sizes %% 2 != 1)) stop("Window sizes must be odd.")

  # basic preparation
  # TODO: field.type="Preciptation" ???
  # TODO: different scores *may* require different scales/thresholds?
  # we use our optimised Rcpp fastSmooth code

  nthresh <- length(thresholds)
  nwin <- length(window_sizes)

#  result <- tibble::tibble(
#    threshold = rep(thresholds, each=nwin),
#    scale     = rep(window_sizes, nthresh),
#    ets       = as.vector(vv$fuzzy$ets),
#    fss       = as.vector(vv$fss$values),
#    hk        = as.vector(vv$multi.event$hk)
#  )
  # Some other scores are probably best added in the same call, so we only calculate fractions once
  # but it may not matter so much...
  list("fss" = score_fss(obfield=obfield, fcfield=fcfield, thresholds=thresholds, window_sizes))
}

#' Run spatial verification for 1 case
#'
#' @param obfield Observation grid.
#' @param fcfield Forecast field
#' @param ... ignored
#' return A 1-row tibble of scores
#' @export
score_sp_bias <- function(obfield, fcfield, ...) {
  ## basic spatial scores that do not require a threshold, scale etc.
  ## so for a given case (date, time, leadtime), every score is a single number.
  ## we store MSE, not RMSE, because eventually we may want to sum over a period, too.
#  VXstats = c("ets", "hk", "f", "bias", "mse")
#  s1 <- SpatialVx::vxstats(obfield, fcfield, which.stats = VXstats)
  dimxy <- prod(dim(obfield))

  # put all together in a tibble
  tibble::tibble(
    bias  = sum(fcfield - obfield)/dimxy,
  )
}

#' Run spatial verification for 1 case
#'
#' @param obfield Observation grid.
#' @param fcfield Forecast field
#' @param ... ignored
#' return A 1-row tibble of scores
#' @export
score_sp_mse <- function(obfield, fcfield, ...) {
  ## basic spatial scores that do not require a threshold, scale etc.
  ## so for a given case (date, time, leadtime), every score is a single number.
  ## we store MSE, not RMSE, because eventually we may want to sum over a period, too.
#  VXstats = c("ets", "hk", "f", "bias", "mse")
#  s1 <- SpatialVx::vxstats(obfield, fcfield, which.stats = VXstats)
  dimxy <- prod(dim(obfield))

  # put all together in a tibble
  tibble::tibble(
    mse   = sum((fcfield - obfield)^2)/dimxy
  )
}

#score_sp_gridded <- function(obfield, fcfield) {
# return(fcfield - obfield)
#}

