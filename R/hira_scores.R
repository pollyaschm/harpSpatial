# various spatial verification scores

#' Run "fuzzy" spatial verification for 1 case
#'
#' @param obfield Observation grid.
#' @param fcfield Forecast field
#' @param thresholds A vector of thresholds
#' @param scales A vector of (odd!) window sizes
#' @return A tibble with columns for threshold, window_size and various scores.
#' @export
hira_scores <- function(score = NULL, obsvect = NULL, indices = NULL,  fcfield = NULL, 
         thresholds = NULL, scales = NULL, stratigies = NULL) {
  # TODO: add score options, plot_func and plot_opt
  # FIXME: you MUST indicate the primary fields (e.g. threshold & scale) !
  score_list <- list(
                     "bias"   = list(fields = c("bias"), "func" = "scores_sp_basic", "plot_func" = "plot_basic"),
                     "mse"    = list(fields = c("mse"),  "func" = "scores_sp_basic", "plot_func" = "plot_basic"),
                     "mae"    = list(fields = c("mae"),  "func" = "scores_sp_basic", "plot_func" = "plot_basic"),
#                 
                     "NACT"    = list(fields = c("hit", "fa", "miss", "cr"), primary = c("threshold", "scale"),
                                      "func" = "scores_sp_neighborhood", "plot_func" = "plot_nact")
#                     , "FSS_p"     = list(fields = c("percentile", "scale", "fss"), "func" = "score_fss", "plot_func" = "plot_fss")
                     )

  # if called without "score", return a list of all scores
  if (is.null(score)) return(score_list)
  else if (!is.element(score, names(score_list))) stop("Unknown score ", score, ".\n")

  # Derive table structure
  # table_structure <- spatial_score_table(score_list[[score]]$fields)
  # if called without "obsvect" and "fcfield", just return the table structure for the given score
  if (is.null(obsvect) && is.null(fcfield)) {
    return(score_list[[score]])
  }

  # FIXME: we may be calling with options that are not recognised/used by the score
#  arglist <- names(as.list(args(score_list[[score]]$func)))
#  message("score function: ", score_list[[score]]$func)
#  message("argument list: ", paste(arglist, collapse=" "))

  myargs <- c(list(obsvect = obsvect, fcfield = fcfield), list(...))
  do.call(score_list[[score]]$func, myargs)
  
}
