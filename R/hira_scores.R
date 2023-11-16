# various spatial verification scores

#' Run "fuzzy" spatial verification for 1 case
#'
#' @param obfield Observation grid.
#' @param fcfield Forecast field
#' @param thresholds A vector of thresholds
#' @param scales A vector of (odd!) window sizes
#' @return A tibble with columns for threshold, window_size and various scores.
#' @export
hira_scores <- function(score = NULL, execute = NULL, ...) {
  # TODO: add score options, plot_func and plot_opt
  # FIXME: you MUST indicate the primary fields (e.g. threshold & scale) !
  # the index here should be comatible with the strategies.
  score_list <- list(
                     "hira_basic"  = list(index = -1, fields = c("bias","mse","mae", "count"), func = "scores_hira_basic"),
                     "hira_me"     = list(index =  0, fields = c("hit", "fa", "miss", "cr"), primary = c("threshold", "scale"), func = "scores_hira"),
                     "hira_pragm"  = list(index =  1, fields = c("bss","bs"), primary = c("threshold", "scale"), func = "scores_hira"),
                     "hira_crss"   = list(index =  2, fields = c("prs","px"), primary = c("threshold", "scale"), func = "scores_hira"),
					 "hira_pph"    = list(index =  3, fields = c("hit", "fa", "miss", "cr"), primary = c("threshold", "scale"), func = "scores_hira")
                     )


  # if called without "score", return a list of all scores
  if (is.null(score)) return(score_list)
  else if (!is.element(score, names(score_list))) stop("Unknown score ", score, ".\n")

  # Derive table structure
  # table_structure <- spatial_score_table(score_list[[score]]$fields)
  # if called without "obsvect" and "fcfield", just return the table structure for the given score
  if (is.null(execute)) {
    return(score_list[[score]])
  }

  # FIXME: we may be calling with options that are not recognised/used by the score
#  arglist <- names(as.list(args(score_list[[score]]$func)))
#  message("score function: ", score_list[[score]]$func)
#  message("argument list: ", paste(arglist, collapse=" "))

  do.call(score_list[[score]]$func, ... )
  
}

##' @export
scores_hira <- function (obsvect, fcvect,...) {
	  
   # Calculate MSE
   mse <- mean((obsvect - fcvect)^2)
   
   # Calculate MAE
   mae <- mean(abs(obsvect - fcvect))
   
   # Calculate Bias
   bias <- mean(obsvect - fcvect)
   
   
   list(basic = data.frame(mse = mse , mae = mae, bias = bias, count = length(obsvect)))

}

 

##' @export
scores_hira  <- function(obsvect, indices, fcfield, thresholds, scales,strategies, ...) {
    scores <- get_hira_scores(obsvect = obsvect,indices=indices,
	  fcfield=fcfield,thresholds=thresholds,scales=scales, strategies=strategies) 
 
 	
}
