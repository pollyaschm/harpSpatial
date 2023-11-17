# various HiRA scores

#' Run "HiRA" spatial verification for 1 case
#'
#' @param score One of HiRA Scores 
#' @param execute If True then fire the score function
#' @return A tibble with columns for obsvect, fcstvect fcfield, threshold, scales ...
#' not exported 
hira_scores <- function(score = NULL, execute = NULL, ...) {
  # TODO: add score options, plot_func and plot_opt
  # FIXME: you MUST indicate the primary fields (e.g. threshold & scale) !
  # the index here should be comatible with the strategies.
  score_list <- list(
                     "hira_bias"   = list(index = -1, fields = c("bias", "count"), primary = c("scale"), func = "scores_hira_basic"),
                     "hira_mse"    = list(index = -1, fields = c("mse", "count"),  primary = c("scale"), func = "scores_hira_basic"),
                     "hira_mae"    = list(index = -1, fields = c("mae", "count"),  primary = c("scale"), func = "scores_hira_basic"),
                     "hira_me"     = list(index =  0, fields = c("hit", "fa", "miss", "cr"), primary = c("threshold", "scale", "count"), func = "scores_hira"),
                     "hira_pragm"  = list(index =  1, fields = c("bss","bs"), primary = c("threshold", "scale", "count"), func = "scores_hira"),
                     "hira_crss"   = list(index =  2, fields = c("prs","px"), primary = c("threshold", "scale"), func = "scores_hira"),
					 "hira_pph"    = list(index =  3, fields = c("hit", "fa", "miss", "cr"), primary = c("threshold", "scale", "count"), func = "scores_hira")
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

#' Run "HiRA" spatial verification for 1 case
#'
#' @param obsvect One of HiRA Scores 
#' @param fcvect If True then fire the score function
#' @return A list of bias, mse and mae.
#' not exported  
scores_hira <- function (obsvect, fcvect,...) {
	  
   # Calculate MSE
   mse <- mean((obsvect - fcvect)^2)
   
   # Calculate MAE
   mae <- mean(abs(obsvect - fcvect))
   
   # Calculate Bias
   bias <- mean(obsvect - fcvect)
   
   
   list(basic = data.frame(mse = mse , mae = mae, bias = bias, count = length(obsvect)))

}

 
#' Run "HiRA" scores for 1 case
#'
#' @param obsvect One of HiRA Scores 
#' @param fcvect If True then fire the score function
#' @return A list of bias, mse and mae.
#' not exported 
scores_hira  <- function(obsvect, indices, fcfield, thresholds, scales,strategies, ...) {
    scores <- get_hira_scores(obsvect = obsvect,indices=indices,
	  fcfield=fcfield,thresholds=thresholds,scales=scales, strategies=strategies) 
 
 	
}

#' Run "HiRA" basic scores for 1 case
#'
#' @param obsvect One of HiRA Scores 
#' @param fcvect If True then fire the score function
#' @return A list of bias, mse and mae.
#' not exported 
scores_hira_basic  <- function(obsvect, indices, fcfield, scales, ...) {
    scores <- get_hira_basic_scores(obsvect = obsvect,indices=indices,
	  fcfield=fcfield,scales=scales) 

}
