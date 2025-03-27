################################################
### SAL function precipitation verification  ###
################################################
# implementation by 
# Daan Degrauwe (RMIB)
# modified for harpSpatial by
# Alex Deckmyn (RMIB)
# March 2019

#' SAL score
#' @param fcfield The forecast field
#' @param obfield The observation field
#' @param thresh_scale Used for threshold of object boundaries: Rmax/threshscale
#' @param threshold A 2-vector containing (precipitation) thresholds for fc and obs.
#' @param same_threshold If TRUE, threshold is a single number for both fields
#' @param maxobj Maximum number of objects to identify
#' @param min_rain The minimum value to be considered. If no grid point is above, SAL is not computed.
#' @param add_objs If TRUE, the object matrices are returned as attributes of the result.
#' @return a data.frame with S, A, L values.
#' @param ... Ignored options.
#' @export
"SAL" <- function(fcfield, obfield,
                  thresh_scale = 15.,
                  min_rain = 0.1,
#                  threshold = c(max(fcfield, na.rm=TRUE), max(obfield, na.rm=TRUE))/threshScale,
                  same_threshold = FALSE,
                  maxobj = 1000,
                  add_objs = FALSE,
                  ...) {
    # check if sizes are equal
    if ( ! all(dim(fcfield) == dim(obfield)) ) {
      stop('Model and Observation precipitation field should have the same dimensions')
    }

    if (missing(min_rain)) min_rain <- harpSpatial_conf$sal_options$min_rain
    if (missing(thresh_scale)) thresh_scale <- harpSpatial_conf$sal_options$thresh_scale
    if (missing(same_threshold)) same_threshold <- harpSpatial_conf$sal_options$same_threshold

    # get max values for object thresholds
    rmax <- c(max(fcfield, na.rm=TRUE), max(obfield, na.rm=TRUE))
    if (any(rmax < min_rain)) {
       message('No precipitation in model and/or observation above given threshold ', min_rain)
       return(data.frame(S = NA_real_, A = NA_real_, L = NA_real_))
    }
    threshold <- rmax/thresh_scale

    # set common threshold, if required
    if (same_threshold) threshold <- rep(max(threshold), 2)
 
    # set NA's in both fields the same
    # FIXME: can SAL calculations really handle missing values?
    fcfield[is.na(obfield)] <- NA
    obfield[is.na(fcfield)] <- NA

    # identify objects
    fc_objects <- sal_identify_objects(fcfield, threshold = threshold[1], maxobj)
    ob_objects <- sal_identify_objects(obfield, threshold = threshold[2], maxobj)

    # combine ob and fc to get S,A,L :
    S <- 2*(fc_objects$stats$s - ob_objects$stats$s) / (fc_objects$stats$s + ob_objects$stats$s)

    A <- 2*(fc_objects$stats$a - ob_objects$stats$a) / (fc_objects$stats$a + ob_objects$stats$a)

    nx <- dim(fcfield)[1]
    ny <- dim(fcfield)[2]
    d <-  sqrt((nx - 1)^2 + (ny - 1)^2)
    L1 <- sqrt((fc_objects$stats$lx - ob_objects$stats$lx)^2 + 
               (fc_objects$stats$ly - ob_objects$stats$ly)^2)
    L2 <- 2 * abs(fc_objects$stats$lr - ob_objects$stats$lr)
    L <- (L1 + L2) / d

    scores <- data.frame("S" = S, "A" = A, "L" = L)
    if (add_objs) {
      attr(scores, "fc_objects") <- fc_objects
      attr(scores, "ob_objects") <- ob_objects
    }
    return(scores)
}

