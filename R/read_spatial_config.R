# Configuration file management

#' Read local harpSpatial configuration files
#' @param config_file The name of the config file to be read.
#'   If NULL, the current configuration is printed.
#' @export
read_spatial_config <- function(config_file = NULL) {
  if (is.null(config_file)) {
    sapply(ls(harpSpatial_conf), function(x) get(x, envir=harpSpatial_conf))
  } else {
    if (!file.exists(config_file)) stop("Config file ", config_file, " not found")
    message("Loading harpSpatial config file ", config_file)
    source(config_file, local=harpSpatial_conf)
    sapply(ls(harpSpatial_conf), function(x) get(x, envir=harpSpatial_conf))
  }
}





