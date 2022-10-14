##' @useDynLib harpSpatial
##' @importFrom Rcpp sourceCpp

# .onAttach <- function(libname, packagename) {
harpSpatial_conf <- NULL

.onLoad <- function(libname, packagename) {
#  print(paste("libname=", libname))
#  print(paste("packagename=", packagename))

  # NOTE: in the future, maybe all into "harpenv"
  # But for now, it is specific to harpSpatial
  if (!is.environment("harpSpatial_conf")) {
    harpSpatial_conf <<- new.env(parent=baseenv()) 
  }

  config_file <- file.path(libname, packagename, "conf/default_conf.R")
#  print(config_file)
  # we always start by defining all default options
  source(config_file, local=harpSpatial_conf)

  # now look for local config
  config_file <- Sys.getenv("HARP_SPATIAL_CONFIG")
  if (config_file != "") {
    harpSpatial_config(config_file)
  }

}
.onUnload <- function (libpath) {
  library.dynam.unload("harpSpatial", libpath)
}
