##' @useDynLib harpSpatial
##' @importFrom Rcpp sourceCpp

# .onAttach <- function(libname, packagename) {
# TODO: in stead of a default config file
#       we could also just have the default values in here.
harpSpatial_conf <- NULL

.onLoad <- function(libname, packagename) {
#  print(paste("libname=", libname))
#  print(paste("packagename=", packagename))

  # NOTE: in the future, maybe all into "harpenv"
  # But for now, it is specific to harpSpatial
  if (!is.environment("harpSpatial_conf")) {
    harpSpatial_conf <<- new.env(parent=baseenv()) 
  }

  # we always start by defining all default options
  config_file <- file.path(libname, packagename, "conf/default_conf.R")
#  print(config_file)
  # it should always exist /after installation/, but roxygen2 can get confused
  if (file.exists(config_file)) source(config_file, local=harpSpatial_conf)

  # now look for local config
  config_file <- Sys.getenv("HARP_SPATIAL_CONFIG")
  if (config_file != "") {
    read_spatial_config(config_file)
  }
 ###################
  # HiRA 
  if (!is.environment("harpSpatial_hira_conf")) {
    harpSpatial_hira_conf <<- new.env(parent=baseenv()) 
  }

  config_file <- file.path(libname, packagename, "conf/default_hira_conf.R")
  
  if (file.exists(config_file)) source(config_file, local=harpSpatial_hira_conf)
   
  # now look for local config
  config_file <- Sys.getenv("HARP_SPATIAL_HIRA_CONFIG")
  if (config_file != "") {
    read_hira_config(config_file)
  }				 

}
.onUnload <- function (libpath) {
  library.dynam.unload("harpSpatial", libpath)
}
