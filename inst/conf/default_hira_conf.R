# Default harpSpatial configuration
###################################

lead_time            = seq(0, 36, 3)
lt_unit              = "h"
members              = NULL
fc_file_path         = ""
fc_file_template     = ""
fc_file_format       = NULL
fc_file_opts         = list()
fc_domain            = NULL
fc_interp_method     = "closest"
fc_accumulation      = NULL

# OBS
ob_file_path         = ""
obsfile_template     = "obstable"
ob_file_format       = NULL
ob_file_opts         = list()
ob_domain            = NULL
ob_interp_method     = "closest"
ob_accumulation      = "15m"

# VERIFICATION DOMAIN
use_mask             = FALSE

# SCORE DETAILS: thresholds etc.
# NOTE: the window_sizes must be n >= 0.
#  The actual boxes have size 2*n+1
window_sizes         = c(0, 1, 2, 4, 8, 12, 20)
thresholds           = c(0.1, 1, 5, 10)

# OUTPUT
sqlite_path          = NULL
sqlite_file          = "harp_hira_scores.sqlite"



