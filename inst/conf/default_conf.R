# Default harpSPatial configuration
lead_time            = seq(0, 36, 3)
lt_unit              = "h"
by                   = "12h"
members              = NULL
fc_file_path         = ""
fc_file_template     = ""
fc_file_format       = "fa"
fc_options           = list()
fc_interp_method     = "closest"
fc_accumulation      = NULL
ob_file_path         = ""
ob_file_template     = ""
ob_file_format       = "hdf5"
ob_options           = list()
ob_interp_method     = "closest"
ob_accumulation      = "15m"
verif_domain         = NULL
use_mask             = FALSE
window_sizes         = c(0, 1, 2, 4, 8, 12, 20)
thresholds           = c(0.1, 1, 5, 10)
sqlite_path          = NULL
sqlite_file          = "harp_spatial_scores.sqlite"
return_data          = TRUE

