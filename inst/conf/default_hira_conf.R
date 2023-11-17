# Default harpSpatial configuration
###################################
				   
						   
lead_time            = seq(0, 36, 3)
lt_unit              = "h"
stations             = harpCore::station_list,
padding_i            = 5,
padding_j            = 5,
scores               = NULL,
obs_path             = ".",	
obsfile_template     = "obstable"		
fcst_model           = NULL,		   
fc_file_path         = ""
fc_file_template     = ""
fc_file_format       = NULL
fc_file_opts         = list()
fc_accumulation      = NULL
window_sizes         = c(0, 1, 2, 4, 8, 12, 20)
thresholds           = c(0.1, 1, 5, 10)
sqlite_path          = NULL
sqlite_file          = "harp_spatial_scores.sqlite"



