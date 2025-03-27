#####################################
# DEFINITION OF VERIFICATION TABLES #
#####################################
# spatial_tables. Column names may not have a space or full stop.
# anyway, we hardcode per table

# return table description for spatial verification score tables
## @param tab a table name ("basic" or "fuzzy")
# @colnames Character vector giving the names of all the columns needed to describe the score (like c("threshold", "scale"), or c("S", "A", "L") 
# not exported
spatial_score_table <- function(template) {
  # NOTE: if we assume that we have different SQLite files for every parameter
  #       we don't really need to add the "prm" column.
  #       BUT: if we ever move to a full SQL database, we might want it anyway.
  # NOTE: we decided to switch fcdate back to unix time, so fctime is no longer needed
  # NOTE: leadtime should be in seconds. Always. Makes it easy to do date calculations.
  standard_fields <- c(
    "model"    = "CHARACTER",
    "prm"      = "CHARACTER",
    "fcdate"   = "REAL",
#    "fctime"   = "REAL",
    "leadtime" = "REAL"
  )
#  score_fields <- structure(rep("REAL", length(score_names)), names=score_names)
  score_fields <- rep("REAL", length(template$fields))
  names(score_fields) <- template$fields
  primary_fields <- rep("REAL", length(template$primary))
  names(primary_fields) <- template$primary

  list(
    fields = c(standard_fields, primary_fields, score_fields),
    primary = c(names(standard_fields), template$primary)
  )
}

# TODO: info table, ...

