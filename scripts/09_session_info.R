################################################################################
# Project : RNA-Seq Differential Expression Analysis Pipeline
# Script  : 09_session_info.R
#
# Author  : Riya
#
# Purpose :
# Save R session information for reproducibility.
################################################################################

# Create output directory if needed
dir.create("docs", showWarnings = FALSE)

# Save complete session information
sink("docs/sessionInfo.txt")
sessionInfo()
sink()

message("Session information saved successfully.")