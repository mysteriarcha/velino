# Script to load libraries required for the analysis of Velino vegetation

# Detach all previously loaded packages
loaded_pkgs <- paste0("package:", names(sessionInfo()$otherPkgs))
if(loaded_pkgs != "package:"){
  invisible(lapply(loaded_pkgs, base::detach, character.only = TRUE, unload = TRUE, force = TRUE))
}

# Check what are the installed packagers by the user:
inst_pkgs <- rownames(installed.packages())

# Install pak package manager (allows to install directly from
# CRAN or github):
if(!("pak" %in% inst_pkgs)) install.packages("pak")

# Make the list of required packages:
pkg_list <-
  list(
    # Basic data manipulation and visualization:
    "tidyverse",  # General data handling
    "magrittr",   # For extra pipe operators
    "gt",         # For tabulation and plotting functions
    
    # Geographic information:
    # "raster",     # Manipulation of raster data
    # "terra",      # Updated raster package
    "sf",         # Geometries
    "tmap",       # Map making
    
    # Data analysis:
    "vegan",      # General vegetation ecology
    "betapart",   # For beta-diversity metrics
    "modEvA",     # For model evaluation
    "AICcmodavg", # idem
    "betareg",    # For beta regression
    "MASS",       # For negative binomial glms
    "DHARMa",     # For spatial autocorrelation
    "qqplotr"     # For qqplot analysis
    
  )

# Install packages that are not in the user's system
invisible(
  lapply(
    pkg_list,
    function(x) if(!(x %in% inst_pkgs)) pak::pkg_install(x)
  )
)

# Load the packages
invisible(
  lapply(
    pkg_list,
    library, 
    character.only = TRUE,
    quietly = TRUE
  )
)

# Clean the environment
rm(list = ls())

print("All packages installed and loaded successfully!")

load_libs <- TRUE