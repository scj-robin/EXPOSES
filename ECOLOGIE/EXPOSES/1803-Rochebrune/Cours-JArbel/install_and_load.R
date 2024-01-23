needed_packages <- c("ggplot2", "hexbin", "viridis", "gridExtra", 
                     "ggpubr", "rgl", "reshape2", "dplyr", "DPpackage", "BNPdensity")
new_packages <- needed_packages[
  !(needed_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) 
  install.packages(new_packages)
lapply(needed_packages, require, character.only = TRUE)