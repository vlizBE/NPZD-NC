 # functions to load packages that are needed for NPZ model (Steven Pint)

install_needed_packages <- function(){ 
not_installed <- setdiff(c("ggplot2", 
                           "xts", 
                           "reshape2", 
                           "lubridate", 
                           "stats", 
                           "plyr", 
                           "dplyr", 
                           "RColorBrewer", 
                           "viridis", 
                           "ggpubr", 
                           "parallel", 
                           "doParallel", 
                           "foreach", 
                           "knitr")
                         , rownames(installed.packages()))

 if (length(not_installed) != 0) {
   install.packages(not_installed)
 } else {
   print("All packages that are needed are installed")
 }
}

