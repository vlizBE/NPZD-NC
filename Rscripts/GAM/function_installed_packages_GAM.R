# functions to load packages that are needed for NPZ model (Steven Pint)

install_needed_packages_GAM <- function(){ 
  not_installed <- setdiff(c("mgcv",
                             "lubridate",
                             "ggplot2",
                             "stats",
                             "xts",
                             "lwdataexplorer",
                             "plyr",
                             "dplyr",
                             "devtools")
                           , rownames(installed.packages()))
  
  if (length(not_installed) != 0) {
    if (lwdataexplorer %in% not_installed){
      install.packages(not_installed[which(not_installed != "lwdataexplorer")])
      devtools::install_github("lifewatch/lwdataexplorer")
      } else {
    install.packages(not_installed)
      }
  } else {
    print("All packages that are needed are installed")
  }
}

