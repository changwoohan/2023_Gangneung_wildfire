############ 1. Packages on CRAN ############

# list of packages to use (in alphabetical order and 5 per row)
list.of.packages = c('dplyr', 'here', 'readr', 'lubridate','ggplot2',
                     'tidyverse', 'readxl', 'zoo', 'gsynth', 'openxlsx',
                     'grid','gridExtra','scales', 'purrr', 'glue',
                     'ggrepel', 'directlabels', 'geomtextpath', 'stringr', 
                     'sqldf', 'panelView', 'tidyr','DataCombine', 'grid')

# check if list of packages is installed. If not, it will install ones not yet installed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos = "https://cloud.r-project.org")

# load packages
lapply(list.of.packages, require, character.only = TRUE)
