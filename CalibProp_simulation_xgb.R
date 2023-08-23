# some error occured when installing sl3... This code is for trying to fix that error
# this chunk is for removing all the packages installed
#
# ip <- as.data.frame(installed.packages())
# head(ip)
# # if you use MRO, make sure that no packages in this library will be removed
# ip <- subset(ip, !grepl("MRO", ip$LibPath))
# # we don't want to remove base or recommended packages either\
# ip <- ip[!(ip[,"Priority"] %in% c("base", "recommended")),]
# # determine the library where the packages are installed
# path.lib <- unique(ip$LibPath)
# # create a vector with all the names of the packages you want to remove
# pkgs.to.remove <- ip[,1]
# head(pkgs.to.remove)
# # remove the packages
# sapply(pkgs.to.remove, remove.packages, lib = path.lib)

rm(list = ls())
graphics.off()
#install/load relevant packages
list.of.packages <- c("tidyverse","caret","xgboost","bootstrap","mgcv","remotes", "ranger", "earth")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, type="binary")
lapply(list.of.packages, library, character.only = TRUE)

source("CalibProp_simulation_ver.2.3.R")

#######################################

K_times <- 5000
my_sample_range <- c(250, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000)
overlap_const <- 1

xgboost_evaluation_results <- sim_ATE_est_with_pi_n_varying(K_times=K_times, n_samples_range = my_sample_range, overlap=overlap_const, outcome_reg_type="xgboost", ps_est_type="xgboost")
xgb_plots <- make_plots_for_results(xgboost_evaluation_results)

xgboost_evaluation_results
xgb_plots
