# This script incorporates the trajectory-fitting scripts
# according to a simple model of exponential growth to carrying capacity.
# It includes the Python scripts used to infer the model parameters
# as well as the model outputs.
# It then combines the trajectory-fitting outputs with other species trajectory information.
source("workflow/analysis/fitColonizationTrajectories/fitColonizationTrajectories.R")
rm(list=ls())
