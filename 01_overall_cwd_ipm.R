###
### Alison C. Ketz 05/31/2023
###
###

###########################################################
### Preliminaries
###########################################################

rm(list = ls())

setwd("~/Documents/ipm/cwd_ipm_sim2")

library(viridis)
library(RColorBrewer)
library(Hmisc)
library(lubridate)
library(readxl)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(xtable)
library(nimble)
library(tidyverse)
# library(dplyr)
library(lattice)
library(foreign)
# library(rgdal)
# library(rgeos)
# library(zoo)
# library(spdep)
library(parallel)
library(doParallel)
library(coda)
# library(INLA)
# library(sf)
# library(terra)
# library(splines)
library(MetBrewer)
library(ggforce)
library(tidyr)
library(ggridges)
library(ggh4x)

###########################################################
### Process number and setting seed
###########################################################

processnum <- 100
set.seed(processnum + 1000)

###########################################################
### Source summary function for posteriors
###########################################################

source("support_functions.R")

###############################################################
# Load "Data" i.e. posterior summary stats from data model fit
###############################################################

source("02_load_all_data_to_run.R")

###############################################################
# Generate simulated data
###############################################################

source("04_generate_data_noremove_icap.R")
# source("04_generate_data_equal_weight_rm.R")
# source("04_generate_data_weighted_icap.R")

###########################################################
### Setup consts etc for running the model
###########################################################

source("05_prelim_survival.R")

source("06_prelim_foi.R")

###########################################################
### Likelihoods
###########################################################

source("07_distributions.R")

# ###########################################################
# ### Functions for Efficient Calculations
# ###########################################################

source("08_calculations.R")

###########################################################
### Run model
###########################################################

source("09_modelcode.R")

###########################################################
### Run model
###########################################################

source("10_run_model.R")
# source("10_run_model_par.R")

###########################################################
### Post processing
###########################################################

source("11_post_process.R")
