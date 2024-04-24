
# commands from: https://github.com/jackscanlan/piperline/blob/main/jack_notes/basc_run.md

# Parse input arguments 
options <- commandArgs(trailingOnly = TRUE)
options

#library(renv)
library(pak)
library(targets)
library(tarchetypes)

source("_targets_packages_nocrew.R")
source("R/functions.R")
source("R/themes.R")

# run pipeline
tar_make(script = "_targets_nocrew.R")
