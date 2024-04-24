
# commands from: https://github.com/jackscanlan/piperline/blob/main/jack_notes/basc_run.md

#library(renv)
library(pak)
library(targets)
library(tarchetypes)

source("_targets_packages.R")
source("R/functions.R")
source("R/themes.R")

# run pipeline
tar_make(script = "_targets_nocrew.R")
