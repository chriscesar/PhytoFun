# 00_generate_phyto_and_zoop_data.R ####
# script to call other scripts to produce the phytoplankton and zooplankton
# data for further study

ld_pkgs <- c("tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

# source("R/00Newsetup.R")
# Imports phyo data from data exported from BIOSYS via BOXI

tic("01_DataImportandFormat_v2")
source("R/01_DataImportandFormat_v2.R")
# assigns Aphia IDs and appends carbon values
# exports taxa missing carbon data
toc(log=TRUE)

tic("02_AppendLifeforms")
source("R/02_AppendLifeforms.R")
# appends lifeform data from PLET Master (v7)
toc(log=TRUE)

tic("03_JoinPhytoZoops")
source("R/03_JoinPhytoZoops.R")
# imports zooplankton data & homogenises variable names across phyto & zoop data
# binds into a single (long) data frame & appends site metadata
# homogenises carbon content units
toc(log=TRUE)