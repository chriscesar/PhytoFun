## 00Newsetup.R

# 000.init.setup.R ####
## set up project, join site data to shapefiles

# load packages ####
ld_pkgs <- c("tidyverse","tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

tic("DATA IMPORT: Downloading, importing & loading data");print("Downloading, importing & loading data")
# load data source locations ####
source("R/dataFolders.R")

df_phyto0 <- readxl::read_xlsx("data/Phyto_2000_2025.xlsx",sheet = "PhyoData",
                               guess_max = 2000)
saveRDS(janitor::clean_names(df_phyto0), file = "outputs/Phyto_2000_2025.Rdat")
toc(log=TRUE)

df_phyto0 <- readRDS(file = "outputs/Phyto_2000_2025.Rdat")
