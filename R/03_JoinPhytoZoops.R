# 03_JoinPhytoZoops.R ####
# import the current zooplankton database and append for analysis

# load packages ####
ld_pkgs <- c("tidyverse","tictoc","stringr")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

source("R/000setup.R")

#load data ####
tic("load data")
## ZOOPS
dfzoop <- read.csv((paste0(zoopfol,"processedData/zoopsAll.csv"))) %>% 
  # dfzoop <- read.csv("data/zoopsAll.csv") %>% 
  mutate(Biosys_short = if_else(
    is.na(BIOSYS.Code),
    NA_character_,
    substr(BIOSYS.Code, 1, nchar(BIOSYS.Code) - 1)
    )
    )
dfzoop$DataSet <- "Zooplankton"

## phyto
dfphyto <- readRDS(file="outputs/Phyto_2000_2025_with_Lifeforms_USE.Rdat") %>% 
  mutate(Biosys_short = if_else(
    is.na(biosys_code),
    NA_character_,
    substr(biosys_code, 1, nchar(biosys_code) - 1)
  )
  )
dfphyto$DataSet <- "Phytoplankton"
toc(log=TRUE)

tic("Trim phyto by sites within zoop data & homogenise variable names")
# Trim phyto by sites within zoop data & homogenise variable names ####

dfphyto_trm <- dfphyto %>% 
  ## keep only BIOSYS sites in phyto which are present in zoops
  semi_join(dfzoop, by = "Biosys_short")

## Prep zoop data ####
# retain 'sensible' variables
zooptrm <- dfzoop %>% dplyr::select(
  sample.date,
  Biosys_short,
  DataSet,
  Aphia.ID,
  Taxa,
  Region, 
  WBID,
  WB,
  LF02,
  Kingdom,
  Phylum,
  Class,
  Order,
  Family,
  Genus,
  Abund_m3,
  mn_carbTot_m3,
  md_carbTot_m3,
  mnCPerIndiv_ug,
  mdCPerIndiv_ug,
  mnlongMaxAxis_mm,
  mdlongMaxAxis_mm
  ) %>% 
  # create 'units' variables as zoops and phyto differ
  dplyr::mutate(
    Abund_m3_units = "Count_per_m3",
    mn_carbTot_m3_units = "tot_C_m3",
    md_carbTot_m3_units = "tot_C_m3",
    mnCPerIndiv_ug_units = "ug",
    mdCPerIndiv_ug_units = "ug",
    mnlongMaxAxis_mm_units = "mm",
    mdlongMaxAxis_mm_units = "mm"
  ) %>% 
  # rename for consistency
  dplyr::rename(
    region = Region,
    wb_id = WBID,
    wb = WB,
    sample_date = sample.date,
    aphia_id = Aphia.ID,
    taxon = Taxa,
    lifeform = LF02,
    Abundance_units = Abund_m3_units,
    mnCarbTot_units = mn_carbTot_m3_units,
    mdCarbTot_units = md_carbTot_m3_units,
    mnCPerIndiv_units = mnCPerIndiv_ug_units,
    mdCPerIndiv_units = mdCPerIndiv_ug_units,
    mnSize_units = mnlongMaxAxis_mm_units,
    mdSize_units = mdlongMaxAxis_mm_units,
    Abundance = Abund_m3,
    mnCarbTot = mn_carbTot_m3,
    mdCarbTot = md_carbTot_m3,
    mnCPerIndiv = mnCPerIndiv_ug,
    mdCPerIndiv = mdCPerIndiv_ug,
    mnSize = mnlongMaxAxis_mm,
    mdSize = mdlongMaxAxis_mm
    ) %>% janitor::clean_names(.) %>% 
  dplyr::select(., -c(mn_c_per_indiv_units,md_c_per_indiv_units,
               mn_c_per_indiv,md_c_per_indiv)) %>% 
  ## reorder columns to match
  dplyr::select(
    region,wb_id,wb,biosys_short,
    sample_date,
    data_set,
    aphia_id, taxon,
    lifeform,
    kingdom, phylum, class, order, family, genus,
    abundance,abundance_units,
    mn_carb_tot,mn_carb_tot_units,
    md_carb_tot,md_carb_tot_units,
    mn_size,mn_size_units,
    md_size,md_size_units
    ) %>% 
  # convert all variables to character to allow easier joining
  mutate(across(everything(), as.character))

## Prep phyto data ####
# retain 'sensible' variables
phytotrm <- dfphyto_trm %>% dplyr::select(
  sample_date,
  Biosys_short,
  DataSet,
  valid_aphia_id,
  name_use,
  region,
  wb_id,
  wb_name,
  phytoplankton_type,
  kingdom,
  phylum,
  class,
  order,
  family,
  genus,
  cells_per_litre_millilitre,
  tot_mn_c_pg_c,
  tot_md_c_pg_c,
  mean_vol_per_cell_um3,
  median_vol_per_cell_um3
  ) %>% 
  # create 'units' variables as zoops and phyto differ
  dplyr::mutate(
    cells_per_litre_millilitre_units = "cells_per_litre",
    tot_mn_c_pg_c_units = "pg_C",
    tot_md_c_pg_c_units = "pg_C",
    mean_vol_per_cell_um3_units = "volume_um3",
    median_vol_per_cell_um3_units = "volume_um3"
  ) %>% 
  # rename for consistency
  dplyr::rename(
    sample_date = sample_date,
    aphia_id = valid_aphia_id,
    aphia_id = valid_aphia_id,
    taxon = name_use,
    region = region,
    wb_id = wb_id,
    wb = wb_name,
    lifeform = phytoplankton_type,
    Abundance = cells_per_litre_millilitre,
    mnCarbTot = tot_mn_c_pg_c,
    mdCarbTot = tot_md_c_pg_c,
    mnSize = mean_vol_per_cell_um3,
    mdSize = median_vol_per_cell_um3,
    abundance_units = cells_per_litre_millilitre_units,
    mn_carb_tot_units = tot_mn_c_pg_c_units,
    md_carb_tot_units = tot_md_c_pg_c_units,
    
    mnSize_units = mean_vol_per_cell_um3_units,
    mdSize_units = median_vol_per_cell_um3_units
    ) %>% janitor::clean_names(.) %>% 
  ## reorder columns to match
  dplyr::select(
    region,wb_id,wb,biosys_short,
    sample_date,
    data_set,
    aphia_id, taxon,
    lifeform,
    kingdom, phylum, class, order, family, genus,
    abundance,abundance_units,
    mn_carb_tot,mn_carb_tot_units,
    md_carb_tot,md_carb_tot_units,
    mn_size,mn_size_units,
    md_size,md_size_units
  ) %>% 
  # convert all variables to character to allow easier joining
  mutate(across(everything(), as.character))

## check names match
table(names(zooptrm) == names(phytotrm))
toc(log=TRUE)

tic("Bind zoops and phyto to single df and export")
# Bind zoops and phyto to single df and export ####
dfall <- bind_rows(zooptrm,phytotrm) %>% 
  mutate(sample_date = as.Date(sample_date),
         abundance = as.numeric(abundance),
         mn_carb_tot = as.numeric(mn_carb_tot),
         md_carb_tot = as.numeric(md_carb_tot),
         mn_size = as.numeric(mn_size),
         md_size = as.numeric(md_size)
         )

# save data
saveRDS(dfall, file = "outputs/ZoopPhytoMatchingBIOSYS.Rdat")
write.csv(dfall, file = "outputs/ZoopPhytoMatchingBIOSYS.csv", row.names = FALSE)
toc(log=TRUE)

unlist(tictoc::tic.log())

# Tidy up ####
rm(list = ls(pattern = "^df"))
rm(list = ls(pattern = "^cb"))
rm(phytotrm, zooptrm, theme_use,ppi, zoopfol)

detach("package:tidyverse", unload=TRUE)
detach("package:tidyr", unload=TRUE)
detach("package:tictoc", unload=TRUE)
detach("package:janitor", unload=TRUE)
detach("package:stringr", unload=TRUE)
