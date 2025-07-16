# 02_AppendLifeforms.R ####

# Append phytoplankton lifeform data to data created in
# 01_DataImportandFormat_v2.R

# 01_DataImportandFormat_v2.R ####
# load packages ####
ld_pkgs <- c("tidyverse","tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

# Load data, join to lifeforms & export ####
tic("Load data")
# Load TS data created in 01_DataImportandFormat_v2.R
df0 <- readRDS("outputs/Phyto_2000_2025_USE.Rdat")

# Load lifeforms table
dflf <- readxl::read_xlsx("data/Masterlist-V7_working_EDIT_CC.xlsx",
                          sheet = "SpeciesInfoUSE") %>% 
  dplyr::filter(PlanktonType != "Fish") %>% 
  dplyr::filter(PlanktonType != "Zooplankton") %>% 
  dplyr::select(AphiaID, SizeClass, "QA Flag", PlanktonType,
                PhytoplanktonType, PhytoplanktonSize, PhytoDepth,
                PhytoFeedingMech, Toxic_Nuisance, PhytoHabitat,
                ProtozoaType, ProtozoaSize,ProtozoaHabitat,ProtozoaFeeding) %>% 
  dplyr::rename("valid_aphia_id" = "AphiaID") %>% 
  mutate(valid_aphia_id = as.character(valid_aphia_id))

# Join lifeform to time series data
df0 %>% 
  left_join(dflf, by = "valid_aphia_id") -> df

df <- janitor::clean_names(df)

saveRDS(df, file = "outputs/Phyto_2000_2025_with_Lifeforms_USE.Rdat")
write.csv(df, file = "outputs/Phyto_2000_2025_with_Lifeforms_USE.csv",
          row.names = FALSE)
toc(log = TRUE)

# Pull out & export list of taxa with their Lifeform values & Carbon values ####
tic("Pull out & export list of taxa with their Lifeform values & Carbon values")

# load carbon data 
dfcarb <- readxl::read_xlsx("data/Phyto carbon PML and EA_v7 Feb 2025.xlsx",
                            sheet = 1) %>% 
  rename("valid_aphia_id" = "Aphia ID")

dfcarb %>% 
  mutate(
    `Volume per cell (µm3)` = str_replace(`Volume per cell (µm3)`, ",", "."), # optional: handle decimal commas
    `Volume per cell (µm3)` = if_else(
      str_detect(`Volume per cell (µm3)`, "^-?\\d*\\.?\\d+$"),
      `Volume per cell (µm3)`,
      NA_character_
    ),
    `Volume per cell (µm3)` = as.numeric(`Volume per cell (µm3)`)
  ) %>% 
  mutate(
    `Carbon per cell (pgC)` = str_replace(`Carbon per cell (pgC)`, ",", "."), # optional: handle decimal commas
    `Carbon per cell (pgC)` = if_else(
      str_detect(`Carbon per cell (pgC)`, "^-?\\d*\\.?\\d+$"),
      `Carbon per cell (pgC)`,
      NA_character_
    ),
    `Carbon per cell (pgC)` = as.numeric(`Carbon per cell (pgC)`)) %>% 
  group_by(valid_aphia_id) %>%
  summarise(mean_vol_per_cell_um3=mean(as.numeric(`Volume per cell (µm3)`),
                                       na.rm = TRUE),
            median_vol_per_cell_um3=median(as.numeric(`Volume per cell (µm3)`),
                                           na.rm = TRUE),
            mean_C_per_cell_pgC=mean(as.numeric(`Carbon per cell (pgC)`),
                                     na.rm = TRUE),
            median_C_per_cell_pgC=median(as.numeric(`Carbon per cell (pgC)`),
                                         na.rm = TRUE)
  ) %>% ungroup() -> dfcarb_summary

df %>% 
  dplyr::select(valid_aphia_id,name_use, size_class:protozoa_feeding) %>% 
  count(pick(everything())) -> taxon_data
  
write.csv(taxon_data, file = "outputs/phyto_taxon_data.csv",
          row.names = FALSE)
saveRDS(taxon_data, file = "outputs/phyto_taxon_data.Rdat")
toc(log = TRUE)

unlist(tictoc::tic.log())

# Tidy up ####
rm(list = ls(pattern = "^df"))
rm(taxon_data)

detach("package:tidyverse", unload=TRUE)
detach("package:tidyr", unload=TRUE)
detach("package:tictoc", unload=TRUE)
