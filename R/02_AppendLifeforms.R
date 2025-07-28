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

df0$sample_dateTime <- df0$sample_date
x <- substr(df0$sample_date,1,10)
x %>% as_tibble()-> x#mutate()
x %>% mutate(date = as.Date(value,)) ->y
df0$sample_date <- y$date

df0 %>%
  relocate(sample_dateTime,.after = sample_time) %>%
  relocate(sample_date) %>% 
  relocate(site_id) %>% 
  relocate(wb_name) %>% 
  relocate(river_basi) %>% 
  arrange(river_basi, wb_name, site_id, sample_date)  -> df0
rm(x,y)
toc(log=TRUE)

tic("extract Biosys code from site_station_name")
## extract BIOSYS codes from site_station_name
# extract Biosys code from site_station_name ###
site_station_name <- unique(df0$site_station_name)
biosys_code <- str_extract(site_station_name, "[A-Z]{3}[0-9]{3}[A-Z]")
# regular expression to extract sections of text with:
# [A-Z]{3} EXACTLY 3 upper case letters (e.g. 'SOL'), followed by
# [0-9]{3} EXACTLY 3 digits (e.g. 001)
# [A-Z] exactly 1 upper case letter (e.g. P)
# example: SOL001P
biosys_codes <- data.frame(site_station_name,biosys_code)

df0 %>% 
  left_join(biosys_codes, by="site_station_name") %>%
  relocate(biosys_code, .after = site_id) -> df0

rm(biosys_codes,biosys_code, site_station_name)
toc(log=TRUE)

tic("Append lifeforms")
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
  group_by(valid_aphia_id, Taxon, `Major Group`) %>%
  summarise(mean_vol_per_cell_um3=mean(as.numeric(`Volume per cell (µm3)`),
                                       na.rm = TRUE),
            median_vol_per_cell_um3=median(as.numeric(`Volume per cell (µm3)`),
                                           na.rm = TRUE),
            mean_C_per_cell_pgC=mean(as.numeric(`Carbon per cell (pgC)`),
                                     na.rm = TRUE),
            median_C_per_cell_pgC=median(as.numeric(`Carbon per cell (pgC)`),
                                         na.rm = TRUE),
            .groups = "drop"
  ) %>% ungroup() -> dfcarb_summary
write.csv(dfcarb_summary, file = "outputs/dfcarbsummary.csv",row.names = FALSE)
# df %>% 
#   dplyr::select(valid_aphia_id,name_use, size_class:protozoa_feeding) %>% 
#   count(pick(everything())) -> taxon_data
#   
# write.csv(taxon_data, file = "outputs/phyto_taxon_data.csv",
#           row.names = FALSE)
# saveRDS(taxon_data, file = "outputs/phyto_taxon_data.Rdat")
toc(log = TRUE)

unlist(tictoc::tic.log())

# Tidy up ####
rm(list = ls(pattern = "^df"))
#rm(taxon_data)

detach("package:tidyverse", unload=TRUE)
detach("package:tidyr", unload=TRUE)
detach("package:tictoc", unload=TRUE)
