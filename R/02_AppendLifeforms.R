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

# Load data and join ####
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

unlist(tictoc::tic.log())
