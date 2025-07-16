# 01_DataImportandFormat_v2.R ####
# load packages ####
ld_pkgs <- c("tidyverse","tictoc", "worrms","purrr")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

# load data ####
# created in 00Newsetup.R
tic("Load data")
df_phyto0 <- readRDS(file = "outputs/Phyto_2000_2025.Rdat") 

unique_names <- unique(df_phyto0$taxa_name)
toc(log=TRUE)

# extract Aphia IDs and append to data ####
tic("Query AphiaID for each name")
aphia_ids <- map_dfr(unique_names, function(taxon) {
  result <- tryCatch(
    {
      # Safe querying
      res <- wm_name2id(name = taxon)
      
      # Ensure a tibble/data frame with expected fields
      if (length(res) == 0) {
        tibble(name = taxon, aphia_id = NA_integer_)
      } else {
        tibble(name = taxon, aphia_id = res)
      }
    },
    error = function(e) {
      tibble(name = taxon, aphia_id = NA_integer_)
    }
  )
})

aphia_ids |> rename("taxa_name" = "name") ->aphia_ids
toc(log=TRUE)

# Generate Worms records ####
tic("Generate Worms record")
# Filter out rows with NA aphia_id
aphia_ids_non_na <- aphia_ids |> filter(!is.na(aphia_id))

# For each aphia_id, retrieve record details
record_info <- aphia_ids_non_na |> 
  mutate(
    wm_data = map(aphia_id, ~ tryCatch(wm_record(.x), error = function(e) NULL))
  ) |>
  unnest_wider(wm_data)

write.csv(record_info,
          file="outputs/record_info1.csv",row.names = FALSE)
toc(log=TRUE)

# Fill in missing Aphia IDs ####
tic("Fill in missing Aphia IDs")
blankAPHIA <- read.csv("outputs/MissingAphiaLookup.csv") #created manually

# For each aphia_id, retrieve record details
record_info_missing <- blankAPHIA |> 
  mutate(
    wm_data = map(APHIA_ID, ~ tryCatch(wm_record(.x), error = function(e) NULL))
  ) |>
  unnest_wider(wm_data)
record_info_missing |> rename("taxa_name" = "Taxon") -> record_info_missing
write.csv(record_info_missing,
          file="outputs/record_info2.csv",row.names = FALSE)

toc(log=TRUE)

# Combine WORMS records ####
tic("Combine WORMS records")

# tidy names
record_info_missing <- janitor::clean_names(record_info_missing) %>% 
  dplyr::select(.,-notes, -flag)
record_info <- janitor::clean_names(record_info)

records <- rbind(record_info,record_info_missing)
write.csv(records,
          file="outputs/record_infoUSE.csv",row.names = FALSE)
toc(log=TRUE)

# Merge to abundance data & write####
tic("Merge to abundance data & write")
left_join(df_phyto0, records, by="taxa_name") -> df_phyto
df_phyto <- df_phyto %>% 
  mutate(name_use = coalesce(valid_name, taxa_name)) %>% 
  relocate(valid_aphia_id, .before = taxa_name) %>% 
  relocate(valid_name, .before = taxa_name) %>% 
  relocate(name_use, .before = valid_name)
# write.csv(df_phyto,
#           file = "outputs/Phyto_2000_2025_Step2.csv", row.names = FALSE)
saveRDS(df_phyto,
        file = "outputs/Phyto_2000_2025_Step2.Rdat")

toc(log = TRUE)

# Import and incorporate carbon data ####
tic("Import and incorporate carbon data")
dfcarb <- readxl::read_xlsx("data/Phyto carbon PML and EA_v7 Feb 2025.xlsx",
                            sheet = 1)
dfcarb %>% rename("valid_aphia_id" = `Aphia ID`) -> dfcarb

# Generate mean and median values for taxa with multiple estimates
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

# Join carbon data to phyto data and multiply by per cell values

df_phyto %>% 
  mutate(valid_aphia_id = as.character(valid_aphia_id)) %>% 
  left_join(dfcarb_summary, by="valid_aphia_id") %>% 
  rename("taxa_reported" = "taxa_name") %>% 
  relocate(taxa_reported, .after = last_col()) %>%
  relocate(valid_aphia_id, .after = last_col()) %>%
  relocate(valid_name, .after = last_col()) %>%
  relocate(name_use, .after = last_col()) %>%
  relocate(cells_per_litre_millilitre, .after = last_col()) %>% 
  mutate(tot_mn_vol_um3 = mean_vol_per_cell_um3*cells_per_litre_millilitre,
         tot_md_vol_um3 = median_vol_per_cell_um3*cells_per_litre_millilitre,
         tot_mn_C_pgC = mean_C_per_cell_pgC*cells_per_litre_millilitre,
         tot_md_C_pgC = median_C_per_cell_pgC*cells_per_litre_millilitre
         ) %>% 
  dplyr::select(-c(sample_type,sample_collector_name, sample_status,
                   analysis_amended_by, analysis_amended_date,url,authority,
                   valid_authority, citation,lsid,match_type,
                   modified)) %>% 
  ungroup() -> df_phyto_out

saveRDS(df_phyto_out, file = "outputs/Phyto_2000_2025_USE.Rdat")
write.csv(df_phyto_out, file = "outputs/Phyto_2000_2025_USE.csv",row.names = FALSE)

toc(log=TRUE)

unlist(tictoc::tic.log())

# Tidy up ####
rm(list=ls(pattern = "^df"))
rm(list=ls(pattern = "^recor"))
rm(list=ls(pattern = "^aph"))
rm(blankAPHIA,unique_names)

detach("package:tidyverse", unload=TRUE)
detach("package:tidyr", unload=TRUE)
detach("package:purrr", unload=TRUE)
detach("package:tictoc", unload=TRUE)
detach("package:worrms", unload=TRUE)
