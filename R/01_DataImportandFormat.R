# 01_DataImportandFormat.R ####
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
aphia_ids %>% rename("taxa_name" = "name") ->aphia_ids
toc(log=TRUE)

# left_join(df_phyto0, aphia_ids, by="taxa_name") -> df_phyto
# df_phyto <- df_phyto %>%
#   relocate(aphia_id, .before = taxa_name)

# Generate Worms record and append to abundance data ####
tic("Generate Worms record and append to abundance data")
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

# Reattach rows with NA aphia_id
aphia_ids_full <- aphia_ids |> 
  left_join(record_info, by = c("aphia_id", "taxa_name"))

left_join(df_phyto0, aphia_ids_full, by="taxa_name") -> df_phyto
df_phyto <- df_phyto %>%
  relocate(aphia_id, .before = taxa_name)
toc(log=TRUE)

# Save data ####
tic("Save data")
saveRDS(janitor::clean_names(df_phyto),
        file = "outputs/Phyto_2000_2025_Step1.Rdat")
write.csv(janitor::clean_names(df_phyto),
          file="outputs/Phyto_2000_2025_Step1.csv",
          row.names = FALSE)
toc(log=TRUE)

# Tidy up ####
rm(list = ls(pattern = "^aphia"))
rm(list = ls(pattern = "^df"))
rm(record_info,unique_names)

detach("package:tidyverse", unload=TRUE)
detach("package:purrr", unload=TRUE)
detach("package:tictoc", unload=TRUE)
detach("package:worrms", unload=TRUE)
