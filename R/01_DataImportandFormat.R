# 01_DataImportandFormat.R ####
# load packages ####
ld_pkgs <- c("tidyverse","tictoc", "worrms","purrr")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

df_phyto0 <- readRDS(file = "outputs/Phyto_2000_2025.Rdat")

unique_names <- unique(df_phyto0$taxa_name)

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
toc()

aphia_ids %>% rename("taxa_name" = "name") ->aphia_ids

left_join(df_phyto0, aphia_ids, by="taxa_name") -> df_phyto
df_phyto <- df_phyto %>% 
  relocate(aphia_id, .before = taxa_name)


~~~
# summarise data

df_phyto <- df_phyto0 %>% 
  mutate(taxName = if_else(
    is.na(taxon_qualifier),
    taxa_name,
    paste0(taxa_name,"_",taxon_qualifier)
  )) %>% 
  select(., -taxa_name, -taxon_qualifier)

