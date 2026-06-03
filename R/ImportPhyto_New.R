# ImportPhyto_New.R ----

## Rebooted code to import and prep raw phytoplankton data from
## BIOSYS extract.  Replaces older versions
## DATE: 03/06/2026
## DATE Last Updated: 
## AUTHOR: Dr Christopher Cesar
## EDITORS:

# load packages & set meta ----
ld_pkgs <- c("tidyverse","tictoc", "worrms","purrr")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

source("R/000setup.R")

# Import data extracted from BIOSYS using BOXI query ----
# tic("DATA IMPORT: Downloading, importing & loading data")
# ## Takes time. Only needs running when data are updated
# ## Can be commented out
# df_phyto0 <- readxl::read_excel("data/Phyto_2000_2026.xlsx",
#                                 sheet = "PhyoData",
#                                 col_types = c(
#                                   "text", "text", "text", "text",
#                                   "text", "text", "date", "text",
#                                   "text", "text", "numeric", "text",
#                                   "text", "text", "numeric", "text",
#                                   "text", "text", "text", "text",
#                                   "text", "text", "text", "text", "text",
#                                   "text", "text", "text", "date",
#                                   "text", "date", "text", "text",
#                                   "text", "text", "text", "text",
#                                   "text", "text", "numeric", "numeric"
#                                   ))
# saveRDS(janitor::clean_names(df_phyto0),
#         file = "outputs/Phyto_raw_extract.Rdat")
# toc(log=TRUE)

tic("Load data")
# Load raw data from Rdat file (quicker processing) ----
df_phyto0 <- readRDS(file = "outputs/Phyto_raw_extract.Rdat")
toc(log=TRUE)

# Format names & add carbon values ----

# from the phyto data and use worrms::wm_name2id to assign Aphia ID values to taxa.
# Remove NA values (unassigned taxa) and map WORMS record info to each Aphia match.
# Export 'NA' Aphia values for manual match ups.
# Combine both records into a single 'records' object.
# Left join the taxon record info to the phytoplankton abundance data.
# Write out objects.
# + Import carbon-by-taxon data.  Based on PML & EA estimates. Where a taxon has >1 value, mean and median values are calculated.
# Carbon values appended to phyto data and abundance values multiplied by carbon content per cell values.
# + Summary of taxa with 'missing' carbon values exported.
# + Data exported.

unique_names <- unique(df_phyto0$taxa_name)

## Extract & append Aphia IDs ----
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

## Generate Worms records ----
tic("Generate Worms record")

# Helper function: build classification wide table for a given vector of Aphia IDs
build_classification_wide <- function(aphia_vec,
                                      id_col_name = "aphia_id",
                                      names_prefix = "") {
  # Fetch classification list for each Aphia ID
  class_list <- tibble::tibble(!!id_col_name := aphia_vec) %>%
    dplyr::mutate(
      wm_class = purrr::map(
        !!rlang::sym(id_col_name),
        ~ tryCatch(worrms::wm_classification(.x), error = function(e) NULL)
      )
    ) %>%
    tidyr::unnest(wm_class) %>%
    # Some calls may return NULL; drop those
    dplyr::filter(!is.na(AphiaID) | !is.na(scientificname) | !is.na(rank)) %>%
    dplyr::mutate(rank = tolower(rank)) %>%
    # avoid rare duplicated ranks; keep first occurrence
    dplyr::group_by(!!rlang::sym(id_col_name), rank) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    # Pivot rank -> columns, values = scientificname
    tidyr::pivot_wider(
      id_cols = !!rlang::sym(id_col_name),
      names_from = rank,
      values_from = scientificname,
      names_prefix = names_prefix
    )
  class_list
  }

# Remove rows with NA aphia_id
aphia_ids_non_na <- aphia_ids |> filter(!is.na(aphia_id))

# For each aphia_id, retrieve record details
# record_info <- aphia_ids_non_na |> 
#   mutate(
#     wm_data = map(aphia_id, ~ tryCatch(wm_classification(.x),
#                                        error = function(e) NULL))
#     ) |>
#   unnest_wider(wm_data)

# # For each aphia_id, retrieve record details
# record_info <- build_classification_wide(
#   aphia_ids_non_na$aphia_id,
#   id_col_name = "aphia_id"
# )
# 
# record_info <- left_join(aphia_ids,record_info, by = "aphia_id")
# 
# write.csv(record_info,
#           file="outputs/record_info1_non_na.csv",row.names = FALSE)
toc(log=TRUE)

## Fill in missing Aphia IDs ----
tic("Fill in missing Aphia IDs")
aphia_ids_na <- read.csv("outputs/MissingAphiaLookup.csv") %>%
  #created manually
  janitor::clean_names(.) %>% 
  dplyr::rename(taxa_name = taxon) %>% 
  dplyr::select(taxa_name,aphia_id)

aphia_all <- rbind(aphia_ids_non_na,aphia_ids_na)

toc(log=TRUE)

tic("Extract classification info from WORMS")
## Extract classification info from WORMS ####
record_info <- build_classification_wide(
  aphia_all$aphia_id,
    id_col_name = "aphia_id"
  )

record_info <- left_join(aphia_all,record_info, by = "aphia_id")
write.csv(record_info, "outputs/phyto_taxon_info.csv",row.names = FALSE)

## sanity check: are all names in our phyto data in the record_info object?
table(df_phyto0$taxa_name %in% record_info$taxa_name)

### Quick tidy
rm(aphia_all,aphia_ids,aphia_ids_na,aphia_ids_non_na)
rm(build_classification_wide)
toc(log = TRUE)

tic("Import carbon data & append")
## Import carbon data & append to aphia ID info ----
dfcarb_summary <- readxl::read_xlsx("data/PhytoCarbon_blankFill.xlsx",
                                    sheet = "out") %>% 
  # remove NA values and convert to numeric for joining
  dplyr::filter(valid_aphia_id != "NA") %>% 
  dplyr::mutate(aphia_id = as.numeric(valid_aphia_id)) %>% 
  dplyr::select(-valid_aphia_id)

aphia_carb <- dplyr::left_join(record_info,
                 dfcarb_summary,
                 by = "aphia_id"
                 )

# append taxonomy & carbon data to phyto sample data
df_phyto0 %>% 
  left_join(.,aphia_carb, by = "taxa_name") %>%# View()
  dplyr::select(-c(
    "...2",region,area,water_body,water_body_type,
    catchment,wfd_waterbody_id,analysis_amended_by,
    analysis_amended_date
    )) %>% 
  janitor::clean_names(.) %>% 
  # reorder columns
  relocate(river_basi,wb_id,wb_name,wb_type,wb_hm_designation,
             wb_typology,surveillan,
             distance_to_wb,eastings,northings) %>% 
  relocate(
    parent_name_usage_id,
    aphia_id,
    kingdom,
    subkingdom,
    infrakingdom,
    phylum,
    phylum_division,
    subphylum,
    subphylum_subdivision,
    infraphylum,
    parvphylum,
    gigaclass,
    superclass,
    class,
    subclass,
    infraclass,
    subterclass,
    superorder,
    order,
    suborder,
    superfamily,
    family,
    subfamily,
    genus,
    subgenus,
    species,
    forma,
    variety,
    taxa_name,
    cells_per_litre_millilitre,
    colonies_per_litre_millilitre,
    mean_vol_per_cell_um3_use,
    median_vol_per_cell_um3_use,
    mean_c_per_cell_pg_c_use,
    median_c_per_cell_pg_c_use,
    carbon_value_match,
    .after=last_col()
  ) %>% 
  dplyr::mutate(data_set = "Phytoplankton")-> df_phyto

## tidy up
rm(df_phyto0, aphia_carb,dfcarb_summary,record_info)
toc(log=TRUE)

tic("extract Biosys code from site_station_name")
## extract BIOSYS codes from site_station_name
# extract Biosys code from site_station_name ###
site_station_name <- unique(df_phyto$site_station_name)

df_phyto %>% dplyr::select(site_station_name,site_id) %>% 
  dplyr::mutate(biosys_code = str_extract(site_station_name,
                                   "[A-Z]{3}[0-9]{3}[A-Z]")) %>% 
  dplyr::mutate(tmp = paste0(site_station_name,
                             site_id)) %>% 
  dplyr::select(tmp, biosys_code) %>% 
  dplyr::distinct()->biosys_wip

# regular expression to extract sections of text with:
# [A-Z]{3} EXACTLY 3 upper case letters (e.g. 'SOL'), followed by
# [0-9]{3} EXACTLY 3 digits (e.g. 001)
# [A-Z] exactly 1 upper case letter (e.g. P)
# example: SOL001P

df_phyto %>% 
  dplyr::mutate(tmp = paste0(site_station_name,
                             site_id)) %>% 
  left_join(., biosys_wip, by="tmp") %>%
  relocate(biosys_code, .after = site_id) %>% 
  dplyr::select(-tmp) -> df_phyto

rm(biosys_wip, site_station_name)
toc(log=TRUE)

tic("Append lifeforms")
# Append life form data ----
lf_tmp <- readxl::read_xlsx("data/Masterlist-V7_working_EDIT_CC.xlsx",
                            sheet = "SpeciesInfoUSE") %>% 
  dplyr::filter(PlanktonType != "Fish") %>% 
  dplyr::filter(PlanktonType != "Zooplankton") %>% 
  dplyr::select(AphiaID, SizeClass, "QA Flag", PlanktonType,
                PhytoplanktonType, PhytoplanktonSize, PhytoDepth,
                PhytoFeedingMech, Toxic_Nuisance, PhytoHabitat,
                ProtozoaType, ProtozoaSize,ProtozoaHabitat,ProtozoaFeeding,
                PhytoLF) %>% 
  dplyr::rename("valid_aphia_id" = "AphiaID") %>% 
  mutate(valid_aphia_id = as.character(valid_aphia_id),
         aphia_id = as.numeric(valid_aphia_id))

# Join lifeform to time series data
df_phyto %>% 
  left_join(.,lf_tmp, by = "aphia_id") -> df_phyto

df_phyto <- janitor::clean_names(df_phyto)
rm(lf_tmp)
toc(log = TRUE)


### NEXT:
# script: '03_JoinPhytoZoops.R

### NEXT:
# script: '05_homogeniseUnits.R'

