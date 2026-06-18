# ImportPhyto_New.R ----

## Rebooted code to import and prep raw phytoplankton data from
## BIOSYS extract.  Replaces older versions
## DATE: 03/06/2026
## DATE Last Updated: 18/06/2026
## AUTHOR: Dr Christopher Cesar
## EDITORS:

# load packages & set meta ----
ld_pkgs <- c("tidyverse","tictoc", "worrms","purrr")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

source("R/000setup.R")
# remove unnecessary items from setup code
rm(list = ls(pattern = "^cb"));rm(ppi,theme_use)

# 00100 Import data extracted from BIOSYS using BOXI query ----
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

tictoc::tic("Load data")
# Load raw data from Rdat file (quicker processing) ----
df_phyto0 <- readRDS(file = "outputs/Phyto_raw_extract.Rdat")
tictoc::toc(log=TRUE)

# 00200 Format names & add carbon values ----

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
tictoc::tic("Query AphiaID for each name")
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
tictoc::toc(log=TRUE)

## Generate Worms records ----
tictoc::tic("Generate Worms record")

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
tictoc::toc(log=TRUE)

## Fill in missing Aphia IDs ----
tictoc::tic("Fill in missing Aphia IDs")
aphia_ids_na <- read.csv("outputs/MissingAphiaLookup.csv") %>%
  #created manually
  janitor::clean_names(.) %>% 
  dplyr::rename(taxa_name = taxon) %>% 
  dplyr::select(taxa_name,aphia_id)

aphia_all <- rbind(aphia_ids_non_na,aphia_ids_na)

tictoc::toc(log=TRUE)

tictoc::tic("Extract classification info from WORMS")
## Extract classification info from WORMS ####
record_info <- build_classification_wide(
  aphia_all$aphia_id,
    id_col_name = "aphia_id"
  )

record_info <- left_join(aphia_all,record_info, by = "aphia_id")
write.csv(record_info, "outputs/phyto_taxon_info.csv",row.names = FALSE)

## sanity check: are all names in our phyto data in the record_info object?
print("####################################################"); print("Do all names in phyto data match?"); print(paste0("######    ",table(df_phyto0$taxa_name %in% record_info$taxa_name)," of ",nrow(df_phyto0)));print("####################################################")

### Quick tidy
rm(aphia_all,aphia_ids,aphia_ids_na,aphia_ids_non_na)
rm(build_classification_wide)
tictoc::toc(log = TRUE)

tictoc::tic("Import carbon data & append")
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
    taxon_qualifier,
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
tictoc::toc(log=TRUE)

tictoc::tic("extract Biosys code from site_station_name")
## extract BIOSYS codes from site_station_name
# extract Biosys code from site_station_name ###
site_station_name <- unique(df_phyto$site_station_name)

df_phyto %>% dplyr::select(site_station_name,site_id) %>% 
  dplyr::mutate(biosys_code = str_extract(site_station_name,
                                   "[A-Z]{3}[0-9]{3}[A-Z]")) %>% 
  dplyr::mutate(tmp = paste0(site_station_name,
                             site_id)) %>% 
  dplyr::select(tmp, biosys_code) %>% 
  dplyr::distinct() %>% 
  mutate(biosys_code_short = if_else(
    is.na(biosys_code),
    NA_character_,
    substr(biosys_code, 1, nchar(biosys_code) - 1)
  )
  ) ->biosys_wip

# regular expression to extract sections of text with:
# [A-Z]{3} EXACTLY 3 upper case letters (e.g. 'SOL'), followed by
# [0-9]{3} EXACTLY 3 digits (e.g. 001)
# [A-Z] exactly 1 upper case letter (e.g. P)
# example: SOL001P

df_phyto %>% 
  dplyr::mutate(tmp = paste0(site_station_name,
                             site_id)) %>% 
  left_join(., biosys_wip, by="tmp") %>%
  relocate(biosys_code,biosys_code_short, .after = site_id) %>% 
  dplyr::select(-tmp) -> df_phyto

rm(biosys_wip, site_station_name)
tictoc::toc(log=TRUE)

tictoc::tic("Append lifeforms")
# 00300 Append life form data ----
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
tictoc::toc(log = TRUE)

## 00310 Explicitly state units used ----
tictoc::tic("Explicitly state units")
## keep 'old' version for error-checking
df_phyto0_old <- df_phyto

# modify 'live' version
df_phyto %>% 
  # remove rows with no abundance info
  filter(!if_all(c(cells_per_litre_millilitre,
                   colonies_per_litre_millilitre),is.na)) %>% 
  # specify type of abundance
  dplyr::mutate(abundance_units = dplyr::case_when(
    !is.na(
      cells_per_litre_millilitre) & is.na(
        colonies_per_litre_millilitre) ~ "cells per litre",
    is.na(
      cells_per_litre_millilitre) & !is.na(
        colonies_per_litre_millilitre) ~ "colonies per litre",
    is.na(cells_per_litre_millilitre) & is.na(
      colonies_per_litre_millilitre) ~ "none",
    !is.na(
      cells_per_litre_millilitre) & !is.na(
        colonies_per_litre_millilitre) ~ "error"
    )
    ) %>%
  relocate(abundance_units,
           .after=colonies_per_litre_millilitre) %>% 
  # consolidate cells and colonies into single variable
  dplyr::mutate(abundance = dplyr:: case_when(
    is.na(cells_per_litre_millilitre) ~ as.numeric(
      colonies_per_litre_millilitre),
    is.na(colonies_per_litre_millilitre) ~ as.numeric(
      cells_per_litre_millilitre),
    TRUE ~ NA
  )) %>% 
  dplyr::relocate(abundance,
                  .after=colonies_per_litre_millilitre) %>% 
  # if there are values for both cells & colonies, use 'cells' value
  dplyr::mutate(abundance = dplyr::case_when(
    abundance_units == "error" ~ cells_per_litre_millilitre,
    abundance_units != "error" ~ abundance
  )) %>% 
  # ... and update units variable
  dplyr::mutate(abundance_units = case_when(
    abundance_units == "error" ~ "cells per litre",
    abundance_units != "error" ~ abundance_units
  )) %>% 
  # remove cells and colonies variables
  dplyr::select(-cells_per_litre_millilitre,
                -colonies_per_litre_millilitre) %>% #names() #%>% 
  # mean/median cell volumes
  dplyr::rename(
    cell_vol_ind_mean = mean_vol_per_cell_um3_use,
    cell_vol_ind_median = median_vol_per_cell_um3_use
         ) %>% 
  ## convert to numeric
  dplyr::mutate(cell_vol_ind_mean = as.numeric(cell_vol_ind_mean),
                cell_vol_ind_median = as.numeric(cell_vol_ind_median)) %>% 
  # units for cell volume
  # dplyr::mutate(cell_vol_units = "um3") %>% #names()
  dplyr::mutate(cell_vol_ind_units = "um3 per individual") %>% 
  # dplyr::relocate(cell_vol_units, .after = cell_vol_ind_median) %>%
  dplyr::relocate(cell_vol_ind_units,.after = cell_vol_ind_median) %>%
  # carbon per cell units
  dplyr::rename(c_per_ind_mean = mean_c_per_cell_pg_c_use,
                c_per_ind_median = median_c_per_cell_pg_c_use) %>% #names()
  ## convert to numeric
  dplyr::mutate(
    c_per_ind_mean = as.numeric(c_per_ind_mean),
    c_per_ind_median = as.numeric(c_per_ind_median)
    ) %>% 
  # dplyr::mutate(c_units = "pg") %>% 
  dplyr::mutate(c_per_ind_units = "pg C per individual") %>% 
  dplyr::relocate(c_per_ind_units,
                  .after = c_per_ind_median) %>% #names()
  # calculate volume per sample, using abundances
  dplyr::mutate(
    cell_vol_total_mean = abundance*cell_vol_ind_mean,
    cell_vol_total_median = abundance*cell_vol_ind_median,
    ) %>% 
  dplyr::relocate(cell_vol_total_mean, cell_vol_total_median,
                  .after = cell_vol_ind_median) %>% #names()
  dplyr::mutate(cell_vol_total_sample = "um3 in sample") %>% 
  # calculate carbon per sample using abundances
  dplyr::mutate(
    c_total_mean= abundance*c_per_ind_mean,
    c_total_median= abundance*c_per_ind_median
  ) %>% #names()
  dplyr::mutate(c_total_sample_units = "pg C per litre") %>% 
  dplyr::relocate(c_total_mean, c_total_median,c_total_sample_units,
                  .after = c_per_ind_median) -> df_phyto
tictoc::toc(log=TRUE)

### NEXT:
# script: '03_JoinPhytoZoops.R
# 00400 Load, format & append zooplankton data ----
## load zoop data

tictoc::tic("load & format zoop data")
df_zoop <- read.csv((paste0(zoopfol,"processedData/zoopsAll.csv"))) %>% 
  # dfzoop <- read.csv("data/zoopsAll.csv") %>% 
  mutate(Biosys_short = if_else(
    is.na(BIOSYS.Code),
    NA_character_,
    substr(BIOSYS.Code, 1, nchar(BIOSYS.Code) - 1)
  )
  ) %>% 
  janitor::clean_names()
df_zoop$data_set <- "Zooplankton"

zoop_meta <- readxl::read_xlsx(paste0(zoopfol,
                                      "processedData/260608_MBA_Returns_Amalgamated_USE.xlsx"),
                               sheet = "SiteMeta") %>% 
  dplyr::select(-Region...8) %>% 
  dplyr::rename(Region = Region...11) %>% 
  janitor::clean_names(.)

## homogenise variable names ----
### Zoops ----
zooptrm <- df_zoop %>% dplyr::as_tibble() %>% 
  # remove NA net volumes
  dplyr::filter(!is.na(net_volume_sampled_m3)) %>%
  dplyr::select(
    pot_number,
    sample_date,
    biosys_short,
    data_set,
    aphia_id,
    taxa,
    region,
    wbid,
    wb,
    eastings,
    northings,
    lf02,
    kingdom:display_name,
    abund_m3,mn_carb_tot_m3_ug,
    md_carb_tot_m3_ug,
    mn_c_per_indiv_ug,
    md_c_per_indiv_ug,
    mnlong_max_axis_mm,
    mdlong_max_axis_mm
  ) %>% 
  # add 'empty' columns to match those in phyto data
  dplyr::mutate(
    subphylum_subdivision = NA_character_,
    phylum_division = NA_character_,
    forma = NA_character_,
    variety = NA_character_
  ) %>% 
  dplyr::relocate(forma,variety,.after = species) %>% 
  dplyr::relocate(phylum_division,.after = phylum) %>% 
  dplyr::relocate(subphylum_subdivision,.after = subphylum) %>% 
  ## create 'units' columns
  dplyr::mutate(
    abund_m3_units = "count per m3",carb_tot_units = "ugC per m3",
    c_per_ind_units = "ugC per inividual",
    mnlong_max_axis_mm_units = "mm",
    mdlong_max_axis_mm_units = "mm"
  ) %>% 
  #rename for consistency
  dplyr::rename(
    sample_id = pot_number,
    wb_id = wbid,taxon = taxa,
    lifeform = lf02,
    abundance = abund_m3,
    abundance_units = abund_m3_units,
    mn_carbTot = mn_carb_tot_m3_ug,
    md_carbTot = md_carb_tot_m3_ug,
    mn_carbInd = mn_c_per_indiv_ug,
    md_carbInd = md_c_per_indiv_ug,
    mn_size = mnlong_max_axis_mm,
    md_size = mdlong_max_axis_mm,
    mn_size_units = mnlong_max_axis_mm_units,
    md_size_units = mdlong_max_axis_mm_units
  ) %>% #names()
  janitor::clean_names() %>% 
  ## reorder columns to match
  dplyr::select(
    sample_id,
    region,
    wb_id,
    wb,
    eastings,
    northings,
    biosys_short,
    sample_date,
    data_set,
    aphia_id,
    taxon,
    lifeform,
    kingdom:subspecies,# kingdom,phylum,class,order,family,genus,
    abundance,abundance_units,
    mn_carb_tot,,
    md_carb_tot,carb_tot_units,
    mn_carb_ind,,
    md_carb_ind,c_per_ind_units,
    mn_size,mn_size_units,
    md_size,md_size_units
    ) %>% 
  dplyr::mutate(
    c_per_ind_units = "ugC per inividual"
    ) %>% 
  dplyr::relocate(c_per_ind_units, .after = md_carb_ind) %>% 
  # convert all variables to character to allow easier joining
  mutate(across(everything(), as.character)) %>% 
  dplyr::mutate(
    size_units = "mm per individual"
  ) %>% dplyr::select(
    -c(mn_size_units,md_size_units)
    ) %>% 
  as_tibble()
tictoc::toc(log=TRUE)

tictoc::tic("Format phyto data")
### Phyto ----
## Prep phyto data ####
# retain 'sensible' variables
df_phyto %>% #names()
  dplyr::select(
  sample_id,
  sample_date,
  biosys_code_short,
  data_set,
  aphia_id,
  taxa_name,
  river_basi,
  wb_id,
  wb_name,
  eastings,
  northings,
  phyto_lf,
  kingdom:variety,
  abundance,
  abundance_units,
  c_total_mean,
  c_total_median,
  c_total_sample_units,
  c_per_ind_mean,
  c_per_ind_median,
  c_per_ind_units,
  cell_vol_ind_mean,
  cell_vol_ind_median,
  cell_vol_ind_units,
  ) %>% 
  # rename for consistency
  dplyr::rename(
    region = river_basi,
    wb = wb_name,
    biosys_short = biosys_code_short,
    carb_tot_units = c_total_sample_units,
    mn_carb_ind = c_per_ind_mean,
    md_carb_ind = c_per_ind_median,
    taxon = taxa_name,
    lifeform = phyto_lf,
    mn_carb_tot = c_total_mean,
    md_carb_tot = c_total_median,
    mn_size = cell_vol_ind_mean,
    md_size = cell_vol_ind_median,
    size_units = cell_vol_ind_units,
    ) %>%
  ## reorder columns to match
  dplyr::select(
    sample_id,region,wb_id,wb,
    eastings,
    northings,biosys_short,
    sample_date,
    data_set,
    aphia_id,
    taxon,
    lifeform,
    kingdom:variety,
    abundance,abundance_units,
    mn_carb_tot,
    md_carb_tot,
    carb_tot_units,
    mn_carb_ind,
    md_carb_ind,
    c_per_ind_units,
    mn_size,
    md_size,
    size_units,
  ) %>% 
  # add 'empty' columns to match zoops
  dplyr::mutate(
    infraorder = NA_character_,
    parvorder = NA_character_,
    tribe = NA_character_,
    section = NA_character_,
    subsection = NA_character_,
    subspecies = NA_character_,
    ) %>% 
dplyr::relocate(
  infraorder,
  parvorder,
  .after = suborder
  ) %>% 
  dplyr::relocate(section,
                  subsection,
                  subspecies,
                  .after = variety) %>% 
  dplyr::relocate(
    tribe, .after = subfamily) %>% 
  # convert all variables to character to allow easier joining
  mutate(across(everything(), as.character)) -> phytotrm
tictoc::toc(log=TRUE)

tictoc::tic("Format joined data")
# sense check: Do names match?
print("DO ALL NAMES MATCH?");table(names(zooptrm) == names(phytotrm))

# Join to a single df & format ----
df_all <- rbind(phytotrm, zooptrm)

## Quick tidy
rm(df_phyto, df_zoop, phytotrm, zooptrm, zoop_meta, df_phyto0_old)

tictoc::toc(log = TRUE)

## Format numeric variables to numeric ----
### create backup version for auditing ----
tictoc::tic("Homogenise units")
df_all0 <- df_all

df_all %>% 
  dplyr::mutate(
    abundance = as.numeric(abundance),
    mn_carb_tot = as.numeric(mn_carb_tot),
    md_carb_tot = as.numeric(md_carb_tot),
    mn_carb_ind = as.numeric(mn_carb_ind),
    md_carb_ind = as.numeric(md_carb_ind),
    mn_size = as.numeric(mn_size),
    md_size = as.numeric(md_size),
    ) %>% 
  ### NEXT: Convert to consistent units using
  # script: '03_JoinPhytoZoops.R'
  #### Mean & Median CARBON PER INDIVIDUAL AS ugC ----
  dplyr::mutate(
    mn_carb_ind_as_ugC = case_when(
      c_per_ind_units == "ugC per inividual" ~ mn_carb_ind,
      c_per_ind_units == "pg C per individual" ~ mn_carb_ind/1e6,
      TRUE ~ NA_real_
      ),
    md_carb_ind_as_ugC = case_when(
      c_per_ind_units == "ugC per inividual" ~ md_carb_ind,
      c_per_ind_units == "pg C per individual" ~ md_carb_ind/1e6,
      TRUE ~ NA_real_
      )) %>% 
  #### ABUNDANCE PER M3 ----
  dplyr::mutate(
    abundance_m3 = case_when(
      abundance_units == "count per m3" ~ abundance,
      abundance_units == "cells per litre" ~ abundance*1000,
      abundance_units == "colonies per litre" ~ abundance*1000,
      TRUE ~ NA_real_
      )
    ) %>% 
  ## multiply carbon per individual by abundance/m3
  ##to show carbon per m3
  dplyr::mutate(
    mn_carb_ugC_per_m3 = abundance_m3*mn_carb_ind_as_ugC,
    md_carb_ugC_per_m3 = abundance_m3*md_carb_ind_as_ugC) -> df_all

tictoc::toc(log = TRUE)

#= Assign label region names for WBs ----
tictoc::tic("Factorise regions & WBs")
df_all %>% 
  # homogenise RBD names
  mutate(region = case_when(
    region == "Anglian" ~ "Anglian",
    region == "Humber" ~ "Humber",
    region == "NEast" ~ "Northumbria",
    region == "North West" ~ "North West",
    region == "Northumbria" ~ "Northumbria",
    region == "NWest" ~ "North West",
    region == "Severn" ~ "Severn",
    region == "Solway Tweed" ~ "Solway Tweed",
    region == "South East" ~ "South East",
    region == "South West" ~ "South West",
    region == "Southern" ~ "South East",
    region == "SWest" ~ "South West",
    region == "Thames" ~ "Thames",
    TRUE ~ NA_character_
  ))  %>% 
  mutate(
    region = factor(region, levels = c(
      "Northumbria","Humber", "Anglian",
      "Thames","South East","South West",
      "Severn","North West","Solway Tweed"
      )
      )
    ) %>% 
  dplyr::mutate(
    region_lab = case_when(
      region == "North West" ~ "NW",
      region == "Anglian" ~ "Ang",
      region == "South East" ~ "SE",
      region == "Thames" ~ "Thm",
      region == "South West" ~ "SW",
      region == "Northumbria" ~ "Nmb",
      region == "Solway Tweed" ~ "SolTw",
      region == "Humber" ~ "Hmb",
      region == "Severn" ~ "Sev",
      TRUE ~ NA_character_
      )
    ) %>% 
  dplyr::mutate(
    region_lab = factor(region_lab,
                        levels = c(
                          "Nmb", "Hmb", "Ang",
                          "Thm", "SE", "SW",
                          "Sev", "NW", "SolTw"
                          )
                        )) -> df_all

## Factorise water body names by region ----
northumbria_levels <- c(
  "Northumberland North",
  "Holy Island & Budle Bay",
  "Farne Islands to Newton Haven",
  "Northumberland South",
  "BLYTH (N)",
  "TYNE",
  "Tyne and Wear",
  "WEAR",
  "TEES",
  "Tees Coastal"
  )

humber_levels <- c(
  "Yorkshire North",
  "ESK (E)",
  "Yorkshire South",
  "HUMBER LOWER",
  "HUMBER MIDDLE",
  "HUMBER UPPER"
  )

anglian_levels <- c(
  "Wash Outer",
  "WASH INNER",
  "Lincolnshire",
  "Lincs Offshore",
  "Norfolk North",
  "Norfolk East",
  "BURE & WAVENEY & YARE & LOTHING",
  "Suffolk",
  "ORWELL",
  "STOUR (ESSEX)",
  "Harwich Approaches",
  "Essex",
  "BLACKWATER",
  "Blackwater Outer",
  "WITHAM",
  "GREAT OUSE"
  )

thames_levels <- c(
  "Thames Coastal North",
  "THAMES LOWER",
  "THAMES MIDDLE",
  "THAMES UPPER",
  "MEDWAY",
  "SWALE",
  "Murston Lakes"
  )

se_levels <- c(
  "STOUR (KENT)",
  "Kent North",
  "Whitstable Bay",
  "Thames Coastal South",
  "Kent South",
  "ROTHER",
  "Sussex East",
  "Sussex",
  "OUSE",
  "ADUR",
  "ARUN",
  "CUCKMERE",
  "PAGHAM HARBOUR",
  "CHICHESTER HARBOUR",
  "LANGSTONE HARBOUR",
  "PORTSMOUTH HARBOUR",
  "Solent",
  "SOUTHAMPTON WATER",
  "BEAULIEU RIVER",
  "LYMINGTON",
  "WESTERN YAR",
  "NEWTOWN RIVER",
  "Isle of Wight East",
  "EASTERN YAR",
  "MEDINA"
  )


sw_levels <- c(
  "Dorset / Hampshire",
  "Weymouth Bay",
  "Portland Harbour",
  "Fleet Lagoon",
  "Lyme Bay East",
  "Lyme Bay West",
  "Devon South",
  "TEIGN",
  "EXE",
  "Tor Bay",
  "DART",
  "Salcombe Harbour",
  "KINGSBRIDGE",
  "PLYMOUTH TAMAR",
  "Plymouth Sound",
  "Plymouth Coast",
  "LOOE",
  # "Fowey",  # placeholder if appears later
  "St Austell",
  "CARRICK ROADS INNER",
  "Carrick Roads Outer",
  "Fal / Helford",
  "HELFORD",
  "Penzance",
  "Lands End to Trevose Head",
  "Cornwall South",
  "Cornwall North",
  "CAMEL",
  "Barnstaple Bay",
  "TAW / TORRIDGE",
  "Bridgwater Bay",
  "Bristol Channel Inner South",
  "Bristol Channel Outer South",
  "AVON",
  "POOLE HARBOUR"
  )

severn_levels <- c(
  "SEVERN MIDDLE",
  "SEVERN UPPER",
  "BRISTOL AVON"
  )

nw_levels <- c(
  "Solway Outer South",
  "Cumbria",
  "LEVEN",
  "Morecambe Bay",
  "LUNE",
  "WYRE",
  "RIBBLE",
  "KENT",
  "Mersey Mouth",
  "MERSEY"
  )

solway_levels <- c(
  "TWEED",
  "SOLWAY"
  )

## Create a single factor across all regions to put water bodies in the 
## correct order across all regions
wb_levels <- c(
  northumbria_levels,
  humber_levels,
  anglian_levels,
  thames_levels,
  se_levels,
  sw_levels,
  severn_levels,
  nw_levels,
  solway_levels
)

rm(northumbria_levels,
   humber_levels,
   anglian_levels,
   thames_levels,
   se_levels,
   sw_levels,
   severn_levels,
   nw_levels,
   solway_levels)

df_all <- df_all %>%
  mutate(
    wb = factor(wb, levels = wb_levels)
    )
rm(wb_levels)
tictoc::toc(log = TRUE)

# Classify date variable ----
tictoc::tic("Classify dates")
df_all %>% 
  dplyr::mutate(sample_date = lubridate::date(sample_date)
                ) -> df_all
tictoc::toc(log = TRUE)

# Export data ----
tictoc::tic("Export data")

saveRDS(
  df_all,
  file = "outputs/ZoopPhyto_long_USE.Rdat"
  )

write.csv(
  df_all,
  file = "outputs/ZoopPhyto_long_USE.csv",
  row.names = FALSE
  )

tictock::toc(log=TRUE)

unlist(tictoc::tic.log())

# Tidy up ####
rm(list = ls(pattern = "^df"))
rm(ppi,GISfol,theme_use,unique_names,zoopfol)

detach("package:tidyverse", unload=TRUE)
detach("package:tidyr", unload=TRUE)
detach("package:tictoc", unload=TRUE)
detach("package:worrms", unload=TRUE)
detach("package:purrr", unload=TRUE)
