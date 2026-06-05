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

tic("Load data")
# Load raw data from Rdat file (quicker processing) ----
df_phyto0 <- readRDS(file = "outputs/Phyto_raw_extract.Rdat")
toc(log=TRUE)

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
print("####################################################");print("Do all names in phyto data match?");print(paste0("######    ",table(df_phyto0$taxa_name %in% record_info$taxa_name)));print("####################################################")

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
toc(log=TRUE)

tic("Append lifeforms")
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
toc(log = TRUE)

## 00310 Explicitly state units used ----
tictoc::tic("Explicitly state units")
df_phyto %>% 
  # remove rows with no abundance info
  filter(!if_all(c(cells_per_litre_millilitre,
                   colonies_per_litre_millilitre),is.na)) %>% 
  # specify type of abundance
  dplyr::mutate(abundance_units = dplyr::case_when(
    !is.na(
      cells_per_litre_millilitre) & is.na(
        colonies_per_litre_millilitre) ~ "cells_per_l",
    is.na(
      cells_per_litre_millilitre) & !is.na(
        colonies_per_litre_millilitre) ~ "colonies_per_l",
    is.na(cells_per_litre_millilitre) & is.na(
      colonies_per_litre_millilitre) ~ "none",
    !is.na(
      cells_per_litre_millilitre) & !is.na(
        colonies_per_litre_millilitre) ~ "error"
    )) %>%
  relocate(abundance_units, .after=colonies_per_litre_millilitre) %>% 
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
    abundance_units == "error" ~ "cells_per_l",
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
  dplyr::mutate(cell_vol_units = "um3") %>% #names()
  dplyr::relocate(cell_vol_units, .after = cell_vol_ind_median) %>%
  # carbon per cell units
  dplyr::rename(c_per_ind_mean = mean_c_per_cell_pg_c_use,
                c_per_ind_median = median_c_per_cell_pg_c_use) %>% #names()
  ## convert to numeric
  dplyr::mutate(
    c_per_ind_mean = as.numeric(c_per_ind_mean),
    c_per_ind_median = as.numeric(c_per_ind_median)
    ) %>% 
  dplyr::mutate(c_units = "pg") %>% 
  dplyr::relocate(c_units, .after = c_per_ind_median) %>% #names()
  # calculate volume per sample, using abundances
  dplyr::mutate(
    cell_vol_total_mean = abundance*cell_vol_ind_mean,
    cell_vol_total_median = abundance*cell_vol_ind_median,
    ) %>% 
  dplyr::relocate(cell_vol_total_mean, cell_vol_total_median,
                  .after = cell_vol_ind_median) %>% #names()
  # calculate carbon per sample using abundances
  dplyr::mutate(
    c_total_mean= abundance*c_per_ind_mean,
    c_total_median= abundance*c_per_ind_median
  ) %>% #names()
  dplyr::relocate(c_total_mean, c_total_median,
                  .after = c_per_ind_median) -> df_phyto
toc(log=TRUE)

##################################
##################################
## FROM HERE####
##################################
##################################


### NEXT:
# script: '03_JoinPhytoZoops.R
# 00400 Load, format & append zooplankton data ----
## load zoop data
tic("load zoop data")
df_zoop <- read.csv((paste0(zoopfol,"processedData/zoopsAll.csv"))) %>% 
  # dfzoop <- read.csv("data/zoopsAll.csv") %>% 
  mutate(Biosys_short = if_else(
    is.na(BIOSYS.Code),
    NA_character_,
    substr(BIOSYS.Code, 1, nchar(BIOSYS.Code) - 1)
  )
  )
df_zoop$data_set <- "Zooplankton"

zoop_meta <- readxl::read_xlsx(paste0(zoopfol,
                                      "processedData/260114_MBA_Returns_Amalgamated_USE.xlsx"),
                               sheet = "SiteMeta") %>% 
  dplyr::select(-Region...8) %>% 
  dplyr::rename(Region = Region...11) %>% 
  janitor::clean_names(.)

## homogenise variable names ----
### Zoops ----
zooptrm <- df_zoop %>%
  dplyr::select(
    Pot.Number,
    sample.date,
    Biosys_short,
    data_set,
    Aphia.ID,
    Taxa,
    Region,
    WBID,
    WB,
    Eastings,
    Northings,
    LF02,
    Kingdom, Phylum,Class,Order,Family,Genus,
    Abund_m3,
    mn_carbTot_m3_ug,
    md_carbTot_m3_ug,
    mnCPerIndiv_ug,
    mdCPerIndiv_ug,
    mnlongMaxAxis_mm,
    mdlongMaxAxis_mm
  ) %>% 
  ## create 'units' columns
  dplyr::mutate(
    Abund_m3_units = "Count_per_m3",
    mn_carbTot_units = "ugC_per_m3",
    md_carbTot_units = "ugC_per_m3",
    mnCPerIndiv_ug_units = "ugC_per_inividual",
    mdCPerIndiv_ug_units = "ugC_per_inividual",
    mnlongMaxAxis_mm_units = "mm",
    mdlongMaxAxis_mm_units = "mm"
  ) %>% 
  #rename for consistency
  dplyr::rename(
    sample_id = Pot.Number,
    region = Region,
    wb_id = WBID,
    wb = WB,
    eastings=Eastings,
    northings=Northings,
    sample_date = sample.date,
    aphia_id = Aphia.ID,
    taxon = Taxa,
    lifeform = LF02,
    abundance = Abund_m3,
    abundance_units = Abund_m3_units,
    mn_carbTot = mn_carbTot_m3_ug,
    md_carbTot = md_carbTot_m3_ug,
    mn_carbTot_units = mn_carbTot_units,
    md_carbTot_units = md_carbTot_units,
    mn_carbInd = mnCPerIndiv_ug,
    mn_carbInd_units = mnCPerIndiv_ug_units,
    md_carbInd = mdCPerIndiv_ug,
    md_carbInd_units = mdCPerIndiv_ug_units,
    mn_size = mnlongMaxAxis_mm,
    md_size = mdlongMaxAxis_mm,
    mn_size_units = mnlongMaxAxis_mm_units,
    md_size_units = mdlongMaxAxis_mm_units
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
    kingdom,phylum,class,order,family,genus,
    abundance,abundance_units,
    mn_carb_tot,mn_carb_tot_units,
    md_carb_tot,md_carb_tot_units,
    mn_carb_ind,mn_carb_ind_units,
    md_carb_ind,md_carb_ind_units,
    mn_size,mn_size_units,
    md_size,md_size_units
  ) %>% 
  # convert all variables to character to allow easier joining
  mutate(across(everything(), as.character))

### Phyto ----
## Prep phyto data ####
# retain 'sensible' variables
df_phyto %>% names()
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
  kingdom,phylum,class,order,family,genus,
  cells_per_litre_millilitre,
  tot_mn_c_pg_c_per_l,
  tot_md_c_pg_c_per_l,
  mean_c_per_cell_pg_c,
  median_c_per_cell_pg_c,
  mean_vol_per_cell_um3,
  median_vol_per_cell_um3,
  mean_c_per_cell_pg_c,
  median_c_per_cell_pg_c
) %>% names()
  ## create 'units' columns
  dplyr::mutate(
    cells_per_litre_millilitre_units = "cells_per_litre",
    tot_mn_c_pg_c_per_l_units = "pg_C_per_litre",
    tot_md_c_pg_c_per_l_units = "pg_C_per_litre",
    mean_c_per_cell_pg_c_units = "pg_C_per_cell",
    median_c_per_cell_pg_c_units = "pg_C_per_cell",
    mean_vol_per_cell_um3_units = "individual_cell_volume_um3",
    median_vol_per_cell_um3_units = "individual_cell_volume_um3"
  ) %>% 
  # rename for consistency
  dplyr::rename(
    region = river_basi,
    wb_id = wb_id,
    wb = wb_name,
    biosys_short = biosys_code_short,
    sample_date=sample_date,
    data_set = DataSet,
    aphia_id = valid_aphia_id,
    taxon = name_use,
    lifeform = phyto_lf,
    kingdom = kingdom, phylum = phylum, class = class,
    order = order, family = family, genus = genus,
    abundance = cells_per_litre_millilitre,
    abundance_units = cells_per_litre_millilitre_units,
    mn_carb_tot = tot_mn_c_pg_c_per_l,
    mn_carb_tot_units = tot_mn_c_pg_c_per_l_units,
    md_carb_tot = tot_md_c_pg_c_per_l,
    md_carb_tot_units = tot_md_c_pg_c_per_l_units,
    mn_carb_ind = mean_c_per_cell_pg_c,
    mn_carb_ind_units = mean_c_per_cell_pg_c_units,
    md_carb_ind = median_c_per_cell_pg_c,
    md_carb_ind_units = median_c_per_cell_pg_c_units,
    mn_size = mean_vol_per_cell_um3,
    mn_size_units = mean_vol_per_cell_um3_units,
    md_size = median_vol_per_cell_um3,
    md_size_units = median_vol_per_cell_um3_units
  ) %>% 
  janitor::clean_names() %>% 
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
    kingdom,phylum,class,order,family,genus,
    abundance,abundance_units,
    mn_carb_tot,mn_carb_tot_units,
    md_carb_tot,md_carb_tot_units,
    mn_carb_ind,mn_carb_ind_units,
    md_carb_ind,md_carb_ind_units,
    mn_size,mn_size_units,
    md_size,md_size_units
  ) %>% 
  # convert all variables to character to allow easier joining
  mutate(across(everything(), as.character))

### NEXT:
# script: '05_homogeniseUnits.R'

