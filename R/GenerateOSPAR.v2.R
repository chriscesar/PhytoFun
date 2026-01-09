# GenerateOSPAR.v2.R ####
# amalgamate data for OSPAR submission

# In it's current form, we convert the majority of phytoplankton taxon names into
# a form with a recognised Aphia ID (RLIST = "ERID".
# Taxa labelled as 'phytoplankton' have no Aphia ID, so get flagged as RLIST = "NT"

# load packages ####
ld_pkgs <- c("tidyverse","tictoc","sf","units")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log
tictoc::tic("TOTAL TIME")
source("R/helperFunctions.R")

# tictoc::tic("Load_data_exported_from_BOXI_Query")
# df0 <- readxl::read_xlsx("data/Phyto_2000_2025.xlsx",sheet = "PhyoData",
#                                guess_max = 2000) %>% 
#   janitor::clean_names()
# saveRDS(janitor::clean_names(df_phyto0), file = "outputs/ospar_Phyto_2000_2025.raw.Rdat")
# tictic::toc(log = TRUE)

df0 <- readRDS("outputs/ospar_Phyto_2000_2025.raw.Rdat")

# create new object with appropriate variable names ####
names(df0)

## Convert grid refs ####
tictoc::tic("Convert grid refs")
# refs <- convert_multiple_osgb(unique(df0$ngr)) %>% rename(ngr=grid_ref)
# saveRDS(refs, file = "outputs/ospar_gridRefs.raw.Rdat")
refs <- readRDS(file = "outputs/ospar_gridRefs.raw.Rdat")
tictoc::toc(log=TRUE)

tictoc::tic("Update taxon names & Aphia IDs")
# Convert taxon names refs and update to 'correct' Aphia ID ####
taxa <- unique(df0$taxa_name)

# extract Aphia IDs  ####
tic("Query AphiaID for each name")
aphia_ids <- map_dfr(taxa, function(taxon) {
  result <- tryCatch(
    {
      # Safe querying
      res <- worrms::wm_name2id(name = taxon)
      
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
toc(log=TRUE)

## extract unmatched aphia id taxa & export for matching ####
# aphia_NA <- aphia_ids %>% dplyr::filter(is.na(aphia_id))
# write.csv(aphia_NA, file = "outputs/ospar_missingAphiaRaw.csv",row.names = FALSE)
# rm(aphia_NA)
### Import missing Aphia values ####
blankAPHIA <- read.csv("outputs/MissingAphiaLookupImport.csv") #created manually
#

tictoc::tic("retrieve record details: non-na")
# Filter out rows with NA aphia_id
aphia_ids_non_na <- aphia_ids %>% 
  dplyr::filter(!is.na(aphia_id)) %>% 
  dplyr::mutate(
    RLIST = "ERID"
  )

# For each aphia_id, retrieve record details
record_info_non_na <- aphia_ids_non_na |> 
  dplyr::mutate(
    wm_data = purrr::map(aphia_id, ~ tryCatch(worrms::wm_record(.x),
                                              error = function(e) NULL))
  ) |>
  tidyr::unnest_wider(wm_data)
tictoc::toc(log=TRUE)

tictoc::tic("retrieve record details: na")
# For each aphia_id, retrieve record details
record_info_missing <- blankAPHIA |> 
  mutate(
    wm_data = purrr::map(APHIA_ID, ~ tryCatch(worrms::wm_record(.x),
                                              error = function(e) NULL))
  ) |>
  tidyr::unnest_wider(wm_data)
record_info_missing |> rename("taxa_name" = "Taxon") -> record_info_missing
tictoc::toc(log=TRUE)

# Combine WORMS records ####
tic("Combine WORMS records")

# tidy names
record_info_missing0 <- janitor::clean_names(record_info_missing) %>% 
  dplyr::select(.,-notes, -flag)
record_info_non_na0 <- janitor::clean_names(record_info_non_na) %>% dplyr::rename(taxa_name=name)

records <- rbind(record_info_non_na0,record_info_missing0) %>%
  dplyr::rename(reported_taxon = taxa_name)
write.csv(records,
          file="outputs/ospar_taxon_records.csv",row.names = FALSE)
toc(log=TRUE)

tictoc::toc(log=TRUE)

tictoc::tic("Pull out required taxon variables & append to data")
dfout <- df0

# tx <- as_tibble(dfout$taxa_name) %>% dplyr::rename(SPECI = value) %>% 
#   left_join(records %>%
#               dplyr::select(reported_taxon,rlist) %>%
#               dplyr::rename(SPECI = reported_taxon,
#                             RLIST = rlist),by = SPECI)
tx <- dfout$taxa_name %>% 
  tibble::as_tibble() %>%  
  dplyr::rename(speci = value) %>% 
  dplyr::left_join(
    records %>% 
      dplyr::select(reported_taxon, rlist,aphia_id,aphia_id_2,valid_aphia_id,
                    scientificname, valid_name) %>% 
      dplyr::rename(speci = reported_taxon,
                    rlist = rlist),
    by = "speci"
  ) %>% 
  dplyr::rename(
    SPECI = speci,
    RLIST = rlist
  )

dfout$reported_Name <- tx$SPECI
dfout$name_Type <- tx$RLIST
dfout$Aphia_ID <- tx$aphia_id
dfout$Aphia_ID_2 <- tx$aphia_id_2
dfout$SPECI <- tx$valid_aphia_id
dfout$Name_scientific <- tx$scientificname
dfout$Name_valid <- tx$valid_name
dfout$RLIST_old <- tx$RLIST
# rm(tx,records, record_info_missing,record_info_missing0, 
#    record_info_non_na, record_info_non_na0,
#    aphia_ids_non_na,aphia_ids, blankAPHIA)
tictoc::toc(log=TRUE)

tictoc::tic("Sort & join NGR data")

ngr <- as_tibble(dfout$ngr) %>% dplyr::rename(ngr=value) %>% 
  dplyr::left_join(.,refs, by = "ngr")

## append latlong values
dfout$LATIT <- round(ngr$latitude,5) # rounded to 5 dp
dfout$LONGI <- round(ngr$longitude,5) # rounded to 5 dp
tictoc::toc(log = TRUE)

tictoc::tic("Create data for export")
dfout %>% 
  ### add/calculate variables
  dplyr::mutate(
    # Reporting laboratory; 32 = EA Head Office
    # RLABO = 32,
    RLABO = "AWUK", # updated/corrected with EA PBO office (KFH)
    # Monitoring year
    MYEAR = lubridate::year(sample_date),
    # Ship or platform code
    # SHIPC = "Coastal Survey Vessel",
    SHIPC = "AA31", # changed to 'generic ship' code
    # Cruise identifier (series of sampling occasions)
    # If CRUIS is empty set year_month as  default value
    CRUIS = paste0(sprintf("%02d", month(sample_date)),"_",year(sample_date)),
    SMPNO = sample_id,
    SDATE = sample_date |>
      as.POSIXct(tz = "UTC") |>
      format("%Y%m%d"),
    # STNNO needs to be unique for each sampling event
    # Station identification /Sampling event ID
    # using BIOSYS Site ID
    # STNNO = site_id, # Needs to be unique for each sample, so append BIOSYS & Date
    ############ APPEND SMPNO to end 
    # STNNO = paste0(site_id,SDATE), #generates non-unique values (replicate samples?)
    # STNNO = paste0("ID",#ADDED VARIABLE TO ENSURE EXCEL READS AS TEXT (OTHERWISE, 0-9 CHARACTERS BEYOND LENGTH=15 DISPLAYED AS ZERO)
    #                as.character(site_id),
    #                as.character(SDATE),
    #                as.character(SMPNO)),
    STNNO = paste0(site_id,"_",SDATE),
    DTYPE = "PP",
    PARAM = "ABUNDNR",
    MUNIT = "nrcells/l",
    ALABO = "AWUK", # updated/corrected with EA PBO office (KFH)
    SLABO = "AWUK", # updated/corrected with EA PBO office (KFH)
    # RLIST = "ERID",
    RLIST = case_when(
      is.na(SPECI) ~ "NT",
      TRUE ~ "ERID"
    ),
    ####################### DONE #2 ############################
    ######### AGGREGATE BY SPECI NAME (to avoid duplicates)   ##
    ####################### DONE #2 ############################
    SPECI = as.character(SPECI),
    SPECI = case_when(
      is.na(SPECI) ~ "phytoplankton",
      TRUE ~ SPECI
    ),
    ##########################################################
    SMTYP = case_when(
      sample_method == "MARINE: Block Net" ~"BPL", ##Block Net method data entered in error
      sample_method == "MARINE: Bottle Sample" ~"BPL",
      sample_method == "MARINE: Integrated Hose Sample" ~"HOS",
      sample_method == "MARINE: Surface water sample" ~"BPL"
    ),
    PURPM = "B~S~T~E",
    ## Create 'empty' columns
    POSYS = NA_character_,
    # STATN = site_station_name,
    ## remove special characters from site names
    ### Change length station name to biosys site code ##
    STATN = site_id,
    # STATN = str_replace_all(site_station_name,"[^A-Za-z0-9 ]", ""),
    
    WADEP = NA_character_,
    EDATE = NA_character_,
    STIME = NA_character_,
    ATIME = NA_character_,
    ETIME = NA_character_,
    MNDEP = 0, # changed to zero to reflect vertical dragging up of samples
    MXDEP = 0.2, # following MB's email
    NOAGG = NA_character_,
    FNFLA = NA_character_,
    FINFL = NA_character_,
    SMVOL = NA_character_,
    WIRAN = NA_character_,
    CLMET = NA_character_,
    FLVOL = NA_character_,
    NPORT = NA_character_,
    SFLAG = NA_character_,
    STRID = NA_character_,
    SIZCL = NA_character_,
    SIZRF = NA_character_,
    MAGNI = NA_character_,
    COEFF = NA_character_,
    TRPHY = NA_character_,
    STAGE = NA_character_,
    VFLAG = NA_character_,
    QFLAG = NA_character_,
    ###################### DONE #1 ###################################
    VALUE = cells_per_litre_millilitre, #remove NA VALUES!!!!!!!!!! ##
    ###################### DONE #1 ###################################
    CPORT = NA_character_,
    SDVOL = NA_character_,
    REFSK = NA_character_,
    METST = NA_character_,
    METFP = NA_character_,
    METPT = NA_character_,
    METCX = NA_character_,
    METPS = NA_character_,
    METOA = NA_character_,
    FORML = NA_character_,
    ACCRD = NA_character_,
    ACORG = NA_character_,
    MESHS = NA_character_,
    SAREA = NA_character_,
    SPEED = NA_character_,
    PDMET = NA_character_,
    SPLIT = NA_character_,
    DURAT = NA_character_,
    ICCOD = NA_character_,
    SUBST = NA_character_,
    DEPOS = NA_character_,
    PCNAP = NA_character_,
    PRSUB = NA_character_,
    MATRX = NA_character_,
    WLTYP = NA_character_,
    MSTAT = NA_character_,
    MPROG = "CEMP~NATL"# flagged as NATIONAL MONITORING #added CEMP
  ) %>% 
  dplyr::select(
    RLABO,
    MYEAR,
    SHIPC,
    CRUIS,
    STNNO,
    LATIT,
    LONGI,
    # POSYS,
    ################### DONE #3 ###########################
    STATN, ### shorten STATION NAMES                     ##
    ################### DONE #3 ###########################
    # WADEP,
    SDATE,
    # EDATE,
    # STIME,
    # ATIME,
    # ETIME,
    DTYPE,
    SMPNO,
    MNDEP,
    MXDEP,
    # NOAGG,
    # FNFLA,
    # FINFL,
    # SMVOL,
    # WIRAN,
    # CLMET,
    # FLVOL,
    # NPORT,
    RLIST,
    SPECI,
    # Name_valid,
    # reported_Name,
    # RLIST_old,
    # SFLAG,
    # STRID,
    # SIZCL,
    # SIZRF,
    # MAGNI,
    # COEFF,
    # TRPHY,
    # STAGE,
    PARAM,
    MUNIT,
    # VFLAG,
    # QFLAG,
    VALUE,
    # CPORT,
    # SDVOL,
    ALABO,
    # REFSK,
    # METST,
    # METFP,
    # METPT,
    # METCX,
    # METPS,
    # METOA,
    # FORML,
    # ACCRD,
    # ACORG,
    SLABO,
    SMTYP,
    # MESHS,
    # SAREA,
    # SPEED,
    # PDMET,
    # SPLIT,
    # DURAT,
    # ICCOD,
    # SUBST,
    # DEPOS,
    # PCNAP,
    # PRSUB,
    # MATRX,
    # WLTYP,
    # MSTAT,
    PURPM,
    MPROG,
    ## testing only
    # site_id
  ) -> dfout_exp
tictoc::toc(log = TRUE)

# housekeeping & tidy up ####
dfout_exp %>% 
  ## arrange for consistent outputs ####
dplyr::arrange(SDATE,SMPNO,SPECI) %>% 
  ## 1) ----------- Remove NA values in VALUE ####
  dplyr::filter(.,!is.na(VALUE)) %>% 
  ## 2) ----------- Sum repeated taxa within samples
  dplyr::group_by(dplyr::across(!VALUE)) %>% 
  dplyr::summarise(VALUE = sum(VALUE), .groups = "drop") %>% ungroup() %>%
  dplyr::relocate(VALUE,.after = MUNIT ) %>%
  ## 3) ----------- Shorten station names and remove symbols
  # guide at:
  # https://www.ices.dk/data/Documents/Station/StationUserGuide.zip
  # suggests string =< 50 characters.
  # 3.1 remove special characters
  dplyr::mutate(STATN = stringr::str_replace_all(STATN,
                                                 "[^[:alnum:] ]", "")
                ) %>% 
  # 3.2 ensure double spaces are removed
  dplyr::mutate(STATN = stringr::str_squish(STATN)) %>% 
  # 3.3 truncate station names to 50 characters
  dplyr::mutate(STATN = substr(STATN, 1,50)
                ) -> df_exp_tweak

tictoc::tic("Write data")
write.csv(
  x = df_exp_tweak,
  # x = dfout_exp,
  file = paste0("outputs/OSPAR/",format(Sys.Date(), format="%Y%m%d"),"_OSPAR_Phyto_export_MASTER.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8" #ensure UTF-8 encoding
  )

saveRDS(
  df_exp_tweak,
  # dfout_exp,
  file = paste0("outputs/OSPAR/",format(Sys.Date(), format="%Y%m%d"),"_OSPAR_Phyto_export_MASTER.Rdat"),
  )
tictoc::toc(log = TRUE)

# tictoc::tic("Create Aphia_Clean")
# # Create Aphia_Clean version: retain only those taxa recorded as RLIST_old == "ERID"
# dfout_exp %>% 
#   dplyr::filter(., RLIST_old == "ERID") %>% 
#   # remove unneeded columns
#   dplyr::select(
#     -c(
#       Name_valid,
#       reported_Name,
#       RLIST_old
#     )
#   ) -> dfout_AphiaClean
# tictoc::tic("Write data")
# write.csv(x = dfout_AphiaClean,
#           file = "outputs/OSPAR_Phyto_export_Aphia.csv",
#           row.names = FALSE)
# saveRDS(dfout_AphiaClean,file = "outputs/OSPAR_Phyto_export_Aphia.Rdat")
# tictoc::toc(log = TRUE)
# tictoc::toc(log = TRUE)
# 
# tictoc::tic("Create NT version")
# # Create Aphia_Clean version: retain only those taxa recorded as RLIST_old == "NT"
# dfout_exp %>% 
#   dplyr::filter(., RLIST_old == "NT") %>% 
#   # remove unneeded columns
#   dplyr::select(
#     -c(
#       Name_valid,
#       reported_Name,
#       RLIST_old
#     )
#   ) -> dfout_NT
# tictoc::tic("Write data")
# write.csv(x = dfout_NT,
#           file = "outputs/OSPAR_Phyto_export_NT.csv",
#           row.names = FALSE)
# saveRDS(dfout_NT,file = "outputs/OSPAR_Phyto_export_NT.Rdat")
# tictoc::toc(log = TRUE)

tictoc::toc(log = TRUE)
tictoc::toc(log=TRUE)
unlist(tictoc::tic.log())
