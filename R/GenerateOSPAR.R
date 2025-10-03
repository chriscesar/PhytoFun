
# GenerateOSPAR.R ####
# amalgamate data for OSPAR submission

# load packages ####
ld_pkgs <- c("tidyverse","tictoc","sf","units")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

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
aphia_ids_non_na <- aphia_ids |> dplyr::filter(!is.na(aphia_id)) %>% 
  dplyr::mutate(
    RLIST = "PEG_BVOL"
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
      dplyr::select(reported_taxon, rlist) %>% 
      dplyr::rename(speci = reported_taxon,
                    rlist = rlist),
    by = "speci"
  ) %>% 
  dplyr::rename(
    SPECI = speci,
    RLIST = rlist
  )

dfout$SPECI <- tx$SPECI
dfout$RLIST <- tx$RLIST
rm(tx,records, record_info_missing,record_info_missing0, 
   record_info_non_na, record_info_non_na0,
   aphia_ids_non_na,aphia_ids, blankAPHIA)
tictoc::toc(log=TRUE)

tictoc::tic("Sort & join NGR data")

ngr <- as_tibble(dfout$ngr) %>% dplyr::rename(ngr=value) %>% 
  dplyr::left_join(.,refs, by = "ngr")

## append latlong values
dfout$LATIT <- ngr$latitude
dfout$LONGI <- ngr$longitude
tictoc::toc(log = TRUE)

tictoc::tic("Create data for export")
dfout %>% 
  ### add/calculate variables
  dplyr::mutate(
    # Reporting laboratory; 32 = EA Head Office
    RLABO = 32,
    # Monitoring year
    MYEAR = lubridate::year(sample_date),
    # Ship or platform code
    SHIPC = "Coastal Survey Vessel",
    # Cruise identifier (series of sampling occasions)
    # If CRUIS is empty set year_month as  default value
    CRUIS = paste0(sprintf("%02d", month(sample_date)),"_",year(sample_date)),
    # Station identification /Sampling event ID
    # using BIOSYS Site ID
    STNNO = site_id,
    SDATE = sample_date |>
      as.POSIXct(tz = "UTC") |>
      format("%Y%m%d"),
    DTYPE = "PP",
    SMPNO = sample_id,
    PARAM = "ABUNDNR",
    MUNIT = "nrcells/l",
    ALABO = 32,
    SLABO = 32,
    SMTYP = case_when(
      sample_method == "MARINE: Block Net" ~"BPL", ##Block Net method data entered in error
      sample_method == "MARINE: Bottle Sample" ~"BPL",
      sample_method == "MARINE: Integrated Hose Sample" ~"HOS",
      sample_method == "MARINE: Surface water sample" ~"BPL"
      ),
    PURPM = "B~S~T~E",
    ## Create 'empty' columns
    POSYS = NA_character_,
    STATN = site_station_name,
    WADEP = NA_character_,
    EDATE = NA_character_,
    STIME = NA_character_,
    ATIME = NA_character_,
    ETIME = NA_character_,
    MNDEP = NA_character_,
    MXDEP = NA_character_,
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
    VALUE = cells_per_litre_millilitre,
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
    MPROG = "NATL"# flagged as NATIONAL MONITORING
    ) %>% 
  dplyr::select(
    RLABO,
    MYEAR,
    SHIPC,
    CRUIS,
    STNNO,
    LATIT,
    LONGI,
    POSYS,
    STATN,
    WADEP,
    SDATE,
    EDATE,
    STIME,
    ATIME,
    ETIME,
    DTYPE,
    SMPNO,
    MNDEP,
    MXDEP,
    NOAGG,
    FNFLA,
    FINFL,
    SMVOL,
    WIRAN,
    CLMET,
    FLVOL,
    NPORT,
    SPECI,
    RLIST,
    SFLAG,
    STRID,
    SIZCL,
    SIZRF,
    MAGNI,
    COEFF,
    TRPHY,
    STAGE,
    PARAM,
    MUNIT,
    VFLAG,
    QFLAG,
    VALUE,
    CPORT,
    SDVOL,
    ALABO,
    REFSK,
    METST,
    METFP,
    METPT,
    METCX,
    METPS,
    METOA,
    FORML,
    ACCRD,
    ACORG,
    SLABO,
    SMTYP,
    MESHS,
    SAREA,
    SPEED,
    PDMET,
    SPLIT,
    DURAT,
    ICCOD,
    SUBST,
    DEPOS,
    PCNAP,
    PRSUB,
    MATRX,
    WLTYP,
    MSTAT,
    PURPM,
    MPROG
    ) -> dfout
tictoc::toc(log = TRUE)

tictoc::tic("Write data")
write.csv(x = dfout,
          file = "outputs/OSPAR_Phyto_export.csv",
          row.names = FALSE)
tictoc::toc(log = TRUE)

unlist(tictoc::tic.log())