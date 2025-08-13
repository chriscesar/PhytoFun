# PullWimsData.R ####

# Investigate Solent data

# load packages ####
ld_pkgs <- c("tidyverse","tictoc","stringr", "arrow")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

source("R/000setup.R")

# load data ####
tic("load data")
dfsummary <- read.csv(file = "outputs/Phyto_by_sample.csv")

dfwimsSites <- readxl::read_xlsx("data/all.saline.sites_WBs_EastNorth.xlsx",
                                 sheet = "SalineSitesJoined",
                                 guess_max = 10000
                                 )
toc(log=TRUE)

# Pull data from PARQUET files ####
tic("Pull data from PARQUET files")
# estuarine
wims.est <- open_dataset("C:/WIMS_cache/WIMS_est")
  # 'O:/National_Hydroecology/WIMS_cache/WIMS_freshwater')

# sea
wims.sea <- open_dataset("C:/WIMS_cache/WIMS_sea")
toc(log = TRUE)

# just show the data structure
# glimpse(wims.est)

# check which regions are there
# tic()
# wims.est |> select(wims_region) |> collect() |> distinct(wims_region)
# wims.sea |> select(wims_region) |> collect() |> distinct(wims_region) 
# toc()

# Which WIMS sites contain Phyto data?

tic("Which samples contain Phyto data?")
# est
est_phyto_site <- wims.est |>
  dplyr::filter(MEAS_DETERMINAND_CODE == '4574') |> collect() |>
  dplyr::filter(MEAS_RESULT != 0) |>
  dplyr::select(SAMP_SMPT_USER_REFERENCE, res) |>
  dplyr::distinct()

# sea
sea_phyto_site <- wims.sea |>
  dplyr::filter(MEAS_DETERMINAND_CODE == '4574') |> collect() |>
  dplyr::filter(MEAS_RESULT != 0) |>
  dplyr::select(SAMP_SMPT_USER_REFERENCE, res) |>
  dplyr::distinct()
toc()

sites <- rbind(est_phyto_site, sea_phyto_site)
# rm(est_phyto_site, sea_phyto_site)



# ----------------------------------------------------------------
# read data for identified sites

site.list <- data.frame(WIMSCode = sites$SAMP_SMPT_USER_REFERENCE,
                            Region = sites$res)

site.list <- 
  site.list %>% 
  mutate(notation = paste(Region, WIMSCode, sep = '-'))

tic()
hms.data.est <- 
  wims.est |>
  filter(notation %in% site.list$notation) %>%
  filter(time_period != 1) %>% 
  filter(DATE_TIME > "2000-01-01 00:00:00") %>% 
  # & MEAS_DETERMINAND_CODE %in%  c('0076', '0180')) %>%
  collect()

hms.data.sea <- 
  wims.sea |>
  filter(notation %in% site.list$notation) %>%
  filter(time_period != 1) %>% 
  filter(DATE_TIME > "2000-01-01 00:00:00") %>% 
  # filter(time_period != 1 & MEAS_DETERMINAND_CODE %in%  c('0076', '0180')) %>%
  collect()

wimsdat <- bind_rows(hms.data.est, hms.data.sea) |>
  mutate(
    sample_date = as_date(
      parse_date_time(
        as.character(DATE_TIME),
        orders = c("ymd HMS", "ymd")
      )
    )
  )

toc()

# Join sites to WB ####
tic("Join sites to WB")
wimsdat %>% 
  dplyr::left_join(.,dfwimsSites %>% dplyr::select(notation,WBID,WBName,Type,
                                                   EA_AREA_NA,RIVER_BASI),
                   by = "notation") %>% 
  dplyr::filter(!is.na(WBName)) %>% 
  # create WBID_date variable
  dplyr::mutate(code = paste0(WBID,"_",sample_date)) %>% 
  dplyr::mutate(code = str_remove_all(code, "-")) -> wimsdat_wb
toc(log=TRUE)

# create WBID_date variable for plankton data ####
tic("Create WBID_date variable for plankton data")
dfsummary %>% 
  dplyr::mutate(code = paste0(wbid, "_", sample_date)) %>% 
  dplyr::mutate(code = str_remove_all(code, "-")) -> dfsummary_wb

# Trim WIMS data by those '.$code' values in plankton data ####
tic("Trim WIMS data by those '.$code' values in plankton data")
logics <- wimsdat_wb$code %in% dfsummary_wb$code

wimsdat_wb[logics,] -> wimsdat_wb_trm
rm(logics)
toc(log = TRUE)

# Widen trimmed wims data ####
tic("Widen trimmed wims data")
## append units to descriptors
wimsdat_wb_trm %>%
  dplyr::mutate(dete_units = paste0(DETE_SHORT_DESC,"_",UNIT_SHORT_DESC)) %>% 
  #remove now-obsolete info
  dplyr::select(-c(
    DETE_SHORT_DESC,
    UNIT_SHORT_DESC
  )) %>% 
  # append qualifiers to result
  dplyr::mutate(
    result_use = if_else(
      MEAS_SIGN %in% c("<", ">"),
      paste0(MEAS_SIGN, MEAS_RESULT),
      as.character(MEAS_RESULT)
      )) %>% 
  dplyr::select(-c(MEAS_SIGN, MEAS_RESULT)) %>% 
  # Pivot to wider
  pivot_wider(
    names_from = dete_units,
    values_from = result_use
    ) -> wimsdat_wb_trm_w


unlist(tictoc::tic.log())

