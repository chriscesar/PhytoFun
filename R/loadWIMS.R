# loadWIMS.R ####
# load data from WIMS ####

dfwimsSites <- readxl::read_xlsx("data/all.saline.sites_WBs_EastNorth.xlsx",
                                 sheet = "SalineSitesJoined",
                                 guess_max = 10000
                                 )

# Pull data from PARQUET files ####
tic("Pull data from PARQUET files")
# estuarine
wims.est <- open_dataset("C:/WIMS_cache/WIMS_est")
# 'O:/National_Hydroecology/WIMS_cache/WIMS_freshwater')

# sea
wims.sea <- open_dataset("C:/WIMS_cache/WIMS_sea")
toc(log = TRUE)

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
rm(est_phyto_site, sea_phyto_site)

site.list <- data.frame(WIMSCode = sites$SAMP_SMPT_USER_REFERENCE,
                        Region = sites$res)

rm(sites)

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
rm(hms.data.est,hms.data.sea, site.list)
toc(log=TRUE)

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
rm(dfwimsSites,wimsdat)

rm(wims.est,wims.sea)
toc(log=TRUE)

# # Trim WIMS data by those '.$code' values in plankton data ####
# tic("Trim WIMS data by those '.$code' values in plankton data")
# logics <- wimsdat_wb$code %in% dfsummary_wb$code
# 
# wimsdat_wb[logics,] -> wimsdat_wb_trm
# rm(logics)
# toc(log = TRUE)
