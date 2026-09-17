# loadWIMS_New.R ####
# load data from WIMS ####

# load packages ####
ld_pkgs <- c("tidyverse","arrow","tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)

# load meta data ----
source("R/000folders.R")

# load salty sites ----
# dfwimsSites <- read.csv(paste0(
#   wims_meta_fol,
#   "WIMS_TraC_EngWal_C3.csv")
#   )
wims_cache <- "C:/WIMS_cache/"
dfwimsSites <- readxl::read_xlsx(
  paste0("C:/WIMS_cache/",
         "Saline_Sites_WIMS_use.xlsx")
  )

# Pull data from PARQUET files ####
tic("Pull data from PARQUET files")
# # estuarine
# wims.est <- open_dataset(paste0(wims_cache,"WIMS_est"))
# # 'O:/National_Hydroecology/WIMS_cache/WIMS_freshwater')
# 
# # sea
# wims.sea <- open_dataset(paste0(wims_cache,"WIMS_sea"))
wims <- open_dataset(paste0(wims_cache,"WIMS_estuarine_sea_water"))
toc(log = TRUE)

tic("Which samples contain Phyto data?")
# # est
# est_phyto_site <- wims.est |>
#   dplyr::filter(MEAS_DETERMINAND_CODE == '4574') |> collect() |>
#   dplyr::filter(MEAS_RESULT != 0) |>
#   dplyr::select(SAMP_SMPT_USER_REFERENCE, res) |>
#   dplyr::distinct()
# 
# # sea
# sea_phyto_site <- wims.sea |>
#   dplyr::filter(MEAS_DETERMINAND_CODE == '4574') |> collect() |>
#   dplyr::filter(MEAS_RESULT != 0) |>
#   dplyr::select(SAMP_SMPT_USER_REFERENCE, res) |>
#   dplyr::distinct()

# all saline phyto
sites <- wims %>% 
  dplyr::filter(MEAS_DETERMINAND_CODE == '4574') |> collect() |>
  dplyr::filter(MEAS_RESULT != 0) |>
  dplyr::select(SAMP_SMPT_USER_REFERENCE, res) |>
  dplyr::distinct()
toc()

sites <- rbind(est_phyto_site, sea_phyto_site)
# rm(est_phyto_site, sea_phyto_site)

site.list <- data.frame(WIMSCode = sites$SAMP_SMPT_USER_REFERENCE,
                        Region = sites$res)

rm(sites)

site.list <- 
  site.list %>% 
  mutate(notation = paste(Region, WIMSCode, sep = '-'))

tic()
# hms.data.est <- 
#   wims.est |>
#   filter(notation %in% site.list$notation) %>%
#   filter(time_period != 1) %>% 
#   # filter(DATE_TIME > "2000-01-01 00:00:00") %>% 
#   # & MEAS_DETERMINAND_CODE %in%  c('0076', '0180')) %>%
#   collect()
# 
# hms.data.sea <- 
#   wims.sea |>
#   filter(notation %in% site.list$notation) %>%
#   filter(time_period != 1) %>% 
#   # filter(DATE_TIME > "2000-01-01 00:00:00") %>% 
#   # filter(time_period != 1 & MEAS_DETERMINAND_CODE %in%  c('0076', '0180')) %>%
#   collect()

hms.data <- 
  wims |>
  filter(notation %in% site.list$notation) %>%
  # filter(time_period != 1) %>% 
  # filter(DATE_TIME > "2000-01-01 00:00:00") %>% 
  # & MEAS_DETERMINAND_CODE %in%  c('0076', '0180')) %>%
  collect()

# wimsdat <- bind_rows(hms.data.est, hms.data.sea) |>
#   mutate(
#     sample_date = as_date(
#       parse_date_time(
#         as.character(DATE_TIME),
#         orders = c("ymd HMS", "ymd")
#       )
#     )
#   )

wimsdat <- hms.data %>% 
  mutate(
    sample_date = as_date(
      parse_date_time(
        as.character(DATE_TIME),
        orders = c("ymd HMS", "ymd")
      )
    )
  )

# rm(hms.data.est,hms.data.sea, site.list)
rm(site.list)
toc(log=TRUE)

###############################
###############################
###############################
###############################
###############################
# Join sites to WB ####
tic("Join sites to WB")
wimsdat %>% 
  dplyr::left_join(.,dfwimsSites %>% dplyr::select(notation,WB_ID,WB_Name,
                                                   WB_Type,
                                                   RBD),
                   by = "notation") %>% #names()
  dplyr::filter(!is.na(WB_Name)) %>% # TO CHECK: WHY ARE WE LOSTING SITES? ####
  # create WBID_date variable
  dplyr::mutate(code = paste0(WB_ID,"_",sample_date)) %>% 
  dplyr::mutate(code = str_remove_all(code, "-")) -> wimsdat_wb
rm(dfwimsSites,wimsdat)

rm(wims.est,wims.sea)
toc(log=TRUE)

# arrange regions as factor in clockwise from north east ####
wimsdat_wb %>% 
  mutate(RIVER_BASI = factor(RBD,
                             levels = c(
                               "Anglian",
                               "Humber",
                               "Thames",
                               "Severn",
                               "Northumbria",
                               "Solway Tweed",
                               "North West",
                               "Dee",
                               "South East",
                               "South West",
                               "Western Wales"
                               ))) -> wimsdat_wb
# write to RData and csv ####
tic("Write export")
saveRDS(wimsdat_wb, file="outputs/wims_all_New.RDat")
write.csv(wimsdat_wb, file="outputs/wims_all_New.csv",row.names = FALSE)
toc(log=TRUE)
# # Trim WIMS data by those '.$code' values in plankton data ####
# tic("Trim WIMS data by those '.$code' values in plankton data")
# logics <- wimsdat_wb$code %in% dfsummary_wb$code
# 
# wimsdat_wb[logics,] -> wimsdat_wb_trm
# rm(logics)
# toc(log = TRUE)
