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



# 
# 
# 
# 
# 
# 
# 
# 
# 
# # other data processing
# # save(hms.data, file = '~/A_RStudio_projects/H4/peripheralH4/Harmonised and Ospar/hms.test.data.rda')
# 
# # =====================================
# # annual count of orthophosphate data
# library(sf)
# library(ggplot2)
# 
# load("O:/National_Hydroecology/WIMS_cache/all.wims.sites.rda")
# load("O:/WRTS/Hydroecology/GIS/shapefiles_base_zip/areas.pf3.ms.wsea.rda")
# 
# areas <- 
#   areas.pf3.ms.wsea |>
#   filter(SEAWARD == 'No') |>
#   st_set_crs(value = 27700)
# 
# st_crs(areas) <- 27700
# 
# 
# tic()
# orthop <- 
#   wims.p |> 
#   filter(MEAS_DETERMINAND_CODE == '0180') |> 
#   collect()
# toc()
# 
# orthop2 <- 
#   orthop |>
#   mutate(year = year(DATE_TIME))
# 
# orthop.site.year <- 
#   orthop2 |>
#   group_by(notation, year) |>
#   tally()
# 
# orthop.site.year2 <-
#   orthop.site.year |>
#   left_join(all.wims.sites |> 
#               select(notation, 
#                      SMPT_TYPE, 
#                      SMPT_EASTING, 
#                      SMPT_NORTHING, 
#                      ARE_DESC)) |>
#   filter(!is.na(SMPT_EASTING)) |>
#   st_as_sf(coords = c('SMPT_EASTING', 'SMPT_NORTHING')) |>
#   st_set_crs(value = 27700)
# 
# orthop.site.year2 |>
#   filter(year >= 1980 & year <= 1999) |>
#   #  filter(year >= 2000 & year <= 2019 )
#   #  filter(year == 1980) |>
#   ggplot() +
#   geom_sf(data = areas) +
#   geom_sf(aes(colour = n), size = 0.5) +
#   scale_colour_viridis_c() +
#   facet_wrap(~year) + 
#   coord_sf(datum = sf::st_crs(27700)) +
#   theme_bw()
# 
# ggsave()
# 
# # -------------------------------
# # Data for Joe McGovern marine Ireland
# library(readxl)
# site.list.hms <- 
#   read_excel("O:/NCES Team/Evidence Synthesis/R/Harmonised and Ospar/HMS metadata.xlsx", 
#              sheet = 'sites') |>
#   mutate(notation = paste(Region, WIMSCode, sep = '-'))
# 
# my.dets <- c('0076', '9857', '4127', '0180', '0116', '0061', '7675', 
#              '9821', '7867', '9822', '0165')
# 
# tic()
# hms.data <- 
#   wims.p |>
#   filter(notation %in% site.list.hms$notation & MEAS_DETERMINAND_CODE %in%  my.dets) %>%
#   collect()
# toc()
# 
# # excel too large
# #writexl::write_xlsx(hms.data, path = 'Harmonised and Ospar/Data for Joe McGovern 2025-04-30.xlsx')
# # nanoparquet::write_parquet(hms.data, file = 'Harmonised and Ospar/Data for Joe McGovern 2025-04-30.parquet')
# 
# # ------------------------------
# # Data for CEFAS - site and det lists in one spreadsheet
# library(readxl)
# library(stringr)
# library(lubridate)
# 
# my.filename <- 
#   "O:/National_Hydroecology/WIMS_cache/scripts/EA_Legacy_data_request_Hassan_CEFAS.xlsx"
# 
# site.list.ospar <- 
#   read_excel(my.filename, 
#              sheet = 'Riverine Sites')
# det.list.ospar <-
#   read_excel(my.filename, 
#              sheet = 'Determinands') |>
#   mutate(det.code = str_replace(dtr_id, '#', ''))
# 
# ospar.data <- 
#   wims.p |>
#   filter(notation %in% site.list.ospar$site_id & MEAS_DETERMINAND_CODE %in%  det.list.ospar$det.code) |>
#   collect()
# 
# ospar.data <-
#   ospar.data |>
#   mutate(det.code.desc = paste(DETE_SHORT_DESC, MEAS_DETERMINAND_CODE, sep = ': '),
#          DATE = as.Date(DATE_TIME),
#          DATE_DEC = round(decimal_date(floor_date(DATE, 'month')), 2), 
#          YEAR = year(DATE), 
#          MONTH = month(DATE), 
#          network = 'OSPAR'
#   )
# 
# monthly.site.det <-
#   ospar.data |>
#   group_by(network, notation, DATE_DEC, det.code.desc) |>
#   summarise(n = n()) |>
#   ungroup()
# 
# # now use script X_plot sampling effort in HMS/OSPAR folder
# 
# nanoparquet::write_parquet(ospar.data, 
#                            file = 'O:/NCES Team/Evidence Synthesis/R/Harmonised and Ospar/OSPAR sampling effort/OSPAR_DATA_HASSAN_2025-08-05.parquet')
# writexl::write_xlsx(
#   ospar.data, path = 'O:/NCES Team/Evidence Synthesis/R/Harmonised and Ospar/OSPAR sampling effort/OSPAR_DATA_HASSAN_2025-08-05.xlsx')
# 

