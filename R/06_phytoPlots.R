# 06_phytoPlots.R ####
# load packages ####
ld_pkgs <- c("tidyverse","tictoc","stringr")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

tic("Load data")

source("R/000setup.R")

## load phyto data
dfphyto <- readRDS(file="outputs/Phyto_2000_2025_with_Lifeforms_USE.Rdat") %>% 
  # created in script '02_AppendLifeforms.R'
  mutate(Biosys_short = if_else(
    is.na(biosys_code),
    NA_character_,
    substr(biosys_code, 1, nchar(biosys_code) - 1)
  )
  ) %>% 
  ### change river_basi to a factor
  dplyr::mutate(river_basi = factor(river_basi,
                                    levels = c(
                                      "Northumbria",
                                      "Humber",
                                      "Anglian",
                                      "Thames",
                                      "South East",
                                      "South West",
                                      "Severn",
                                      "North West",
                                      "Solway Tweed"
                                      ))) %>% 
  dplyr::mutate(wb_name = factor(wb_name,
                                 levels = unique(wb_name[order(river_basi)])
                                 ),
                Biosys_short = factor(biosys_code,
                                      levels = unique(biosys_code[order(river_basi)]))
                )

dfphyto$DataSet <- "Phytoplankton"

dfphyto %>% 
  dplyr::select(
    river_basi,
    wb_name,
    site_id,
    biosys_code,
    sample_date,
    site_station_name,
    area,
    water_body_type,
    name_use,
    phyto_lf,
    mean_vol_per_cell_um3,
    cells_per_litre_millilitre,
    tot_mn_c_pg_c_per_l
    ) -> dfphyto_trm

toc(log=TRUE)

# Summarise data by sampling event & plot ####
tic("Summarise data by sampling event & plot")

# how many samples should a WB have?
number <- 100

# dfphyto_trm %>% 
#   dplyr::select(-c(
#     name_use,
#     phyto_lf
#   )) %>% 
#   dplyr::filter(cells_per_litre_millilitre>0) %>% 
#   dplyr::group_by(across(c(
#     -cells_per_litre_millilitre,
#     -tot_mn_c_pg_c_per_l,
#     -mean_vol_per_cell_um3
#     ))) %>% 
#   dplyr::summarise(cells_per_litre_millilitre = sum(cells_per_litre_millilitre,
#                                                     na.rm = TRUE),
#                    tot_mn_c_pg_c_per_l = sum(tot_mn_c_pg_c_per_l,
#                                              na.rm = TRUE),
#                    mean_vol_per_cell_um3 = mean(mean_vol_per_cell_um3,
#                                                 na.rm = TRUE),
#                    .groups = "drop") %>% 
#   dplyr::mutate(yday = lubridate::yday(sample_date)) %>% 
#   dplyr::group_by(wb_name) %>%
#   dplyr::filter(n() >= number) %>%
#   dplyr::ungroup() %>%
#   dplyr::filter(cells_per_litre_millilitre >1) %>% 
#   dplyr::group_by(across(c(
#     -cells_per_litre_millilitre,
#     -tot_mn_c_pg_c_per_l,
#     -mean_vol_per_cell_um3
#     ))) %>% 
#   dplyr::summarise(
#     weighted_mean_vol = weighted.mean(mean_vol_per_cell_um3,
#                                       cells_per_litre_millilitre,
#                                       na.rm = TRUE)) %>% 
#   dplyr::filter(sample_date>"2008-01-01") %>%

dfphyto_trm %>%
  dplyr::select(-c(name_use, phyto_lf)) %>%
  dplyr::filter(cells_per_litre_millilitre > 0) %>%
  dplyr::group_by(across(c(
    -cells_per_litre_millilitre,
    -tot_mn_c_pg_c_per_l,
    -mean_vol_per_cell_um3
  ))) %>%
  dplyr::summarise(
    cells_per_litre_millilitre = sum(cells_per_litre_millilitre, na.rm = TRUE),
    tot_mn_c_pg_c_per_l = sum(tot_mn_c_pg_c_per_l, na.rm = TRUE),
    mean_vol_per_cell_um3 = mean(mean_vol_per_cell_um3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(yday = lubridate::yday(sample_date)) %>%
  dplyr::semi_join(
    dfphyto_trm %>%
      dplyr::filter(!is.na(sample_date)) %>%
      dplyr::distinct(wb_name, sample_date) %>%
      dplyr::count(wb_name, name = "n_dates") %>%
      dplyr::filter(n_dates >= number),
    by = "wb_name"
  ) %>%
  dplyr::filter(cells_per_litre_millilitre > 1) %>%
  dplyr::group_by(across(c(
    -cells_per_litre_millilitre,
    -tot_mn_c_pg_c_per_l,
    -mean_vol_per_cell_um3
  ))) %>%
  dplyr::summarise(
    weighted_mean_vol = weighted.mean(
      mean_vol_per_cell_um3,
      cells_per_litre_millilitre,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::filter(sample_date > "2008-01-01") %>% 
  ggplot(., aes(
    # x = yday,
    x = sample_date,
    # y = log10(tot_mn_c_pg_c_per_l+1)
    y = log10(weighted_mean_vol)
    )
    ) +
  geom_point(
    alpha = 0.25
  ) +
  facet_wrap(.~wb_name) +
  geom_smooth(
    method = "lm",
    # method="gam",
    se = FALSE,
    # formula = y ~ s(x, bs = "ds",k=10),# gam for date
    # formula = y ~ s(x, bs = "cc",k=10),# gam for day of year
    )
