# Seasonal_Dig.R ----

# Set up ----
## load packages & set meta ----
ld_pkgs <- c("tidyverse","tictoc","purrr","mgcv")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

source("R/000setup.R")
source("R/helperFunctions.R")

## load GAM functions ----
# site_peak_distance_devt
source("R/functions/classify_gam_smooths.R")

# classify_gam_smooths
source("https://raw.githubusercontent.com/Defra-Data-Science-Centre-of-Excellence/ES_Restoration/refs/heads/main/classify_gam_smooths.R?token=GHSAT0AAAAAAEM5ZQ5MJEBUHYU7TXEZECHU2XIWXRQ")

# add_site_points_devt
source("https://raw.githubusercontent.com/Defra-Data-Science-Centre-of-Excellence/ES_Restoration/refs/heads/main/add_site_points_devt.R?token=GHSAT0AAAAAAEM5ZQ5M4SRCQOE4MBYAP2Z62XILCXA")

# func_plot_smth_deriv_rug_devt
source("https://raw.githubusercontent.com/Defra-Data-Science-Centre-of-Excellence/ES_Restoration/refs/heads/main/func_plot_smth_deriv_rug_devt.R?token=GHSAT0AAAAAAEM5ZQ5NSED22QPBSEXO4PZG2XILDLQ")

## load data ----
df0 <- readRDS("outputs/ZoopPhyto_long_USE.Rdat")
names(df0)

df0 %>% tidyr::as_tibble() %>%
  dplyr::filter(data_set == "Phytoplankton") %>% 
  dplyr::mutate(
    yday = lubridate::yday(sample_date),
    year = lubridate::year(sample_date),
    month = lubridate::month(sample_date)
    ) %>% #names()
  dplyr::select(
    sample_id,river_basin_district,wb,sample_date,
    yday, year, month, taxon, lifeform,abundance
  ) %>% 
  dplyr::group_by(dplyr::across(
    -abundance
    )
    ) %>% 
  dplyr::summarise(
    abundance = sum(abundance),
    .groups = "drop"
    ) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    lb = case_when(
      stringr::str_starts(taxon, "Phaeocystis") ~ "Phaeocystis",
      lifeform == "Diatom" ~ "Diatom",
      lifeform == "Dinoflagellate" ~ "Dinoflagellate",
      TRUE ~ NA
    )
  ) %>% 
  dplyr::filter(!(is.na(lb))) %>% 
  dplyr::select(-c(taxon,lifeform)) %>% 
  dplyr::group_by(dplyr::across(-abundance)) %>% 
  dplyr::summarise(
    abundance = sum(abundance),
    .groups = "drop"
  ) %>%
  ungroup() %>% 
  tidyr::pivot_wider(
    names_from = lb,
    values_from = abundance,values_fill = 0
    ) %>% 
  dplyr::filter(river_basin_district != "0") %>% 
  dplyr::filter(river_basin_district != "Western Wales") %>% 
  droplevels() -> df_plot
  
df_plot %>% 
  # dplyr::filter(river_basin_district == "North West") %>%
  # dplyr::filter(river_basin_district == "South West") %>%
  dplyr::filter(river_basin_district == "Anglian") %>%
  droplevels() %>% 
  ggplot2::ggplot(.,
                  aes(
                    x = yday,
                    y = log10(Diatom+1)
                  ))+
  geom_point(alpha = 0.4)+
  geom_smooth(
              method = "gam",
              formula = y ~ s(x, k = 3),
              se=FALSE,
              # show.legend = FALSE
              )+
  facet_wrap(.~year,scales = "free_y")

# run model by year ----
df_plot %>% 
  # dplyr::filter(river_basin_district == "North West") %>%
  # dplyr::filter(river_basin_district == "South West") %>%
  dplyr::filter(river_basin_district == "Anglian") %>%
  droplevels() %>% 
  dplyr::filter(year != 2026) %>% 
  dplyr::filter(year > 2006) %>%
  dplyr::filter(year == 2021) -> tmp
  
mod <- mgcv::gam(Diatom ~ s(yday,bs = "cr", k = 3), data = tmp)
classify_gam_smooths(mod)


#########
safe_gam <- purrr::safely(
  function(data) {
    mgcv::gam(
      Diatom ~ s(yday, bs = "cr", k = 3),
      data = data
    )
  }
)

mods <- df_plot %>%
  filter(
    year != 2026,
    year > 2006
  ) %>%
  droplevels() %>%
  group_by(river_basin_district, year) %>%
  nest() %>%
  mutate(
    fit = map(data, safe_gam),
    model = map(fit, "result"),
    error = map(fit, "error")
  ) %>%
  ungroup()

mods <- mods %>%
  mutate(
    summary = map(model, ~ if (is.null(.x)) NULL else summary(.x))
    )

mods <- mods %>%
  mutate(
    classification = map(
      model,
      ~ if (is.null(.x)) NULL else classify_gam_smooths(.x)
    )
  )
