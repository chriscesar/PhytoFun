# Seasonal_Dig.R ----

# Set up ----
## load packages & set meta ----
ld_pkgs <- c("tidyverse","tictoc","purrr","mgcv","seas","broom",
             "forcats","scales")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

source("R/000setup.R")
source("R/helperFunctions.R")

## load GAM functions ----
# site_peak_distance_devt
source("R/functions/classify_gam_smooths.R")

# # classify_gam_smooths
# source("https://raw.githubusercontent.com/Defra-Data-Science-Centre-of-Excellence/ES_Restoration/refs/heads/main/classify_gam_smooths.R?token=GHSAT0AAAAAAEM5ZQ5MJEBUHYU7TXEZECHU2XIWXRQ")
# 
# # add_site_points_devt
# source("https://raw.githubusercontent.com/Defra-Data-Science-Centre-of-Excellence/ES_Restoration/refs/heads/main/add_site_points_devt.R?token=GHSAT0AAAAAAEM5ZQ5M4SRCQOE4MBYAP2Z62XILCXA")
# 
# # func_plot_smth_deriv_rug_devt
# source("https://raw.githubusercontent.com/Defra-Data-Science-Centre-of-Excellence/ES_Restoration/refs/heads/main/func_plot_smth_deriv_rug_devt.R?token=GHSAT0AAAAAAEM5ZQ5NSED22QPBSEXO4PZG2XILDLQ")

## load data ----
df0 <- readRDS("outputs/ZoopPhyto_long_USE.Rdat")
names(df0)

# df0 %>% tidyr::as_tibble() %>%
#   dplyr::filter(data_set == "Phytoplankton") %>% 
#   dplyr::mutate(
#     yday = lubridate::yday(sample_date),
#     year = lubridate::year(sample_date),
#     month = lubridate::month(sample_date)
#     ) %>% #names()
#   dplyr::select(
#     sample_id,river_basin_district,wb,sample_date,
#     yday, year, month, taxon, lifeform,abundance
#   ) %>% 
#   
#   # Calculate total abundance of ALL phytoplankton in each sample
#   dplyr::group_by(
#     sample_id, river_basin_district, wb, sample_date,
#     yday, year, month
#   ) %>%
#   dplyr::mutate(
#     Phytoplankton = sum(abundance, na.rm = TRUE)
#   ) %>%
#   dplyr::ungroup() %>%
#   
#   # Existing taxon-level aggregation
#   dplyr::group_by(dplyr::across(
#     -c(abundance, Phytoplankton)
#   )) %>%
#   # 
#   # dplyr::group_by(dplyr::across(
#   #   -abundance
#   #   )
#   #   ) %>% 
#   dplyr::summarise(
#     abundance = sum(abundance),
#     .groups = "drop"
#     ) %>% 
#   dplyr::ungroup() %>% 
#   dplyr::mutate(
#     lb = case_when(
#       stringr::str_starts(taxon, "Phaeocystis") ~ "Phaeocystis",
#       lifeform == "Diatom" ~ "Diatom",
#       lifeform == "Dinoflagellate" ~ "Dinoflagellate",
#       lifeform == "Protozoa" ~ "Protozoa",
#       TRUE ~ NA
#     )
#   ) %>% 
#   dplyr::filter(!(is.na(lb))) %>% 
#   dplyr::select(-c(taxon,lifeform)) %>% 
#   dplyr::group_by(dplyr::across(-abundance)) %>% 
#   dplyr::summarise(
#     abundance = sum(abundance),
#     .groups = "drop"
#   ) %>%
#   ungroup() %>% 
#   tidyr::pivot_wider(
#     names_from = lb,
#     values_from = abundance,values_fill = 0
#     ) %>% 
#   dplyr::filter(river_basin_district != "0") %>% 
#   dplyr::filter(river_basin_district != "Western Wales") %>% 
#   droplevels() -> df_plot

df0 %>%
  tidyr::as_tibble() %>%
  dplyr::filter(data_set == "Phytoplankton") %>%
  dplyr::mutate(
    yday = lubridate::yday(sample_date),
    year = lubridate::year(sample_date),
    month = lubridate::month(sample_date)
  ) %>%
  dplyr::select(
    sample_id, river_basin_district, wb, sample_date,
    yday, year, month, taxon, lifeform, abundance
  ) %>%
  
  # Calculate total abundance of ALL phytoplankton in each sample
  dplyr::group_by(
    sample_id, river_basin_district, wb, sample_date,
    yday, year, month
  ) %>%
  dplyr::mutate(
    Phytoplankton = sum(abundance, na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  
  # Existing taxon-level aggregation
  dplyr::group_by(dplyr::across(
    -c(abundance, Phytoplankton)
  )) %>%
  dplyr::summarise(
    abundance = sum(abundance, na.rm = TRUE),
    Phytoplankton = dplyr::first(Phytoplankton),
    .groups = "drop"
  ) %>%
  
  dplyr::mutate(
    lb = dplyr::case_when(
      stringr::str_starts(taxon, "Phaeocystis") ~ "Phaeocystis",
      lifeform == "Diatom" ~ "Diatom",
      lifeform == "Dinoflagellate" ~ "Dinoflagellate",
      lifeform == "Protozoa" ~ "Protozoa",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(lb)) %>%
  dplyr::select(-c(taxon, lifeform)) %>%
  
  dplyr::group_by(dplyr::across(
    -c(abundance, Phytoplankton)
  )) %>%
  dplyr::summarise(
    abundance = sum(abundance, na.rm = TRUE),
    Phytoplankton = dplyr::first(Phytoplankton),
    .groups = "drop"
  ) %>%
  
  tidyr::pivot_wider(
    names_from = lb,
    values_from = abundance,
    values_fill = 0
  ) %>%
  dplyr::filter(
    river_basin_district != "0",
    river_basin_district != "Western Wales"
  ) %>%
  droplevels() -> df_plot

# generate data & model summaries ----
in_var <- c(
  "Phytoplankton",
  "Diatom",
  "Dinoflagellate",
  "Phaeocystis",
  "Protozoa"
  )

# Safely fit a GAM to one basin-year-taxon dataset
safe_gam <- purrr::safely(
  function(data) {
    mgcv::gam(
      abundance ~ s(yday, bs = "cr", k = 3),
      data = data
      )
    }
  )

# Find the highest observed abundance and its sample date
get_max <- function(data) {
  
  valid_data <- data %>%
    dplyr::filter(
      !is.na(abundance),
      !is.na(sample_date)
    )
  
  if (nrow(valid_data) == 0L) {
    return(
      dplyr::tibble(
        max_date = as.Date(NA),
        max_val = NA_real_
      )
    )
  }
  
  valid_data %>%
    dplyr::slice_max(
      order_by = abundance,
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::transmute(
      max_date = sample_date,
      max_val = abundance
    )
  }

# Create one model per basin, year and taxon
mods <- df_plot %>%
  dplyr::filter(
    year != 2026,
    # year > 2006
  ) %>%
  dplyr::select(
    river_basin_district,
    year,
    sample_date,
    yday,
    dplyr::all_of(in_var)
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(in_var),
    names_to = "taxon",
    values_to = "abundance"
  ) %>%
  droplevels() %>%
  dplyr::group_by(
    river_basin_district,
    year,
    taxon
  ) %>%
  tidyr::nest() %>%
  dplyr::mutate(
    fit = purrr::map(data, safe_gam),
    model = purrr::map(fit, "result"),
    error = purrr::map(fit, "error"),
    maximum = purrr::map(data, get_max)
    ) %>%
  tidyr::unnest(maximum) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    summary = purrr::map(
      model,
      ~ if (is.null(.x)) NULL else summary(.x)
      ),
    classification = purrr::map(
      model,
      ~ if (is.null(.x)) NULL else classify_gam_smooths(.x)
      )
    )

get_gam_peak <- function(model, data) {
  
  # Return NA if GAM failed
  if (is.null(model)) {
    return(
      dplyr::tibble(
        gam_peak_yday = NA_real_,
        gam_peak_date = as.Date(NA),
        gam_peak_abundance = NA_real_
      )
    )
  }
  
  # Range of days actually sampled
  yday_range <- range(data$yday, na.rm = TRUE)
  
  # Fine prediction grid
  newdata <- dplyr::tibble(
    yday = seq(
      from = yday_range[1],
      to   = yday_range[2],
      by   = 0.1
    )
  )
  
  # GAM predictions
  newdata$predicted_abundance <- predict(
    model,
    newdata = newdata,
    type = "response"
  )
  
  # Find predicted maximum
  peak <- newdata %>%
    dplyr::slice_max(
      order_by = predicted_abundance,
      n = 1,
      with_ties = FALSE
    )
  
  # Get year for conversion back to a date
  yr <- unique(data$sample_date |> lubridate::year())
  
  dplyr::tibble(
    gam_peak_yday = peak$yday,
    
    gam_peak_date = as.Date(
      paste0(yr[1], "-01-01")
    ) + round(peak$yday) - 1,
    
    gam_peak_abundance = peak$predicted_abundance
  )
}

mods <- mods %>%
  dplyr::mutate(
    gam_peak = purrr::map2(
      model,
      data,
      get_gam_peak
    )
  ) %>%
  tidyr::unnest(gam_peak)

# mods %>% 
#   dplyr::select(river_basin_district,year,taxon,max_date,max_val) %>% 
#   as.data.frame(.) %>% dplyr::as_tibble() %>% 
#   dplyr::mutate(max_date_use = dplyr::case_when(
#     max_val == 0 ~ NA_Date_,
#     max_val != 0 ~ max_date,
#     TRUE ~ NA
#     )
#     ) %>% 
#   dplyr::select(-max_val) %>% 
#   dplyr::mutate(
#     yday = lubridate::yday(max_date_use)) -> dat

dat <- mods %>%
  dplyr::select(
    river_basin_district,
    year,
    taxon,
    max_date,
    max_val,
    gam_peak_date,
    gam_peak_yday,
    gam_peak_abundance
  ) %>%
  dplyr::filter(
    max_val != 0
  ) %>% 
  dplyr::mutate(season = seas::mkseas(max_date,width = "DJF")) %>% 
  dplyr::mutate(yday = lubridate::yday(max_date)) %>% 
  dplyr::mutate(taxon = dplyr::case_when(
    taxon == "Phytoplankton" ~ "PHYTOPLANKTON",
    taxon != "Phytoplankton" ~ taxon
  )) %>% 
  dplyr::mutate(taxon = factor(taxon, levels = c(
    "PHYTOPLANKTON",
    "Diatom",
    "Dinoflagellate",
    "Phaeocystis",
    "Protozoa"
  )))

# dat %>% 
#   ggplot2::ggplot(.,
#                   aes(
#                     x = yday,
#                     y = taxon,
#                     colour = year
#                   )) +
#   geom_point() +
#   facet_wrap(. ~ river_basin_district) +
#   theme(
#     strip.text = element_text(face = 2),
#     axis.text = element_text(face = 2),
#   )

alpha <- 0.4

dat %>% 
  dplyr::filter(river_basin_district != "Severn") %>% 
  
  ggplot(., aes(x = year,
                y= yday,
                # y= gam_peak_yday,
                )) +
  ## add season rectangles
  # # spring (MAM)
  # geom_rect(xmin = -Inf, xmax = Inf, ymin = 60, ymax = 151, alpha = 0.05,fill = "#B8D8A0")+
  # # summer (JJA)
  # geom_rect(xmin = -Inf, xmax = Inf, ymin = 152, ymax = 243, alpha = 0.05,fill = "#F3B8C8")+
  # # autumn (SON)
  # geom_rect(xmin = -Inf, xmax = Inf, ymin = 244, ymax = 334, alpha = 0.05,fill = "#6D9A93")+
  # #winter 1(D)
  # geom_rect(xmin = -Inf, xmax = Inf, ymin = 335, ymax = 365, alpha = 0.05,fill = "#9CCEE8")+
  # #winter 2(JF)
  # geom_rect(xmin = -Inf, xmax = Inf, ymin = 1, ymax = 59, alpha = 0.05, fill = "#9CCEE8")+
  ## Seasonal background bands
  # spring (MAM)
  ggplot2::annotate(
    "rect",
    xmin = -Inf, xmax = Inf, ymin = 60, ymax = 151,
    alpha = alpha, fill = "#B8D8A0"
    ) +
  # summer (JJA)
  ggplot2::annotate(
    "rect",
    xmin = -Inf, xmax = Inf, ymin = 152, ymax = 243,
    alpha = alpha, fill = "#F3B8C8"
    ) +
  # autumn (SON)
  ggplot2::annotate(
    "rect",
    xmin = -Inf, xmax = Inf, ymin = 244, ymax = 334,
    alpha = alpha, fill = "#6D9A93"
    ) +
  # winter (December)
  ggplot2::annotate(
    "rect",
    xmin = -Inf, xmax = Inf, ymin = 335, ymax = 365,
    alpha = alpha, fill = "#9CCEE8"
    ) +
  # winter (January-February)
  ggplot2::annotate(
    "rect",
    xmin = -Inf, xmax = Inf, ymin = 1, ymax = 59,
    alpha = alpha, fill = "#9CCEE8"
    ) +
  geom_point(
    # aes(colour = season),
    # alpha = 0.6,
    size = 2.5,
    show.legend = FALSE
    ) +
  geom_smooth(
    method = "gam",linewidth = 1.5,
    formula = y ~ s(x, k = 6),
    se = FALSE
  ) +
  facet_grid(taxon ~ river_basin_district) +
  labs(
    x = "Year",
    y = "Day of annual maximum"
  ) +
  theme_bw()+
  labs(
    title = "Day of the year in which the maximum annual abundances of phytoplankton taxa were observed"
  )+
  theme(
    plot.title = element_text(face = 2),
    strip.text = element_text(face = 2),
    axis.title.x = element_blank(),
    axis.text.y = element_text(face = 2),
    axis.text.x = element_text(face = 2, size = 12),
    axis.title = element_text(face = 2),
    )

ggsave(plot = get_last_plot(),
       filename = paste0(output,
                         "202611_Phenology/",
                         "phyto_annual_maximum_timing_by_rbd.png"),
       width = 18, height = 8, units = "in"
       )

### heatmap ----
dat %>% 
  dplyr::filter(year>2006) %>% 
  ggplot(., aes(
    x = year,
    y = forcats::fct_rev(taxon),
    fill = yday
    )
    ) +
  geom_tile() +
  facet_wrap(~ river_basin_district, ncol = 1) +
  scale_fill_viridis_c(
    option = "turbo",
    name = "Day of year"
  ) +
  labs(
    x = "Year",
    y = NULL
  ) +
  theme_bw()+
  theme(
    plot.title = element_text(face = 2),
    strip.text = element_text(face = 2),
    axis.text = element_text(face = 2),
    axis.title = element_text(face = 2),
  )

ggsave(plot = get_last_plot(),
       filename = paste0(output,
                         "202611_Phenology/",
                         "phyto_peak_timing_heatmap_by_rbd.png"),
       width = 18, height = 8, units = "in"
       )

### change relative to baseline ----
dat |>
  dplyr::filter(river_basin_district != "Severn") %>% 
  dplyr::group_by(river_basin_district, taxon) |>
  dplyr::mutate(
    mean_yday = mean(yday, na.rm = TRUE),
    anomaly = yday - mean_yday
  ) |>
  dplyr::ungroup() %>%
  ggplot(., aes(year, anomaly)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_point() +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 6)
  ) +
  facet_grid(taxon ~ river_basin_district) +
  labs(
    title = "Changes in the timing of annual phytoplankton abundance peaks",
    subtitle = "Annual peak timing relative to the long-term mean for each taxon and River Basin District",
    x = "Year",
    y = "Peak timing anomaly (days)",
    caption = paste0("Zero values represent the mean year-day occurrence of annual abundance maxima for each taxon in a given River Basin District\n",
                     "Positive trends suggest bloom peaks occurring later in the year\n",
                     "Negative trends suggest bloom peaks occurring earlier in the year")
  ) +
  theme_bw()+
  theme(
    plot.title = element_text(face = 2,size = 14),
    plot.subtitle = element_text(face = 2,size = 12),
    plot.caption = element_text(face = 2, size = 11),
    strip.text.x = element_text(face = 2, size=12),
    strip.text.y = element_text(face = 2, size=9),
    # strip.text.y.right = element_text(face = 2, size=14),
    axis.text = element_text(face = 2, size = 11),
    axis.title.y = element_text(face = 2,size = 12),
    axis.title.x = element_blank(),
  )

ggsave(plot = get_last_plot(),
       filename = paste0(output,
                         "202611_Phenology/",
                         "phyto_peak_timing_anomalies_by_rbd.png"),
       width = 18, height = 8, units = "in"
       )

### days per decade ----

trends <- dat |>
  filter(!is.na(yday)) |>
  group_by(river_basin_district, taxon) |>
  reframe(
    tidy(lm(yday ~ year, data = pick(everything())))
  ) |>
  filter(term == "year") |>
  mutate(days_per_decade = estimate * 10)

### Combine distribution and temporal trend ----
ggplot(dat, aes(forcats::fct_rev(taxon), yday, colour = river_basin_district)) +
  geom_boxplot() +
  coord_flip() +
  labs(
    x = NULL,
    y = "Day of annual maximum"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = 2,size = 14),
    plot.subtitle = element_text(face = 2,size = 12),
    plot.caption = element_text(face = 2, size = 11),
    strip.text = element_text(face = 2, size=13),
    # strip.text.y.right = element_text(face = 2, size=14),
    axis.text = element_text(face = 2, size = 11),
    axis.title.y = element_text(face = 2,size = 12),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
  )

######
# RECOMMENDED(?) ----
## Prepare the data --------------------------------------------------------

plot_dat <- dat %>%
  filter(
    !is.na(year),
    !is.na(yday),
    !is.na(taxon),
    !is.na(river_basin_district)
  ) %>%
  mutate(
    year = as.integer(year),
    taxon = factor(taxon,levels = c("PHYTOPLANKTON","Diatom","Dinoflagellate","Phaeocystis","Protozoa")),
    river_basin_district = factor(river_basin_district)
  )

# Approximate first day of each month in a non-leap year
month_breaks <- c(1, 32, 60, 91, 121, 152,
                  182, 213, 244, 274, 305, 335)

month_labels <- month.abb

## fig 1 ----
(fig1 <- ggplot(
  plot_dat,
  aes(
    # x = fct_reorder(taxon, yday, .fun = median, na.rm = TRUE),
    x = taxon,
    y = yday
  )
) +
  geom_boxplot(
    aes(fill = taxon),
    width = 0.65,
    alpha = 0.65,
    outlier.shape = NA,
    show.legend = FALSE
  ) +
  geom_jitter(
    # aes(colour = taxon),
    width = 0.15,
    height = 0,
    alpha = 0.45,
    size = 1.3,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ river_basin_district,
    ncol = 3
  ) +
  scale_fill_manual(values = c("Diatom" = "#c1fba4",
                               "Dinoflagellate" = "#ffd6e0",
                               "Phaeocystis" = "#90f1ef",
                               "Protozoa" = "#ffef9f"
                               ))+
  scale_y_continuous(
    breaks = month_breaks,
    labels = month_labels,
    limits = c(1, 366),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    title = "Seasonal timing of annual phytoplankton peaks",
    subtitle = paste(
      "Boxplots show variation among years;",
      "points are individual annual maxima"
    ),
    x = NULL,
    y = "Timing of annual maximum"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = 2),
    plot.subtitle = element_text(face = 2),
    axis.text.x = element_text(
      # angle = 45,
      hjust = 0.5,
      face=2,
      size = 11
    ),
    axis.text.y = element_text(
      face=2
    ),
    axis.title = element_text(face=2),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold"),
    plot.title.position = "plot"
  ))

ggsave(plot = get_last_plot(),
       filename = paste0(output,
                         "202611_Phenology/",
                         "phyto_peak_timing_distribution_by_rbd.png"),
       width = 18, height = 8, units = "in"
       )

(fig1_horizontal <- fig1 +
  coord_flip())


## fig 2 ----
(fig2 <- ggplot(
  plot_dat,
  aes(x = year, y = yday)
) +
  geom_point(
    aes(colour = taxon),
    alpha = 0.65,
    size = 1.7,
    show.legend = FALSE
  ) +
  geom_smooth(
    aes(colour = taxon, fill = taxon),
    method = "gam",
    formula = y ~ s(x, k = 5),
    se = TRUE,
    linewidth = 0.8,
    alpha = 0.15,
    show.legend = FALSE
  ) +
  facet_grid(
    rows = vars(taxon),
    cols = vars(river_basin_district)
  ) +
  scale_x_continuous(
    breaks = breaks_pretty(n = 5)
  ) +
  # scale_colour_manual(values = c("Diatom" = "#c1fba4",
  #                                "Dinoflagellate" = "#ffd6e0",
  #                                "Phaeocystis" = "#90f1ef",
  #                                "Protozoa" = "#ffef9f"
  #                                ))+
  scale_y_continuous(
    breaks = month_breaks,
    labels = month_labels,
    limits = c(1, 366),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    title = "Annual phytoplankton peak timing through time",
    subtitle = paste(
      "Points are observed annual phytoplankton maxima;",
      "lines and ribbons show temporal smooths and 95% confidence intervals"
    ),
    x = "Year",
    y = "Timing of annual maximum"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face=2),
    plot.subtitle = element_text(face=2),
    axis.text = element_text(face=2, size = 11),
    axis.title.y = element_text(face=2,size = 13),
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92"),
    strip.text.x = element_text(
      face=2,
      hjust = 0.5,
      size = 13,
    ),
    strip.text.y = element_text(
      angle = 90*3,
      size = 10,
      face = "bold"
    ),
    plot.title.position = "plot"
  ))

ggsave(plot = get_last_plot(),
       filename = paste0(output,
                         "202611_Phenology/",
                         "phyto_peak_timing_trends_gam_by_rbd.png"),
       width = 18, height = 8, units = "in"
       )

#### fig 2 - linear ----
(fig2_linear <- ggplot(
  plot_dat,
  aes(x = year, y = yday)
  ) +
  geom_point(
    aes(colour = taxon),
    alpha = 0.65,
    size = 1.7,
    show.legend = FALSE
  ) +
  geom_smooth(
    aes(colour = taxon, fill = taxon),
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.8,
    alpha = 0.15,
    show.legend = FALSE
  ) +
  facet_grid(
    rows = vars(taxon),
    cols = vars(river_basin_district)
  ) +
  scale_x_continuous(
    breaks = breaks_pretty(n = 5)
  ) +
  scale_y_continuous(
    breaks = month_breaks,
    labels = month_labels,
    limits = c(1, 366)
  ) +
  labs(
    title = "Annual phytoplankton peak timing through time",
    subtitle = paste(
      "Points are observed annual phytoplankton maxima;",
      "lines and ribbons show linear trends and 95% confidence intervals"
      ),
    x = "Year",
    y = "Timing of annual maximum"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face=2),
    plot.subtitle = element_text(face=2),
    axis.text = element_text(face=2, size = 11),
    axis.title.y = element_text(face=2,size = 13),
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92"),
    strip.text.x = element_text(
      face=2,
      hjust = 0.5,
      size = 13,
    ),
    strip.text.y = element_text(
      angle = 90*3,
      size = 10,
      face = "bold"
    ),
    plot.title.position = "plot"
  ))

ggsave(plot = get_last_plot(),
       filename = paste0(output,
                         "202611_Phenology/",
                         "phyto_peak_timing_trends_linear_by_rbd.png"),
       width = 18, height = 8, units = "in"
       )

## fig 3 - Estimated change in peak timing ----
# Fit a linear time trend to each taxon x RBD group -----------------------
trend_estimates <- plot_dat %>%
  group_by(river_basin_district, taxon) %>%
  filter(n_distinct(year) >= 3) %>%
  mutate(
    year_centred = year - mean(year, na.rm = TRUE)
  ) %>%
  group_modify(
    ~ tidy(
      lm(yday ~ year_centred, data = .x),
      conf.int = TRUE
    )
  ) %>% #View()
  ungroup() %>%
  filter(term == "year_centred") %>%
  transmute(
    river_basin_district,
    taxon,
    estimate_per_year = estimate,
    lower_per_year = conf.low,
    upper_per_year = conf.high,
    days_per_decade = estimate * 10,
    lower_per_decade = conf.low * 10,
    upper_per_decade = conf.high * 10,
    p_value = p.value,
    trend = case_when(
      upper_per_decade < 0 ~ "Earlier",
      lower_per_decade > 0 ~ "Later",
      TRUE ~ "Uncertain"
    )
  )

(fig3 <- ggplot(
  trend_estimates,
  aes(
    x = days_per_decade,
    y = fct_rev(taxon),
    colour = trend
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.6
  ) +
  geom_errorbarh(
    aes(
      xmin = lower_per_decade,
      xmax = upper_per_decade
    ),
    height = 0.18,
    linewidth = 0.65
  ) +
  geom_point(size = 2.6) +
  facet_wrap(
    ~ river_basin_district,
    ncol = 3
  ) +
  scale_colour_manual(
    values = c(
      # "Earlier" = "#2166AC",
      "Earlier" = "blue",
      "Uncertain" = "grey55",
      "Later" = "#B2182B"
    ),
    breaks = c("Earlier", "Uncertain", "Later")
  ) +
  labs(
    title = "Estimated change in phytoplankton peak timing",
    subtitle = paste(
      "Points are linear trends and horizontal bars are 95% confidence intervals;",
      "negative values indicate earlier peaks"
    ),
    x = "Change in peak timing (days per decade)",
    y = NULL,
    colour = "Direction"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(face=2),
    legend.title = element_text(face=2),
    plot.title = element_text(face=2),
    plot.subtitle = element_text(face=2),
    axis.text = element_text(face=2, size = 11),
    axis.title.y = element_text(face=2,size = 13),
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.background = element_rect(fill = "grey92"),
    strip.text.x = element_text(
      face=2,
      hjust = 0.5,
      size = 13,
      ),
    strip.text.y = element_text(
      angle = 90*3,
      size = 12,
      face = "bold"
      ),
    plot.title.position = "plot"
    ))

ggsave(plot = get_last_plot(),
       filename = paste0(output,
                         "202611_Phenology/",
                         "phyto_peak_timing_trends_per_decade_by_rbd.png"),
       width = 18, height = 8, units = "in"
       )
