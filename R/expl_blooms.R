# expl_blooms.R ----

# Exploration of the timings and co-occurences of blooms in phytoplankton
# taxa of specific ecological interest
## DATE: 13/07/2026
## DATE Last Updated: 
## AUTHOR: Dr Christopher Cesar
## EDITORS:

# load packages ----
ld_pkgs <- c("tidyverse","ggthemes","tictoc","viridis")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log
source("R/000setup.R")

tictoc::tic("Load data"); print("Load data")
df0 <- readRDS("outputs/ZoopPhyto_long_USE.Rdat") %>% 
  # convert date from character to date
  dplyr::mutate(sample_date = lubridate::date(sample_date)) %>% 
  # create year-day variable
  dplyr::mutate(yday = lubridate::yday(sample_date)) %>% 
  # create season variable
  dplyr::mutate(seas = seas::mkseas(sample_date, width = "DJF")) %>% 
  # drop zoops
  dplyr::filter(data_set == "Phytoplankton")
tictoc::toc(log = TRUE)

# PLOTS -----
tictoc::tic("Calculate sum by sample and plot")
print("Calculate sum by sample and plot")
## Calculate sum by sample and plot ----
df0 %>% #names()
  dplyr::select(
    sample_id, sample_date,yday, seas,
    region, wb_id, wb,
    # lifeform, taxon,
    abundance, abundance_units
    ) %>% #names()
    dplyr::group_by(across(-abundance)) %>% 
  dplyr::summarise(abundance = sum(abundance, na.rm = TRUE),
                   .groups = "drop") %>% ungroup() %>% 
  # create 'bloom' flag
  dplyr::mutate(bloom = case_when(
    abundance >= 1000000 ~ "Bloom",
    abundance < 1000000 ~ "Not bloom",
    TRUE ~ NA_character_
  )) %>% 
  dplyr::filter(abundance > 1) %>% 
  # dplyr::filter(sample_date> "2016-01-01") %>% 
  ggplot2::ggplot(.,
                  aes(
                    x = sample_date,
                    y = log10(abundance)
                  )
                  ) +
  ylim(0,NA) +
  # scale_x_date(
  #   date_breaks = "5 years",
  #   date_labels = "%Y",
  #   limits = as.Date(c("2000-01-01",NA))
  #   ) +
  geom_point(aes(colour = bloom),
             alpha = 0.2,
             show.legend = FALSE
             ) +
  scale_colour_manual(values = c("Bloom" = "red",
                                 "Not bloom" = 1)) +
  facet_wrap(. ~ region) +
  # facet_wrap(. ~ wb) +
  geom_smooth(
    method = "gam",
    # se = FALSE
    ) +
  labs(
    title = "Total phytoplankton densities in EA samples",
    subtitle = "Red cells indicate total sample densities in excess of 10^6 cells per litre",
    caption = "Each point represents the total phytoplankton abundance across all taxa within a sample",
    x = "Date",
    y = "log10 abundance (as cells per litre)",
    ) +
  theme(
    panel.grid.major.x = element_line(colour = "grey80",linetype = 2),
    strip.text = element_text(face = 2),
    axis.text =  element_text(face = 2),
    axis.title = element_text(face = 2),
    plot.title = element_text(face = 2),
    plot.subtitle = element_text(face = 2),
    plot.caption = element_text(face = 2),
  )

tictoc::toc(log=TRUE)

tictoc::tic("Calculate and plot max v1")
print("Calculate and plot max v1")
## Maximum abundance v1 ----
df0 %>% #names()
  dplyr::select(
    sample_id, sample_date,yday, seas,
    region, wb_id, wb,
    lifeform, taxon,
    abundance,
  )  %>% 
  dplyr::filter(!is.na(lifeform)) %>% 
  dplyr::filter(lifeform != "Phytoplankton") %>% 
  group_by(sample_id) %>%
  filter(abundance == max(abundance, na.rm = TRUE)) %>%
  ungroup() %>% 
  # create 'bloom' flag
  dplyr::mutate(bloom = case_when(
    abundance >= 1000000 ~ "Bloom",
    abundance < 1000000 ~ "Not bloom",
    TRUE ~ NA_character_
  )) %>% 
  dplyr::filter(abundance > 1) %>% 
  # dplyr::filter(sample_date> "2016-01-01") %>% 
  ggplot2::ggplot(.,
                  aes(
                    x = sample_date,
                    y = log10(abundance)
                  )
  ) +
  ylim(0,NA) +
  geom_point(aes(colour = bloom),
             alpha = 0.2,
             show.legend = FALSE
  ) +
  scale_colour_manual(values = c("Bloom" = "red",
                                 "Not bloom" = 1)) +
  facet_grid(region ~ lifeform) +
  geom_smooth(method = "gam", se = FALSE) +
  labs(
    title = "Densities of most abundant taxon in EA samples",
    subtitle = "Red cells indicate abundances in excess of 10^6 cells per litre",
    caption = "Each point the abundance of the most abundant taxon within a sample",
    x = "Date",
    y = "log10 abundance (as cells per litre)",
  ) +
  theme(
    panel.grid.major.x = element_line(colour = "grey80",linetype = 2),
    strip.text = element_text(face = 2),
    axis.text =  element_text(face = 2),
    axis.title = element_text(face = 2),
    plot.title = element_text(face = 2),
    plot.subtitle = element_text(face = 2),
    plot.caption = element_text(face = 2),
  )

tictoc::toc(log=TRUE)

tictoc::tic("Calculate and plot max v2")
print("Calculate and plot max v2")
## Maximum abundance v2 ----
n <- 20
df0 %>% #names()
  dplyr::select(
    sample_id, sample_date,yday, seas,
    region, wb_id, wb,
    lifeform, taxon,
    abundance,
  )  %>% 
  dplyr::filter(!is.na(lifeform)) %>% 
  dplyr::filter(lifeform != "Phytoplankton") %>% 
  # find max contributor to abundance in a sample
  group_by(sample_id) %>%
  filter(abundance == max(abundance, na.rm = TRUE)) %>%
  ungroup() %>% 
  # create 'bloom' flag
  dplyr::mutate(bloom = case_when(
    abundance >= 1000000 ~ "Bloom",
    abundance < 1000000 ~ "Not bloom",
    TRUE ~ NA_character_
  )) %>% 
  dplyr::filter(abundance > 1) %>% 
  # remove lifeforms occurring fewer than n times
  dplyr::group_by(lifeform) %>%
  dplyr::filter(dplyr::n() >= n) %>%
  dplyr::ungroup() %>%
  # dplyr::filter(sample_date> "2016-01-01") %>% 
  ggplot2::ggplot(.,
                  aes(
                    x = sample_date,
                    y = log10(abundance),
                    # colour = lifeform,
                    # shape = lifeform,
                    # fill = lifeform
                  )
  ) +
  ylim(0,NA) +
  geom_hline(yintercept = log10(10^6), col = 2, lty = 1)+
  geom_point(
    aes(
      colour = lifeform
      ),
    alpha = 0.1,
    # show.legend = FALSE
  ) +
  scale_shape_manual(values = rep(c(15:18),4))+
  scale_fill_viridis(discrete = TRUE)+
  # scale_linetype(values = c(1:3))+
  # scale_colour_manual(values = c("Bloom" = "red",
  #                                "Not bloom" = 1)) +
  # facet_grid(region ~ lifeform) +
  facet_wrap(. ~ region)+
  geom_smooth(aes(
    group = lifeform,
    colour = lifeform,
    # lty = lifeform,
    ),
    method = "gam", se = FALSE) +
  labs(
    title = "Densities of most abundant taxon in EA samples",
    subtitle = "Red cells indicate abundances in excess of 10^6 cells per litre",
    caption = "Each point the abundance of the most abundant taxon within a sample",
    x = "Date",
    y = "log10 abundance (as cells per litre)",
  ) +
  theme(
    panel.grid.major.x = element_line(colour = "grey80",linetype = 2),
    strip.text = element_text(face = 2),
    axis.text =  element_text(face = 2),
    axis.title = element_text(face = 2),
    plot.title = element_text(face = 2),
    plot.subtitle = element_text(face = 2),
    plot.caption = element_text(face = 2),
    legend.title = element_blank(),
  )

tictoc::toc(log=TRUE)

