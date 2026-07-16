# expl_blooms.R ----

# Exploration of the timings and co-occurences of blooms in phytoplankton
# taxa of specific ecological interest
## DATE: 13/07/2026
## DATE Last Updated: 
## AUTHOR: Dr Christopher Cesar
## EDITORS:

# load packages & data ----
ld_pkgs <- c("tidyverse","ggthemes","tictoc","viridis","ggh4x")
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

# export sample info ----
df0 %>% #names() #%>% 
  dplyr::filter(data_set == "Phytoplankton") %>% 
  dplyr::select(sample_id,
                site_name,
                region,
                wb_id,
                wb,
                eastings,
                northings,
                biosys_short,
                sample_date,
                ) %>% 
  dplyr::distinct() %>% 
  write.csv(., file=paste0(output,"202607_blooms/raw_phyto.csv"),row.names = FALSE)

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

# Maximum taxon ----
## Seasonal peaks v1 ----
tictoc::tic("Seasonal peaks v1")
print("Seasonal peaks v1")
df0 %>%
  dplyr::select(
    sample_id, sample_date,yday, seas,
    region, wb_id, wb,
    lifeform, taxon,
    abundance,
  )  %>% 
  # create month variable
  dplyr::mutate(month = lubridate::month(sample_date,label = TRUE)) %>% 
  dplyr::filter(!is.na(lifeform)) %>% 
  dplyr::filter(lifeform != "Phytoplankton") %>% 
  # keep only Spring (MAM) and late Summer/early Autumn (ASO)
  dplyr::filter(month %in% c("Mar","Apr","May","Aug","Sep","Oct")) %>% 
  dplyr::mutate(seas_use = case_when(
    month %in% c("Mar","Apr","May") ~ "Spr",
    month %in% c("Aug","Sep","Oct") ~ "Sum/Aut",
    TRUE ~ NA_character_
    )
    ) %>% 
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
  ggplot2::ggplot(.,
                  aes(
                    x = sample_date,
                    y = log10(abundance),
                    )
                  ) +
  geom_hline(yintercept = log10(10^6),col="grey",lty=2) +
  ylim(0,NA) +
  
  ## As BLOOM
  geom_point(
    aes(
      colour = bloom
    ),
    alpha = 0.1,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = c("Bloom" = "red",
                                 "Not bloom" = 1)) +
  
  # ## As LF
  # geom_point(
  #   aes(
  #     colour = lifeform
  #   ),
  #   alpha = 0.1,
  #   show.legend = FALSE
  # ) +

  scale_shape_manual(values = rep(c(15:18),4))+
  scale_fill_viridis(discrete = TRUE)+
  ggh4x::facet_wrap2(
    vars(region,seas_use),
    ncol = 6,
    # strip = ggh4x::strip_nested(bleed = FALSE)
    strip = ggh4x::strip_nested(
      bleed = FALSE,
      background_x = list(
        element_rect(fill = "lightblue", colour = "black"), # region
        element_rect(fill = "grey90", colour = "black"), # seas_use
        element_rect(fill = "lightcyan3", colour = "black"), # seas_use
        element_rect(fill = "lightblue", colour = "black"), # region
        element_rect(fill = "grey90", colour = "black"), # seas_use
        element_rect(fill = "lightcyan3", colour = "black"), # seas_use
        element_rect(fill = "lightblue", colour = "black"), # region
        element_rect(fill = "grey90", colour = "black"), # seas_use
        element_rect(fill = "lightcyan3", colour = "black"), # seas_use
        element_rect(fill = "lightgreen", colour = "black"), # seas_use
        element_rect(fill = "lightpink1", colour = "black"), # seas_use
        element_rect(fill = "lightgreen", colour = "black"), # seas_use
        element_rect(fill = "lightpink1", colour = "black"), # seas_use
        element_rect(fill = "lightgreen", colour = "black"), # seas_use
        element_rect(fill = "lightpink1", colour = "black"), # seas_use
        element_rect(fill = "lightgreen", colour = "black"), # seas_use
        element_rect(fill = "lightpink1", colour = "black"), # seas_use
        element_rect(fill = "lightgreen", colour = "black"), # seas_use
        element_rect(fill = "lightpink1", colour = "black"), # seas_use
        element_rect(fill = "lightgreen", colour = "black"), # seas_use
        element_rect(fill = "lightpink1", colour = "black"), # seas_use
        element_rect(fill = "lightgreen", colour = "black"), # seas_use
        element_rect(fill = "lightpink1", colour = "black"), # seas_use
        element_rect(fill = "lightgreen", colour = "black"), # seas_use
        element_rect(fill = "lightpink1", colour = "black"), # seas_use
        element_rect(fill = "lightgreen", colour = "black"), # seas_use
        element_rect(fill = "lightpink1", colour = "black") # seas_use
        )
      )
    ) +
  geom_smooth(
    aes(
      
    ),
  method = "gam", se = FALSE) +
  labs(
    title = "Densities of most abundant taxon in EA samples",
    subtitle = "Red cells indicate abundances in excess of 10^6 cells per litre",
    caption = paste0("Each point represents the abundance of the most abundant taxon within a sample\n",
                     "Plots are faceted by region and whether they were sampled in spring (Mar-May) or late summer/early autumn (Aug-Oct)\n",
                     "Blue lines indicate generalised additive model trends"),
    x = "Date",
    y = "log10 abundance (as cells per litre)",
    ) +
  theme(
    plot.subtitle = element_text(face = 2),
    plot.title = element_text(face=2),
    plot.caption = element_text(face=2,size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2),
    axis.text.x = element_text(face=2,size=12),
    axis.text.y = element_text(face=2,size=9),
    legend.title = element_blank(),
    legend.text = element_text(face=2,size=12),
    strip.text = element_text(face=2),
    strip.background = element_rect(color = 1),
  )
ggsave(plot = get_last_plot(),
       filename = paste0(output,"202607_blooms/spr_summ_max_taxon.png"),
       width = 18, height = 8, units = "in"
       )
tictoc::toc(log=TRUE)

## Seasonal peaks v2 by region ----
tictoc::tic("Seasonal peaks v2 by region ")
print("Seasonal peaks v2 by region ")
for(name in unique(df0$region)) {
  
  message(paste("Generating plot for:", name))

  # subset data
  data_use <- df0 %>% dplyr::filter(region == name) %>% droplevels()
  
  data_use %>% 
    dplyr::select(
    sample_id, sample_date,yday, seas, wb_id, wb,
    lifeform, taxon,
    abundance,
  )  %>% 
  # create month variable
  dplyr::mutate(month = lubridate::month(sample_date,label = TRUE)) %>% 
  dplyr::filter(!is.na(lifeform)) %>% 
  dplyr::filter(lifeform != "Phytoplankton") %>% 
  # # keep only Spring (MAM) and late Summer/early Autumn (ASO)
  # dplyr::filter(month %in% c("Mar","Apr","May","Aug","Sep","Oct")) %>% 
  # dplyr::mutate(seas_use = case_when(
  #   month %in% c("Mar","Apr","May") ~ "Spr",
  #   month %in% c("Aug","Sep","Oct") ~ "Sum/Aut",
  #   TRUE ~ NA_character_
  #   )
  #   ) %>% 
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
  dplyr::filter(month %in% c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct")) %>%
  droplevels() %>% 
  ggplot2::ggplot(.,
                  aes(
                    x = sample_date,
                    y = log10(abundance),
                  )
  ) +
  geom_hline(yintercept = log10(10^6),col="grey",lty=2) +
  ylim(0,NA) +
  
  ## As BLOOM
  geom_point(
    aes(
      colour = bloom
    ),
    alpha = 0.1,
    show.legend = FALSE
  ) +
  facet_wrap(.~month) +
  scale_colour_manual(values = c("Bloom" = "red",
                                 "Not bloom" = 1)) +
  scale_shape_manual(values = rep(c(15:18),4))+
  scale_fill_viridis(discrete = TRUE)+
  ggh4x::facet_wrap2(
    vars(month),
    ncol = 8,
    # strip = ggh4x::strip_nested(bleed = FALSE)
    strip = ggh4x::strip_nested(
      bleed = FALSE,
      background_x = list(
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        
        element_rect(fill = "violetred3", colour = "black"), # JJA
        element_rect(fill = "violetred3", colour = "black"), # JJA
        element_rect(fill = "violetred3", colour = "black"), # JJA
        
        element_rect(fill = "tan1", colour = "black"), # SO
        element_rect(fill = "tan1", colour = "black") # SO
        )
      )
    ) +
  geom_smooth(
    aes(
    ),
    method = "gam", se = FALSE,
    colour = "blue",
    linetype = 1
    ) +
    geom_smooth(
      aes(
      ),
      method = "lm", se = FALSE,
      colour = "#009933",
      linetype = 2
      ) +
  labs(
    title = "Densities of most abundant taxon in EA samples",
    subtitle = paste0(name),
    caption = paste0("Each point represents the abundance of the most abundant taxon within a sample\n",
                     "Red cells indicate abundances in excess of 10^6 cells per litre\n",
                     "Plots are faceted by survey month, with late autumn and winter samples excluded\n",
                     "Blue solid lines indicate generalised additive model trends\n",
                     "Green dashed lines indicate linear model trends"
                     ),
    x = "Date",
    y = "log10 abundance (as cells per litre)",
  ) +
  theme(
    plot.subtitle = element_text(face = 2),
    plot.title = element_text(face=2),
    plot.caption = element_text(face=2,size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2),
    axis.text.x = element_text(face=2,size=12),
    axis.text.y = element_text(face=2,size=9),
    legend.title = element_blank(),
    legend.text = element_text(face=2,size=12),
    strip.text = element_text(face=2),
    strip.background = element_rect(color = 1),
  )
ggsave(plot = get_last_plot(),
       filename = paste0(output,"202607_blooms/monthly_max_",
       vegan::make.cepnames(name),
       "_max_taxon.png"),
       width = 18, height = 8, units = "in"
       )
rm(name, data_use)
}
tictoc::toc(log=TRUE)

## Seasonal peaks v2 by all ----
tictoc::tic("Seasonal peaks v2 all ")
print("Seasonal peaks v2 by all ")

df0 %>% 
  dplyr::select(
      region,sample_id, sample_date,yday, seas, wb_id, wb,
      lifeform, taxon,
      abundance,
    )  %>% 
    # create month variable
    dplyr::mutate(month = lubridate::month(sample_date,label = TRUE)) %>% 
    dplyr::filter(!is.na(lifeform)) %>% 
    dplyr::filter(lifeform != "Phytoplankton") %>% 
    # # keep only Spring (MAM) and late Summer/early Autumn (ASO)
    # dplyr::filter(month %in% c("Mar","Apr","May","Aug","Sep","Oct")) %>% 
    # dplyr::mutate(seas_use = case_when(
    #   month %in% c("Mar","Apr","May") ~ "Spr",
    #   month %in% c("Aug","Sep","Oct") ~ "Sum/Aut",
    #   TRUE ~ NA_character_
    #   )
    #   ) %>% 
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
    dplyr::filter(month %in% c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct")) %>%
    droplevels() %>% 
    ggplot2::ggplot(.,
                    aes(
                      x = sample_date,
                      y = log10(abundance),
                    )
    ) +
    geom_hline(yintercept = log10(10^6),col="grey",lty=2) +
    ylim(0,NA) +
    
    ## As BLOOM
    geom_point(
      aes(
        colour = bloom
      ),
      alpha = 0.1,
      show.legend = FALSE
    ) +
    facet_wrap(.~month+region) +
    scale_colour_manual(values = c("Bloom" = "red",
                                   "Not bloom" = 1)) +
    scale_shape_manual(values = rep(c(15:18),4))+
    scale_fill_viridis(discrete = TRUE)+
    ggh4x::facet_wrap2(
      vars(region,month),
      ncol = 8,
      # strip = ggh4x::strip_nested(bleed = FALSE)
      strip = ggh4x::strip_nested(
        bleed = FALSE,
        background_x = list(
          element_rect(fill = "lightblue", colour = "black"), # region
          element_rect(fill = "lightblue", colour = "black"), # region
          element_rect(fill = "lightblue", colour = "black"), # region
          element_rect(fill = "lightblue", colour = "black"), # region
          element_rect(fill = "lightblue", colour = "black"), # region
          element_rect(fill = "lightblue", colour = "black"), # region
          element_rect(fill = "lightblue", colour = "black"), # region
          element_rect(fill = "lightblue", colour = "black"), # region
          element_rect(fill = "lightblue", colour = "black"), # region
          
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          
          element_rect(fill = "violetred3", colour = "black"), # JJA
          element_rect(fill = "violetred3", colour = "black"), # JJA
          element_rect(fill = "violetred3", colour = "black"), # JJA
          
          element_rect(fill = "tan1", colour = "black"), # SO
          element_rect(fill = "tan1", colour = "black"), # SO
          
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          
          element_rect(fill = "violetred3", colour = "black"), # JJA
          element_rect(fill = "violetred3", colour = "black"), # JJA
          element_rect(fill = "violetred3", colour = "black"), # JJA
          
          element_rect(fill = "tan1", colour = "black"), # SO
          element_rect(fill = "tan1", colour = "black"), # SO
          
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          
          element_rect(fill = "violetred3", colour = "black"), # JJA
          element_rect(fill = "violetred3", colour = "black"), # JJA
          element_rect(fill = "violetred3", colour = "black"), # JJA
          
          element_rect(fill = "tan1", colour = "black"), # SO
          element_rect(fill = "tan1", colour = "black"), # SO
          
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          
          element_rect(fill = "violetred3", colour = "black"), # JJA
          element_rect(fill = "violetred3", colour = "black"), # JJA
          element_rect(fill = "violetred3", colour = "black"), # JJA
          
          element_rect(fill = "tan1", colour = "black"), # SO
          element_rect(fill = "tan1", colour = "black"), # SO
          
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          
          element_rect(fill = "violetred3", colour = "black"), # JJA
          element_rect(fill = "violetred3", colour = "black"), # JJA
          element_rect(fill = "violetred3", colour = "black"), # JJA
          
          element_rect(fill = "tan1", colour = "black"), # SO
          element_rect(fill = "tan1", colour = "black"), # SO
          
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          
          element_rect(fill = "violetred3", colour = "black"), # JJA
          element_rect(fill = "violetred3", colour = "black"), # JJA
          element_rect(fill = "violetred3", colour = "black"), # JJA
          
          element_rect(fill = "tan1", colour = "black"), # SO
          element_rect(fill = "tan1", colour = "black"), # SO
          
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          
          element_rect(fill = "violetred3", colour = "black"), # JJA
          element_rect(fill = "violetred3", colour = "black"), # JJA
          element_rect(fill = "violetred3", colour = "black"), # JJA
          
          element_rect(fill = "tan1", colour = "black"), # SO
          element_rect(fill = "tan1", colour = "black"), # SO
          
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          
          element_rect(fill = "violetred3", colour = "black"), # JJA
          element_rect(fill = "violetred3", colour = "black"), # JJA
          element_rect(fill = "violetred3", colour = "black"), # JJA
          
          element_rect(fill = "tan1", colour = "black"), # SO
          element_rect(fill = "tan1", colour = "black"), # SO
          
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          
          element_rect(fill = "violetred3", colour = "black"), # JJA
          element_rect(fill = "violetred3", colour = "black"), # JJA
          element_rect(fill = "violetred3", colour = "black"), # JJA
          
          element_rect(fill = "tan1", colour = "black"), # SO
          element_rect(fill = "tan1", colour = "black") # SO
        )
      )
    ) +
    geom_smooth(
      aes(
      ),
      method = "gam", se = FALSE,
      colour = "blue",
      linetype = 1
    ) +
    geom_smooth(
      aes(
      ),
      method = "lm", se = FALSE,
      colour = "#009933",
      linetype = 2
    ) +
    labs(
      title = "Densities of most abundant taxon in EA samples",
      # subtitle = paste0(name),
      caption = paste0("Each point represents the abundance of the most abundant taxon within a sample\n",
                       "Red cells indicate abundances in excess of 10^6 cells per litre\n",
                       "Plots are faceted by survey month, with late autumn and winter samples excluded\n",
                       "Blue solid lines indicate generalised additive model trends\n",
                       "Green dashed lines indicate linear model trends"
      ),
      x = "Date",
      y = "log10 abundance (as cells per litre)",
    ) +
    theme(
      plot.subtitle = element_text(face = 2),
      plot.title = element_text(face=2),
      plot.caption = element_text(face=2,size = 12),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face=2),
      axis.text.x = element_text(face=2,size=12),
      axis.text.y = element_text(face=2,size=9),
      legend.title = element_blank(),
      legend.text = element_text(face=2,size=12),
      strip.text = element_text(face=2),
      strip.background = element_rect(color = 1),
    )

ggsave(plot = get_last_plot(),
         filename = paste0(output,"202607_blooms/monthly_max_ALL_REGN_max_taxon.png"),
         width = 18, height = 36, units = "in"
  )

tictoc::toc(log=TRUE)

# Total abundance ----
## Seasonal peaks v1 ----
tictoc::tic("Seasonal peaks v1")
print("Seasonal peaks v1")
df0 %>%
  dplyr::select(
    sample_id, sample_date,yday, seas,
    region, wb_id, wb,
    lifeform, taxon,
    abundance,
  )  %>% 
  # create month variable
  dplyr::mutate(month = lubridate::month(sample_date,label = TRUE)) %>% 
  dplyr::filter(!is.na(lifeform)) %>% 
  dplyr::filter(lifeform != "Phytoplankton") %>% 
  # keep only Spring (MAM) and late Summer/early Autumn (ASO)
  dplyr::filter(month %in% c("Mar","Apr","May","Aug","Sep","Oct")) %>% 
  dplyr::mutate(seas_use = case_when(
    month %in% c("Mar","Apr","May") ~ "Spr",
    month %in% c("Aug","Sep","Oct") ~ "Sum/Aut",
    TRUE ~ NA_character_
  )
  ) %>% 
  # sum abundance in a sample
  dplyr::select(-c(lifeform,taxon)) %>% 
  group_by(across(-abundance)) %>%
  dplyr::summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  ungroup() %>% 
  # create 'bloom' flag
  dplyr::mutate(bloom = case_when(
    abundance >= 1000000 ~ "Bloom",
    abundance < 1000000 ~ "Not bloom",
    TRUE ~ NA_character_
  )) %>% 
  dplyr::filter(abundance > 1) %>%
  ggplot2::ggplot(.,
                  aes(
                    x = sample_date,
                    y = log10(abundance),
                  )
  ) +
  geom_hline(yintercept = log10(10^6),col="grey",lty=2) +
  ylim(0,NA) +
  
  ## As BLOOM
  geom_point(
    aes(
      colour = bloom
    ),
    alpha = 0.1,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = c("Bloom" = "red",
                                 "Not bloom" = 1)) +
  
  # ## As LF
  # geom_point(
  #   aes(
  #     colour = lifeform
  #   ),
  #   alpha = 0.1,
  #   show.legend = FALSE
  # ) +
  
  scale_shape_manual(values = rep(c(15:18),4))+
  scale_fill_viridis(discrete = TRUE)+
  ggh4x::facet_wrap2(
    vars(region,seas_use),
    ncol = 6,
    # strip = ggh4x::strip_nested(bleed = FALSE)
    strip = ggh4x::strip_nested(
      bleed = FALSE,
      background_x = list(
        element_rect(fill = "lightblue", colour = "black"), # region
        element_rect(fill = "grey90", colour = "black"), # seas_use
        element_rect(fill = "lightcyan3", colour = "black"), # seas_use
        element_rect(fill = "lightblue", colour = "black"), # region
        element_rect(fill = "grey90", colour = "black"), # seas_use
        element_rect(fill = "lightcyan3", colour = "black"), # seas_use
        element_rect(fill = "lightblue", colour = "black"), # region
        element_rect(fill = "grey90", colour = "black"), # seas_use
        element_rect(fill = "lightcyan3", colour = "black"), # seas_use
        element_rect(fill = "lightgreen", colour = "black"), # seas_use
        element_rect(fill = "lightpink1", colour = "black"), # seas_use
        element_rect(fill = "lightgreen", colour = "black"), # seas_use
        element_rect(fill = "lightpink1", colour = "black"), # seas_use
        element_rect(fill = "lightgreen", colour = "black"), # seas_use
        element_rect(fill = "lightpink1", colour = "black"), # seas_use
        element_rect(fill = "lightgreen", colour = "black"), # seas_use
        element_rect(fill = "lightpink1", colour = "black"), # seas_use
        element_rect(fill = "lightgreen", colour = "black"), # seas_use
        element_rect(fill = "lightpink1", colour = "black"), # seas_use
        element_rect(fill = "lightgreen", colour = "black"), # seas_use
        element_rect(fill = "lightpink1", colour = "black"), # seas_use
        element_rect(fill = "lightgreen", colour = "black"), # seas_use
        element_rect(fill = "lightpink1", colour = "black"), # seas_use
        element_rect(fill = "lightgreen", colour = "black"), # seas_use
        element_rect(fill = "lightpink1", colour = "black"), # seas_use
        element_rect(fill = "lightgreen", colour = "black"), # seas_use
        element_rect(fill = "lightpink1", colour = "black") # seas_use
      )
    )
  ) +
  geom_smooth(
    aes(
      
    ),
    method = "gam", se = FALSE) +
  labs(
    title = "Phytoplankton densities in EA samples",
    subtitle = "Red cells indicate total phytoplankton abundances in excess of 10^6 cells per litre",
    caption = paste0("Each point represents the total (summed) abundance of phytoplankton within a sample\n",
                     "Plots are faceted by region and whether they were sampled in spring (Mar-May) or late summer/early autumn (Aug-Oct)\n",
                     "Blue lines indicate generalised additive model trends"),
    x = "Date",
    y = "log10 abundance (as cells per litre)",
  ) +
  theme(
    plot.subtitle = element_text(face = 2),
    plot.title = element_text(face=2),
    plot.caption = element_text(face=2,size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2),
    axis.text.x = element_text(face=2,size=12),
    axis.text.y = element_text(face=2,size=9),
    legend.title = element_blank(),
    legend.text = element_text(face=2,size=12),
    strip.text = element_text(face=2),
    strip.background = element_rect(color = 1),
  )
ggsave(plot = get_last_plot(),
       filename = paste0(output,"202607_blooms/spr_summ_total.png"),
       width = 18, height = 8, units = "in"
)
tictoc::toc(log=TRUE)

## Seasonal peaks v2 by all ----
tictoc::tic("Seasonal peaks v2 all ")
print("Seasonal peaks v2 by all ")

df0 %>% 
  dplyr::select(
    region,sample_id, sample_date,yday, seas, wb_id, wb,
    lifeform, taxon,
    abundance,
  )  %>% 
  # create month variable
  dplyr::mutate(month = lubridate::month(sample_date,label = TRUE)) %>% 
  dplyr::filter(!is.na(lifeform)) %>% 
  dplyr::filter(lifeform != "Phytoplankton") %>% 
  # sum abundance in a sample
  dplyr::select(-c(lifeform,taxon)) %>% 
  group_by(across(-abundance)) %>%
  dplyr::summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  ungroup() %>% 
  # create 'bloom' flag
  dplyr::mutate(bloom = case_when(
    abundance >= 1000000 ~ "Bloom",
    abundance < 1000000 ~ "Not bloom",
    TRUE ~ NA_character_
  )) %>% 
  dplyr::filter(abundance > 1) %>%
  dplyr::filter(month %in% c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct")) %>%
  droplevels() %>% 
  ggplot2::ggplot(.,
                  aes(
                    x = sample_date,
                    y = log10(abundance),
                  )
  ) +
  geom_hline(yintercept = log10(10^6),col="grey",lty=2) +
  ylim(0,NA) +
  
  ## As BLOOM
  geom_point(
    aes(
      colour = bloom
    ),
    alpha = 0.1,
    show.legend = FALSE
  ) +
  facet_wrap(.~month+region) +
  scale_colour_manual(values = c("Bloom" = "red",
                                 "Not bloom" = 1)) +
  scale_shape_manual(values = rep(c(15:18),4))+
  scale_fill_viridis(discrete = TRUE)+
  ggh4x::facet_wrap2(
    vars(region,month),
    ncol = 8,
    # strip = ggh4x::strip_nested(bleed = FALSE)
    strip = ggh4x::strip_nested(
      bleed = FALSE,
      background_x = list(
        element_rect(fill = "lightblue", colour = "black"), # region
        element_rect(fill = "lightblue", colour = "black"), # region
        element_rect(fill = "lightblue", colour = "black"), # region
        element_rect(fill = "lightblue", colour = "black"), # region
        element_rect(fill = "lightblue", colour = "black"), # region
        element_rect(fill = "lightblue", colour = "black"), # region
        element_rect(fill = "lightblue", colour = "black"), # region
        element_rect(fill = "lightblue", colour = "black"), # region
        element_rect(fill = "lightblue", colour = "black"), # region
        
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        
        element_rect(fill = "violetred3", colour = "black"), # JJA
        element_rect(fill = "violetred3", colour = "black"), # JJA
        element_rect(fill = "violetred3", colour = "black"), # JJA
        
        element_rect(fill = "tan1", colour = "black"), # SO
        element_rect(fill = "tan1", colour = "black"), # SO
        
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        
        element_rect(fill = "violetred3", colour = "black"), # JJA
        element_rect(fill = "violetred3", colour = "black"), # JJA
        element_rect(fill = "violetred3", colour = "black"), # JJA
        
        element_rect(fill = "tan1", colour = "black"), # SO
        element_rect(fill = "tan1", colour = "black"), # SO
        
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        
        element_rect(fill = "violetred3", colour = "black"), # JJA
        element_rect(fill = "violetred3", colour = "black"), # JJA
        element_rect(fill = "violetred3", colour = "black"), # JJA
        
        element_rect(fill = "tan1", colour = "black"), # SO
        element_rect(fill = "tan1", colour = "black"), # SO
        
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        
        element_rect(fill = "violetred3", colour = "black"), # JJA
        element_rect(fill = "violetred3", colour = "black"), # JJA
        element_rect(fill = "violetred3", colour = "black"), # JJA
        
        element_rect(fill = "tan1", colour = "black"), # SO
        element_rect(fill = "tan1", colour = "black"), # SO
        
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        
        element_rect(fill = "violetred3", colour = "black"), # JJA
        element_rect(fill = "violetred3", colour = "black"), # JJA
        element_rect(fill = "violetred3", colour = "black"), # JJA
        
        element_rect(fill = "tan1", colour = "black"), # SO
        element_rect(fill = "tan1", colour = "black"), # SO
        
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        
        element_rect(fill = "violetred3", colour = "black"), # JJA
        element_rect(fill = "violetred3", colour = "black"), # JJA
        element_rect(fill = "violetred3", colour = "black"), # JJA
        
        element_rect(fill = "tan1", colour = "black"), # SO
        element_rect(fill = "tan1", colour = "black"), # SO
        
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        
        element_rect(fill = "violetred3", colour = "black"), # JJA
        element_rect(fill = "violetred3", colour = "black"), # JJA
        element_rect(fill = "violetred3", colour = "black"), # JJA
        
        element_rect(fill = "tan1", colour = "black"), # SO
        element_rect(fill = "tan1", colour = "black"), # SO
        
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        
        element_rect(fill = "violetred3", colour = "black"), # JJA
        element_rect(fill = "violetred3", colour = "black"), # JJA
        element_rect(fill = "violetred3", colour = "black"), # JJA
        
        element_rect(fill = "tan1", colour = "black"), # SO
        element_rect(fill = "tan1", colour = "black"), # SO
        
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        element_rect(fill = "lightgreen", colour = "black"), # MAM
        
        element_rect(fill = "violetred3", colour = "black"), # JJA
        element_rect(fill = "violetred3", colour = "black"), # JJA
        element_rect(fill = "violetred3", colour = "black"), # JJA
        
        element_rect(fill = "tan1", colour = "black"), # SO
        element_rect(fill = "tan1", colour = "black") # SO
      )
    )
  ) +
  geom_smooth(
    aes(
    ),
    method = "gam", se = FALSE,
    colour = "blue",
    linetype = 1
  ) +
  geom_smooth(
    aes(
    ),
    method = "lm", se = FALSE,
    colour = "#009933",
    linetype = 2
  ) +
  labs(
    title = "Phytoplankton densities in EA samples",
    subtitle = "Red cells indicate total phytoplankton abundances in excess of 10^6 cells per litre",
    caption = paste0("Each point represents the total (summed) abundance of phytoplankton within a sample\n",
                     "Plots are faceted by region and whether they were sampled in spring (Mar-May) or late summer/early autumn (Aug-Oct)\n",
                     "Blue lines indicate generalised additive model trends"),
    x = "Date",
    y = "log10 abundance (as cells per litre)",
  ) +
  theme(
    plot.subtitle = element_text(face = 2),
    plot.title = element_text(face=2),
    plot.caption = element_text(face=2,size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2),
    axis.text.x = element_text(face=2,size=12),
    axis.text.y = element_text(face=2,size=9),
    legend.title = element_blank(),
    legend.text = element_text(face=2,size=12),
    strip.text = element_text(face=2),
    strip.background = element_rect(color = 1),
  )

ggsave(plot = get_last_plot(),
       filename = paste0(output,"202607_blooms/monthly_ALL_REGN_total.png"),
       width = 18, height = 36, units = "in"
)

tictoc::toc(log=TRUE)

## Seasonal peaks v2 by region ----
tictoc::tic("Seasonal peaks v2 by region ")
print("Seasonal peaks v2 by region ")
for(name in unique(df0$region)) {
  
  message(paste("Generating plot for:", name))
  
  # subset data
  data_use <- df0 %>% dplyr::filter(region == name) %>% droplevels()
  
  data_use %>% 
    dplyr::select(
      sample_id, sample_date,yday, seas, wb_id, wb,
      lifeform, taxon,
      abundance,
    )  %>% 
    # create month variable
    dplyr::mutate(month = lubridate::month(sample_date,label = TRUE)) %>% 
    dplyr::filter(!is.na(lifeform)) %>% 
    dplyr::filter(lifeform != "Phytoplankton") %>%
    # sum abundance in a sample
    dplyr::select(-c(lifeform,taxon)) %>% 
    group_by(across(-abundance)) %>%
    dplyr::summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    ungroup() %>% 
    # create 'bloom' flag
    dplyr::mutate(bloom = case_when(
      abundance >= 1000000 ~ "Bloom",
      abundance < 1000000 ~ "Not bloom",
      TRUE ~ NA_character_
    )) %>% 
    dplyr::filter(abundance > 1) %>%
    dplyr::filter(month %in% c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct")) %>%
    droplevels() %>% 
    ggplot2::ggplot(.,
                    aes(
                      x = sample_date,
                      y = log10(abundance),
                    )
    ) +
    geom_hline(yintercept = log10(10^6),col="grey",lty=2) +
    ylim(0,NA) +
    
    ## As BLOOM
    geom_point(
      aes(
        colour = bloom
      ),
      alpha = 0.1,
      show.legend = FALSE
    ) +
    facet_wrap(.~month) +
    scale_colour_manual(values = c("Bloom" = "red",
                                   "Not bloom" = 1)) +
    scale_shape_manual(values = rep(c(15:18),4))+
    scale_fill_viridis(discrete = TRUE)+
    ggh4x::facet_wrap2(
      vars(month),
      ncol = 8,
      # strip = ggh4x::strip_nested(bleed = FALSE)
      strip = ggh4x::strip_nested(
        bleed = FALSE,
        background_x = list(
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          
          element_rect(fill = "violetred3", colour = "black"), # JJA
          element_rect(fill = "violetred3", colour = "black"), # JJA
          element_rect(fill = "violetred3", colour = "black"), # JJA
          
          element_rect(fill = "tan1", colour = "black"), # SO
          element_rect(fill = "tan1", colour = "black") # SO
        )
      )
    ) +
    geom_smooth(
      aes(
      ),
      method = "gam", se = FALSE,
      colour = "blue",
      linetype = 1
    ) +
    geom_smooth(
      aes(
      ),
      method = "lm", se = FALSE,
      colour = "#009933",
      linetype = 2
    ) +
    labs(
      title = "Phytoplankton densities in EA samples",
      subtitle = paste0(name),
      caption = paste0("Each point represents the total (summed) abundance of phytoplankton within a sample\n",
                       "Plots are faceted by region and whether they were sampled in spring (Mar-May) or late summer/early autumn (Aug-Oct)\n",
                       "Blue lines indicate generalised additive model trends\n",
                       "Green dashed lines indicate linear model trends"
                       ),
      x = "Date",
      y = "log10 abundance (as cells per litre)",
    ) +
    theme(
      plot.subtitle = element_text(face = 2),
      plot.title = element_text(face=2),
      plot.caption = element_text(face=2,size = 12),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face=2),
      axis.text.x = element_text(face=2,size=12),
      axis.text.y = element_text(face=2,size=9),
      legend.title = element_blank(),
      legend.text = element_text(face=2,size=12),
      strip.text = element_text(face=2),
      strip.background = element_rect(color = 1),
    )
  ggsave(plot = get_last_plot(),
         filename = paste0(output,"202607_blooms/monthly_total_",
                           vegan::make.cepnames(name),
                           "_total.png"),
         width = 18, height = 8, units = "in"
  )
  rm(name, data_use)
  }
tictoc::toc(log=TRUE)

## Seasonal peaks v2 by region - to Nov ----
tictoc::tic("Seasonal peaks v2 by region - to Nov ")
print("Seasonal peaks v2 by region - to Nov")
for(name in unique(df0$region)) {
  
  message(paste("Generating plot for:", name))
  
  # subset data
  data_use <- df0 %>% dplyr::filter(region == name) %>% droplevels()
  
  data_use %>% 
    dplyr::select(
      sample_id, sample_date,yday, seas, wb_id, wb,
      lifeform, taxon,
      abundance,
    )  %>% 
    # create month variable
    dplyr::mutate(month = lubridate::month(sample_date,label = TRUE)) %>% 
    dplyr::filter(!is.na(lifeform)) %>% 
    dplyr::filter(lifeform != "Phytoplankton") %>%
    # sum abundance in a sample
    dplyr::select(-c(lifeform,taxon)) %>% 
    group_by(across(-abundance)) %>%
    dplyr::summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    ungroup() %>% 
    # create 'bloom' flag
    dplyr::mutate(bloom = case_when(
      abundance >= 1000000 ~ "Bloom",
      abundance < 1000000 ~ "Not bloom",
      TRUE ~ NA_character_
    )) %>% 
    dplyr::filter(abundance > 1) %>%
    dplyr::filter(month %in% c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov")) %>%
    droplevels() %>% 
    ggplot2::ggplot(.,
                    aes(
                      x = sample_date,
                      y = log10(abundance),
                    )
    ) +
    geom_hline(yintercept = log10(10^6),col="grey",lty=2) +
    ylim(0,NA) +
    
    ## As BLOOM
    geom_point(
      aes(
        colour = bloom
      ),
      alpha = 0.1,
      show.legend = FALSE
    ) +
    facet_wrap(.~month) +
    scale_colour_manual(values = c("Bloom" = "red",
                                   "Not bloom" = 1)) +
    scale_shape_manual(values = rep(c(15:18),4))+
    scale_fill_viridis(discrete = TRUE)+
    ggh4x::facet_wrap2(
      vars(month),
      ncol = 9,
      # strip = ggh4x::strip_nested(bleed = FALSE)
      strip = ggh4x::strip_nested(
        bleed = FALSE,
        background_x = list(
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          element_rect(fill = "lightgreen", colour = "black"), # MAM
          
          element_rect(fill = "violetred3", colour = "black"), # JJA
          element_rect(fill = "violetred3", colour = "black"), # JJA
          element_rect(fill = "violetred3", colour = "black"), # JJA
          
          element_rect(fill = "tan1", colour = "black"), # SO
          element_rect(fill = "tan1", colour = "black"), # SO
          element_rect(fill = "tan1", colour = "black") # SO
        )
      )
    ) +
    geom_smooth(
      aes(
      ),
      method = "gam", se = FALSE,
      colour = "blue",
      linetype = 1
    ) +
    geom_smooth(
      aes(
      ),
      method = "lm", se = FALSE,
      colour = "#009933",
      linetype = 2
    ) +
    labs(
      title = "Phytoplankton densities in EA samples",
      subtitle = paste0(name),
      caption = paste0("Each point represents the total (summed) abundance of phytoplankton within a sample\n",
                       "Plots are faceted by region and whether they were sampled in spring (Mar-May) or late summer/early autumn (Aug-Oct)\n",
                       "Blue lines indicate generalised additive model trends\n",
                       "Green dashed lines indicate linear model trends"
      ),
      x = "Date",
      y = "log10 abundance (as cells per litre)",
    ) +
    theme(
      plot.subtitle = element_text(face = 2),
      plot.title = element_text(face=2),
      plot.caption = element_text(face=2,size = 12),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face=2),
      axis.text.x = element_text(face=2,size=12),
      axis.text.y = element_text(face=2,size=9),
      legend.title = element_blank(),
      legend.text = element_text(face=2,size=12),
      strip.text = element_text(face=2),
      strip.background = element_rect(color = 1),
    )
  ggsave(plot = get_last_plot(),
         filename = paste0(output,"202607_blooms/monthly_total_",
                           vegan::make.cepnames(name),
                           "_total_to_Nov.png"),
         width = 18, height = 8, units = "in"
  )
  rm(name, data_use)
}
tictoc::toc(log=TRUE)
