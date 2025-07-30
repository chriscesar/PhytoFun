# 05_homogeniseUnits.R ####

# load packages ####
ld_pkgs <- c("tidyverse","tictoc","stringr", "forcats")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

source("R/000setup.R")

tic("Load data")
# Load data ####
dfall <- readRDS("outputs/ZoopPhytoMatchingBIOSYS.Rdat") #generated in 03_JoinPhytoZoops.R
toc(log=TRUE)

tic("Convert values to consistent units")
# Convert values to consistent units ####
dfall_use <- dfall %>% 
  # remove taxon values and sum total carbon contents by lifeform
  dplyr::select(-c(
    #aphia_id, taxon,
    kingdom, phylum,
    class, order,
    family, genus
    )) %>% 
  ### convert values to the same units
  ### Convert cells per litre to cells per m3
  mutate(
    abundance_m3 = if_else(abundance_units == "cells_per_litre", abundance / 0.001, abundance)
  ) %>%
  dplyr::relocate(abundance_m3, .after = abundance_units) %>% 
  ### Zoop carbon values are in ugC/m3; phyto values are pgC/l
  ### Convert both to ugC/m3:
  #### MASS Conversion: 1 microgram (μg) = 1,000,000 picograms (pg)
  #### SO: 1 pg = 0.000001 μg
  #### VOLUME CONVERSION:  1 cubic metre (m³) = 1,000 litres (L)
  #### SO: 1 L = 0.001 m3
  #### COMBINED CONVERSION:
  #### pg C/L × (0.000001 μg/pg) × (1,000 L/m³) = μg C/m3
  #### Simplifies to: pg C/L × 0.001 = μg C/m3
  #### OR: μg C/m3 = pg C/L ÷ 1,000
  dplyr::mutate(
    mn_carb_tot_ugC_m3 = case_when(
      mn_carb_tot_units == "ugC_per_m3" ~ mn_carb_tot,
      mn_carb_tot_units == "pg_C_per_l" ~ mn_carb_tot/1000,
      is.na(mn_carb_tot) ~ NA_real_,
      TRUE ~ NA_real_
    )
  ) %>% 
  dplyr::mutate(
    md_carb_tot_ugC_m3 = case_when(
      md_carb_tot_units == "ugC_per_m3" ~ md_carb_tot,
      md_carb_tot_units == "pg_C_per_l" ~ md_carb_tot/1000,
      is.na(md_carb_tot) ~ NA_real_,
      TRUE ~ NA_real_
    )
  ) %>% 
  dplyr::relocate(mn_carb_tot_ugC_m3, .after = mn_carb_tot_units) %>% 
  dplyr::relocate(md_carb_tot_ugC_m3, .after = md_carb_tot_units) %>% 
  # Turn Region into a factor
  dplyr::mutate(region = factor(region,
                                levels = c(
                                  "Northumbria","Humber",
                                  "Anglian","Thames",
                                  "South East","South West",
                                  "North West"
                                )))
toc(log=TRUE)

# Plots ####
## ts vs individual carbon contents ####
dfall_use %>% 
  dplyr::select(-c(
    aphia_id, taxon,
    abundance, abundance_units,
    #abundance_m3,
    mn_carb_tot, mn_carb_tot_units,
    md_carb_tot, md_carb_tot_units,
    mn_size,mn_size_units,md_size,md_size_units)) %>% 
  dplyr::filter(sample_date > "2022-06-01") %>% 
  group_by(across(c(-abundance_m3, -mn_carb_tot_ugC_m3,-md_carb_tot_ugC_m3))) %>% 
  summarise(abundance_m3 = sum(abundance_m3,na.rm = TRUE),
            mn_carb_tot_ugC_m3 = sum(mn_carb_tot_ugC_m3, na.rm = TRUE),
            md_carb_tot_ugC_m3 = sum(md_carb_tot_ugC_m3, na.rm = TRUE),
            .groups = "drop") %>% 
  ggplot(., aes(x = sample_date,
                y = log10(mn_carb_tot_ugC_m3+1),
                # colour = data_set,
                shape = data_set,
                fill = data_set
                )) +
  geom_point() +
  scale_fill_manual(values=c("green","blue"))+
  scale_shape_manual(values = c(21,24))+
  facet_wrap(.~region)

## show 'missing' carbon abundances ####
dfall_use %>% 
  dplyr::select(-c(
    abundance, abundance_units,
    #abundance_m3,
    mn_carb_tot, mn_carb_tot_units,
    md_carb_tot, md_carb_tot_units,
    mn_size,mn_size_units,md_size,md_size_units)) %>% 
  dplyr::filter(sample_date > "2022-06-01") %>% 
  dplyr::mutate(yday = lubridate::yday(sample_date)) %>% 
  # group_by(across(c(-abundance_m3, -mn_carb_tot_ugC_m3,-md_carb_tot_ugC_m3))) %>% 
  # summarise(abundance_m3 = sum(abundance_m3,na.rm = TRUE),
  #           mn_carb_tot_ugC_m3 = sum(mn_carb_tot_ugC_m3, na.rm = TRUE),
  #           md_carb_tot_ugC_m3 = sum(md_carb_tot_ugC_m3, na.rm = TRUE)) %>% 
  ## create flag to ID whether or not the taxon has Carbon estimates
  dplyr::mutate(flag = ifelse(is.na(mn_carb_tot_ugC_m3),"No carbon estimate","Carbon estimate")) %>%
  dplyr::mutate(flag2 = paste0(data_set,"_",flag)) %>% #View(.)
  ggplot(., aes(
    x = sample_date,
    # x = yday,
    y = log10(abundance_m3),
    colour = flag2,
    fill = flag2,
    shape = flag2
  )) +
  # geom_point(pch = 21) +
  geom_jitter(
    # pch=21,
    width = 0.75)+
  scale_fill_manual(values=c("green","white","blue","white"))+
  scale_colour_manual(values = c("darkgreen","darkblue","darkgreen","darkblue"))+
  scale_shape_manual(values = c(21,21,24,24))+
  facet_wrap(.~region) +
  labs(title = "Phytoplankton and zooplankton abundances by region",
       subtitle = "Symbols filled based on whether they do (filled), or do not (white) have carbon content estimates")+
  xlab("Date")+
  ylab("log10(Abundance)")+
  theme(
    axis.title = element_text(face=2),
    strip.text = element_text(face=2),
    legend.title = element_blank(),
    title = element_text(face=2)
  )

## ts vs total carbon contents ####
dfall_use %>% 
  # put Regions in 'clockwise' order
  dplyr::mutate(
    region = forcats::fct_relevel(
      region,
      "Northumbria", "Humber", "Anglian", "Thames",
      "South East","South West","North West")) %>%
  dplyr::mutate(rgn_wb = paste0(region,": ",wb)) %>% 
  dplyr::select(-c(
    lifeform,
    abundance, abundance_units,
    #abundance_m3,
    mn_carb_tot, mn_carb_tot_units,
    md_carb_tot, md_carb_tot_units,
    mn_size,mn_size_units,md_size,md_size_units,
    aphia_id,taxon
    )) %>% 
  dplyr::filter(sample_date > "2022-06-01") %>% #names()
  ## sense check
  # dplyr::filter(biosys_short == "SOL001") %>%
  # dplyr::filter(sample_date == "2024-02-02") %>% #View(.)
  group_by(across(c(-abundance_m3, -mn_carb_tot_ugC_m3,-md_carb_tot_ugC_m3))) %>% 
  summarise(abundance_m3 = sum(abundance_m3,na.rm = TRUE),
            mn_carb_tot_ugC_m3 = sum(mn_carb_tot_ugC_m3,na.rm = TRUE),
            md_carb_tot_ugC_m3 = sum(md_carb_tot_ugC_m3, na.rm = TRUE),
            .groups = "drop") %>% #View()
  dplyr::mutate(yday = lubridate::yday(sample_date)) %>% 
  ## reorder sites by region
  dplyr::mutate(biosys_short = factor(biosys_short,
                                      levels = unique(biosys_short[order(region)])
                                      )
                ) %>% 
  ## create label variable
  dplyr::mutate(regn_wb_lbl = paste0(vegan::make.cepnames(region),
                                     "_",
                                     vegan::make.cepnames(wb))) %>%
  dplyr::mutate(regn_wb_biosys_lbl = paste0(vegan::make.cepnames(region),
                                     "_",
                                     vegan::make.cepnames(wb),
                                     "_",
                                     biosys_short)) %>%
  dplyr::mutate(regn_wb_lbl = factor(regn_wb_lbl,
                                      levels = unique(regn_wb_lbl[order(region)])
  )
  ) %>% 
  dplyr::mutate(regn_wb_biosys_lbl = factor(regn_wb_biosys_lbl,
                                      levels = unique(regn_wb_biosys_lbl[order(region)])
  )
  ) %>% 
  ggplot(., aes(
    # x = yday,
    x = sample_date,
    y = log10(mn_carb_tot_ugC_m3+1),
    fill = data_set,
    shape = data_set
                )) +
  geom_point() +
  scale_shape_manual(values = c(21,24))+
  scale_fill_manual(values=c("green","blue"))+
  geom_smooth(method = "loess", se=FALSE, aes(colour = data_set))+
  scale_colour_manual(values=c("green","blue"))+
  facet_wrap(.~regn_wb_biosys_lbl, scale = "free_y") +
  # ggh4x::facet_wrap2(
  #   ~biosys_short,
  #   strip = ggh4x::strip_themed(
  #     background_x = list(
  #       region == "Northumbria"~element_rect(fill = "lightblue"),
  #       region == "Humber"~element_rect(fill = "orange"),
  #       region == "Anglian"~element_rect(fill = "lightgreen"),
  #       region == "Thames"~element_rect(fill = "pink"),
  #       region == "South East"~element_rect(fill = "sienna"),
  #       region == "South West"~element_rect(fill = "grey"),
  #       region == "North West"~element_rect(fill = "purple"),
  #       TRUE~element_rect(fill="grey90")
  #       )
  #     )
  #   )+
  labs(title = "Estimated carbon content in phytoplankton and zooplankton populations")+
  xlab("Date")+
  ylab("log10(Carbon content +1)")+
  theme(
    axis.title = element_text(face=2),
    strip.text = element_text(face=2),
    legend.title = element_blank()
  )

