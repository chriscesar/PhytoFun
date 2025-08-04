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

dfall_use <- dfall %>% 
  # remove taxon values and sum total carbon contents by lifeform
  dplyr::select(-c(
    #aphia_id, taxon,
    kingdom, phylum,
    class, order,
    family, genus
  ))
toc(log=TRUE)

# Plots ####
## ts vs individual carbon contents, grouped by lifeforms ####
dfall_use %>% 
  dplyr::select(region,wb_id,wb,wb_type,biosys_short,
                sample_date, data_set, lifeform,
                abundance_m3, mn_carb_ugC_per_m3,md_carb_ugC_per_m3) %>% 
  dplyr::filter(sample_date > "2022-06-01") %>% #names()
  dplyr::group_by(across(c(
    -abundance_m3,
    -mn_carb_ugC_per_m3,
    -md_carb_ugC_per_m3))) %>% 
  dplyr::summarise(abundance_m3 = sum(abundance_m3, na.rm = TRUE),
                   mn_carb_ugC_per_m3 = sum(mn_carb_ugC_per_m3, na.rm = TRUE),
                   md_carb_ugC_per_m3 = sum(md_carb_ugC_per_m3, na.rm = TRUE),
                   .groups = "drop") %>% 
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
  ggplot(., aes(x = sample_date,
                y = log10(mn_carb_ugC_per_m3+1),
                # colour = data_set,
                shape = data_set,
                fill = data_set
                )) +
  geom_point() +
  scale_fill_manual(values=c("green","blue"))+
  scale_shape_manual(values = c(21,24))+
  facet_wrap(.~regn_wb_biosys_lbl)

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
  dplyr::mutate(flag = ifelse(is.na(mn_carb_ind_as_ugC),"No carbon estimate","Carbon estimate")) %>%
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
  dplyr::mutate(rgn_wb = paste0(region,": ",wb)) %>% 
  dplyr::select(-c(
    lifeform,
    abundance, abundance_units,
    #abundance_m3,
    mn_carb_tot, mn_carb_tot_units,
    md_carb_tot, md_carb_tot_units,
    mn_carb_ind,mn_carb_ind_units,md_carb_ind,
    md_carb_ind_units,mn_carb_ind_as_ugC,md_carb_ind_as_ugC,
    mn_size,mn_size_units,md_size,md_size_units,
    aphia_id,taxon
    )) %>% #names()
  dplyr::filter(sample_date > "2022-06-01") %>% #names()
  group_by(across(c(-abundance_m3,
                    -mn_carb_ugC_per_m3,
                    -md_carb_ugC_per_m3
                    ))) %>% 
  summarise(abundance_m3 = sum(abundance_m3,na.rm = TRUE),
            mn_carb_ugC_per_m3 = sum(mn_carb_ugC_per_m3,na.rm = TRUE),
            md_carb_ugC_per_m3 = sum(md_carb_ugC_per_m3, na.rm = TRUE),
            .groups = "drop") %>% #View()
  dplyr::mutate(yday = lubridate::yday(sample_date)) %>% 
  ## reorder sites by region
  dplyr::mutate(biosys_short = factor(
    biosys_short,
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
  dplyr::mutate(regn_wb_lbl = factor(
    regn_wb_lbl,
    levels = unique(regn_wb_lbl[order(region)])
    )
    ) %>% 
  dplyr::mutate(regn_wb_biosys_lbl = factor(
    regn_wb_biosys_lbl,
    levels = unique(regn_wb_biosys_lbl[order(region)])
    )
    ) %>% 
  ggplot(., aes(
    x = yday,
    #x = sample_date,
    y = log10(md_carb_ugC_per_m3+1),
    fill = data_set,
    shape = data_set
                )) +
  geom_point() +
  scale_shape_manual(values = c(21,24))+
  scale_fill_manual(values=c("green","blue"))+
  geom_smooth(method = "gam",
              # formula = y ~ s(x, bs = "ds",k=10),# for date
              formula = y ~ s(x, bs = "cc",k=10),# for day of year
              se=FALSE, aes(colour = data_set))+
  scale_colour_manual(values=c("green","blue"))+
  facet_wrap(.~regn_wb_biosys_lbl) +
  labs(title = "Estimated carbon content in phytoplankton and zooplankton populations")+
  xlab("Day of year")+
  ylab("log10(Carbon content +1)")+
  theme(
    axis.title = element_text(face=2),
    strip.text = element_text(face=2),
    legend.title = element_blank()
    )

