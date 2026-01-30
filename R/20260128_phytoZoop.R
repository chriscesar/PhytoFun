### 20260128_phytoZoop.R

# load packages ####
ld_pkgs <- c("tidyverse","tictoc","stringr", "forcats", "patchwork","ggh4x")
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

##################################################################
## Subtract Zoops from Phyto within a given month to identify ####
## level of variability in WBs                                ####
##################################################################

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
  ) -> dfuse

## widen by zoops & phyto

# Subtraction ####
## abundance version ####
dfuse %>% 
  select(-c(mn_carb_ugC_per_m3,md_carb_ugC_per_m3)) %>% 
  pivot_wider(.,
              names_from = data_set,
              values_from = abundance_m3
              ) %>% 
  dplyr::mutate(log10_diff = log10(Phytoplankton) - log10(Zooplankton)) %>% 
  dplyr::filter(!is.na(diff)) %>%
  ggplot(.,
         aes(
           # x = yday,
           x = sample_date,
           y = log10_diff,
         )
  ) +
  geom_point(aes(
    colour = region,
    shape = region,
  ),
  size = 3,
  show.legend = FALSE,
  )+
  scale_shape_manual(values = c(15:17,15:17,15))+
  facet_wrap(.~regn_wb_biosys_lbl) +
  # geom_smooth(method="gam", se=FALSE)+
  # geom_smooth(method="loess", se=FALSE)+
  labs(
    y = "Log10 difference",
    x="",
    title = "Difference in phytoplankton and zooplankton abundances",
    subtitle = "Values represent Log10(Phytoplankton abundance)-Log10(Zooplankton abundance)",
  )+
  theme(
    plot.title = element_text(face = 2),
    plot.subtitle = element_text(face = 2),
    axis.title.y = element_text(face = 2),
    strip.text = element_text(face = 2),
    axis.text = element_text(face = 2),
  )

ggsave(plot = last_plot(),
       filename = paste0("outputs/figs/ZP_diff_abund_no_smooth_ts_facet.pdf"),
       width = 18, height = 8, units = "in"
       )

## median carbon version ####
region <- unique(dfuse$region)
# bg_cols <- cbPaletteFill[c(1:4,"grey",6,7)]
bg_cols <- c("#238ACD", "#FFB935", "#23B888", "#B5EAFE",
             "#F07823", "#E894BC","#FFFE77", "#F8C4E1")
bg_cols <- setNames(bg_cols[seq_along(region)], region)

strip_elems <- lapply(region, function(z)
  element_rect(fill = bg_cols[[z]], color = "black", linewidth = 1)
)

dfuse %>% 
  select(-c(abundance_m3,mn_carb_ugC_per_m3)) %>% 
  pivot_wider(.,
              names_from = data_set,
              values_from = md_carb_ugC_per_m3
  ) %>% 
  dplyr::mutate(log10_diff = log10(Phytoplankton) - log10(Zooplankton)) %>% 
  dplyr::filter(!is.na(diff)) %>%
  ggplot(.,
         aes(
           x = sample_date,
           # x = yday,
           y = log10_diff,
           )
         ) +
  geom_point(aes(
    colour = region,
    shape = region,
    ),
    size = 3,
    show.legend = FALSE,
    )+
  scale_shape_manual(values = c(15:17,15:17,15))+
  facet_wrap(.~regn_wb_biosys_lbl) +
  # geom_smooth(method="gam", se=FALSE)+
  # geom_smooth(method="loess", se=FALSE)+
  labs(
    y = "Log10 difference",
    x="",
    title = "Difference in carbon contents between phytoplankton and zooplankton assemblages",
    subtitle = "Values represent Log10(Phytoplankton carbon)-Log10(Zooplankton carbon)",
    )+
  theme(
    plot.title = element_text(face = 2),
    plot.subtitle = element_text(face = 2),
    axis.title.y = element_text(face = 2),
    strip.text = element_text(face = 2),
    axis.text = element_text(face = 2),
  )
ggsave(plot = last_plot(),
       filename = paste0("outputs/figs/ZP_diff_carbon_no_smooth_ts_facet.pdf"),
       width = 18, height = 8, units = "in"
       )

# Division ####
## abundance version ####
dfuse %>% 
  select(-c(mn_carb_ugC_per_m3,md_carb_ugC_per_m3)) %>% 
  pivot_wider(.,
              names_from = data_set,
              values_from = abundance_m3
  ) %>% 
  dplyr::mutate(log10_div = log10(Phytoplankton) / log10(Zooplankton)) %>% 
  dplyr::filter(!is.na(log10_div)) %>%
  ggplot(.,
         aes(
           # x = yday,
           x = sample_date,
           y = log10_div,
         )
  ) +
  geom_point(aes(
    colour = region,
    shape = region,
  ),
  size = 3,
  show.legend = FALSE,
  )+
  scale_shape_manual(values = c(15:17,15:17,15))+
  facet_wrap(.~regn_wb_biosys_lbl) +
  # geom_smooth(method="gam", se=FALSE)+
  # geom_smooth(method="loess", se=FALSE)+
  labs(
    y = "Log10 difference",
    x="",
    title = "Log10 proportion of zooplankton to phytoplankton abundances",
    subtitle = "Values represent Log10(Phytoplankton abundance)/Log10(Zooplankton abundance)",
  )+
  theme(
    plot.title = element_text(face = 2),
    plot.subtitle = element_text(face = 2),
    axis.title.y = element_text(face = 2),
    strip.text = element_text(face = 2),
    axis.text = element_text(face = 2),
  )

ggsave(plot = last_plot(),
       filename = paste0("outputs/figs/ZP_div_abund_no_smooth_ts_facet.pdf"),
       width = 18, height = 8, units = "in"
       )

## median carbon version ####
region <- unique(dfuse$region)
# bg_cols <- cbPaletteFill[c(1:4,"grey",6,7)]
bg_cols <- c("#238ACD", "#FFB935", "#23B888", "#B5EAFE",
             "#F07823", "#E894BC","#FFFE77", "#F8C4E1")
bg_cols <- setNames(bg_cols[seq_along(region)], region)

strip_elems <- lapply(region, function(z)
  element_rect(fill = bg_cols[[z]], color = "black", linewidth = 1)
)

dfuse %>% 
  select(-c(abundance_m3,mn_carb_ugC_per_m3)) %>% 
  pivot_wider(.,
              names_from = data_set,
              values_from = md_carb_ugC_per_m3
  ) %>% 
  dplyr::mutate(log10_diff = log10(Phytoplankton) / log10(Zooplankton)) %>% 
  dplyr::filter(!is.na(diff)) %>%
  ggplot(.,
         aes(
           x = sample_date,
           # x = yday,
           y = log10_diff,
         )
  ) +
  geom_point(aes(
    colour = region,
    shape = region,
  ),
  size = 3,
  show.legend = FALSE,
  )+
  scale_shape_manual(values = c(15:17,15:17,15))+
  facet_wrap(.~regn_wb_biosys_lbl) +
  # geom_smooth(method="gam", se=FALSE)+
  # geom_smooth(method="loess", se=FALSE)+
  labs(
    y = "Log10 difference",
    x="",
    title = "Log10 proportion of zooplankton to phytoplankton carbon content",
    subtitle = "Values represent Log10(Phytoplankton carbon)/Log10(Zooplankton carbon)",
  )+
  theme(
    plot.title = element_text(face = 2),
    plot.subtitle = element_text(face = 2),
    axis.title.y = element_text(face = 2),
    strip.text = element_text(face = 2),
    axis.text = element_text(face = 2),
  )

ggsave(plot = last_plot(),
       filename = paste0("outputs/figs/ZP_div_carbon_no_smooth_ts_facet.pdf"),
       width = 18, height = 8, units = "in"
       )

# Addition ####
## working with TOTAL carbon


## median carbon version ####

region <- unique(dfuse$region)
# bg_cols <- cbPaletteFill[c(1:4,"grey",6,7)]
bg_cols <- c("#238ACD", "#FFB935", "#23B888", "#B5EAFE",
             "#F07823", "#E894BC","#FFFE77", "#F8C4E1")
bg_cols <- setNames(bg_cols[seq_along(region)], region)

strip_elems <- lapply(region, function(z)
  element_rect(fill = bg_cols[[z]], color = "black", linewidth = 1)
)

dfuse %>% 
  select(-c(abundance_m3,mn_carb_ugC_per_m3)) %>% 
  pivot_wider(.,
              names_from = data_set,
              values_from = md_carb_ugC_per_m3
  ) %>% 
  #calculate summed across zoops & phyto
  dplyr::mutate(log10_diff = log10(Phytoplankton)+log10(Zooplankton)) %>% 
  dplyr::mutate(prop_zoops_log10 = log10(Zooplankton)/log10_diff) %>% 
  dplyr::filter(!is.na(log10_diff)) %>%
  ggplot(.,
         aes(
           # x = sample_date,
           x = yday,
           # y = log10_diff,
           y = prop_zoops_log10,
         )
  ) +
  geom_point(aes(
    # colour = region,
    shape = region,
    alpha = prop_zoops_log10,
    colour = prop_zoops_log10,
    size = prop_zoops_log10,
  ),
  # size = 3,
  show.legend = FALSE,
  )+
  scale_shape_manual(values = c(15:17,15:17,15))+
  facet_wrap(.~regn_wb_biosys_lbl) +
  labs(
    y = "Proportion of zooplankton",
    x="",
    title = "Log10 carbon content estimate by Julian day",
    subtitle = "Values represent Log10(Zooplankton carbon)/(Log10(Zooplankton carbon)+Log10(Phytoplankton carbon))\nIncreasing values indicate higher proportions of zooplankton carbon",
    caption =
    "Point shapes indicate EA Region
    Point size and opacity increases with increasing proportion of zooplankton carbon
    Dashed line indicates loess smoother; solid line indicates GAM smoother with a cyclic spline",
    )+
  # geom_smooth(method="gam", formula = y ~ s(x,bs="cp"),
  #             se=FALSE,
  #             colour=2
  # )+
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cc"),
    method.args = list(
      knots = list(x = c(0, 365))
    ),
    se = FALSE,
    colour=2,
  )+

  # geom_smooth(method="gam",
  #             se=FALSE)+
  geom_smooth(method="loess", se=FALSE, colour="darkgreen",linetype = "dashed")+
  theme(
    plot.title = element_text(face = 2),
    plot.caption = element_text(face = 2, size=10),
    plot.subtitle = element_text(face = 2),
    axis.title.y = element_text(face = 2),
    strip.text = element_text(face = 2),
    axis.text = element_text(face = 2),
  )

ggsave(plot = last_plot(),
       filename = paste0("outputs/figs/ZP_sum_carbon_Z_prop_facet.pdf"),
       width = 18, height = 8, units = "in"
       )

dfuse %>% 
  select(-c(abundance_m3,mn_carb_ugC_per_m3)) %>% 
  pivot_wider(.,
              names_from = data_set,
              values_from = md_carb_ugC_per_m3
  ) %>% 
  #calculate summed across zoops & phyto
  dplyr::mutate(log10_diff = log10(Phytoplankton)+log10(Zooplankton)) %>%
  dplyr::mutate(raw_diff = Phytoplankton+Zooplankton) %>% 
  dplyr::mutate(prop_zoops_log10 = log10(Zooplankton)/log10_diff) %>%
  dplyr::mutate(prop_zoops_raw = Zooplankton/raw_diff) %>% 
  dplyr::filter(!is.na(log10_diff)) %>%
  ggplot(.,
         aes(
           # x = sample_date,
           # x = yday,
           # y = log10_diff,
           x = prop_zoops_log10,
           # x = prop_zoops_raw
         )
  ) +
  geom_histogram(colour = 1, fill = "grey60")+
  labs(
    title = "Distribution of zooplankton proportion values",
    subtitle = "Values indicate log10(zooplankton carbon) as a proportion of total log10(plankton carbon)\nHigher values indicate higher proportions of zooplankton",
    x = "Proportion of zooplankton",
    )+
  theme(
    plot.title = element_text(face = 2),
    plot.subtitle = element_text(face = 2),
    axis.title.x = element_text(face = 2),
    axis.title.y = element_blank(),
    strip.text = element_text(face = 2),
    axis.text = element_text(face = 2),
  )
ggsave(plot = last_plot(),
       filename = paste0("outputs/figs/ZP_sum_carbon_Z_prop_hist.pdf"),
       width = 18, height = 8, units = "in"
       )

dfuse %>% 
  select(-c(abundance_m3,mn_carb_ugC_per_m3)) %>% 
  pivot_wider(.,
              names_from = data_set,
              values_from = md_carb_ugC_per_m3
  ) %>% 
  #calculate summed across zoops & phyto
  dplyr::mutate(log10_diff = log10(Phytoplankton)+log10(Zooplankton)) %>%
  dplyr::mutate(raw_diff = Phytoplankton+Zooplankton) %>% 
  dplyr::mutate(prop_zoops_log10 = log10(Zooplankton)/log10_diff) %>%
  dplyr::mutate(prop_zoops_raw = Zooplankton/raw_diff) %>% 
  dplyr::mutate(month = month(sample_date)) %>% 
  dplyr::filter(!is.na(log10_diff)) %>% 
  ggplot(.,aes(
    x=log(Phytoplankton),y=log(Zooplankton),
    colour = as.factor(month)
    ))+
  geom_point()
