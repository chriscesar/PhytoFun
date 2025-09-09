# 05_homogeniseUnits.R ####

# load packages ####
ld_pkgs <- c("tidyverse","tictoc","stringr", "forcats", "patchwork")
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
    #x = yday,
    x = sample_date,
    y = log10(md_carb_ugC_per_m3+1),
    fill = data_set,
    shape = data_set
                )) +
  geom_point() +
  scale_shape_manual(values = c(21,24))+
  scale_fill_manual(values=c("green","blue"))+
  geom_smooth(method = "gam",
              formula = y ~ s(x, bs = "ds",k=10),# for date
              # formula = y ~ s(x, bs = "cc",k=10),# for day of year
              se=FALSE, aes(colour = data_set))+
  scale_colour_manual(values=c("green","blue"))+
  facet_wrap(.~regn_wb_biosys_lbl) +
  labs(title = "Estimated carbon content in phytoplankton and zooplankton populations",
       subtitle = "Plots faceted by Region_Waterbody_Monitoring Site",
       # caption = "Lines indicate generalised additive model predictions with cyclic cubic spline basis function",
       caption = "Lines indicate generalised additive model predictions with Duchon splines",
       # x="Day of year",
       #x="Date",
       y=bquote(bold(log[10][(n+1)]~Carbon~content)))+
  theme(
    axis.title.y = element_text(face=2),
    axis.title.x = element_blank(),
    strip.text = element_text(face=2),
    legend.title = element_blank()
    )


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
    x = data_set,
    y = log10(md_carb_ugC_per_m3+1)
  ))+
  geom_boxplot(outliers = FALSE)+
  labs(title = "Estimated carbon content for phytoplankton and zooplankton populations",
       subtitle = "Each point represents the summed estimated carbon contents in phytoplankton and zooplankton samples",
       # caption = "Lines indicate generalised additive model predictions with cyclic cubic spline basis function",
       # caption = "Lines indicate generalised additive model predictions with Duchon splines",
       # x="Day of year",
       x="Date",
       y=bquote(bold(log[10][(n+1)]~Carbon~content)))+
  geom_jitter(width = 0.3, alpha = 0.2)+
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(face=2)
  )

median(log10(dfall_use$md_carb_ugC_per_m3),na.rm=TRUE)

dfall_use %>% 
  group_by(data_set) %>% 
  summarise(log10(median(md_carb_ugC_per_m3+1,na.rm=TRUE)))

## Mean C ####
dfall_use %>% 
  ## generate month for plotting
  dplyr::mutate(month = lubridate::month(sample_date)) %>% 
  dplyr::select(region,wb_id,wb,wb_type,biosys_short,
                sample_date, month, data_set,
                abundance_m3, mn_carb_ugC_per_m3,md_carb_ugC_per_m3) %>% 
  dplyr::filter(sample_date > "2022-06-01") %>% #names()
  dplyr::group_by(across(c(
    -abundance_m3,
    -mn_carb_ugC_per_m3,
    -md_carb_ugC_per_m3))) %>% 
  dplyr::summarise(abundance_m3 = sum(abundance_m3, na.rm = TRUE),
                   mn_carb_ugC_per_m3 = sum(mn_carb_ugC_per_m3, na.rm = TRUE),
                   md_carb_ugC_per_m3 = sum(md_carb_ugC_per_m3, na.rm = TRUE),
                   .groups = "drop") %>% ungroup() %>% 
  dplyr::select(-sample_date) %>% 
  dplyr::group_by(across(c(
    -abundance_m3,
    -mn_carb_ugC_per_m3,
    -md_carb_ugC_per_m3))) %>% 
  # names()
  dplyr::summarise(abundance_m3 = mean(abundance_m3, na.rm = TRUE),
                   mn_carb_ugC_per_m3 = mean(mn_carb_ugC_per_m3, na.rm = TRUE),
                   md_carb_ugC_per_m3 = mean(md_carb_ugC_per_m3, na.rm = TRUE),
                   .groups = "drop") %>% ungroup() %>% 
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
  # dplyr::filter(data_set=="Phytoplankton") %>% 
  ggplot(., aes(x = month,
                y = mn_carb_ugC_per_m3,
                # colour = data_set,
                shape = data_set,
                fill = data_set
  )) +
  geom_point() +
  scale_fill_manual(values=c("green","blue"))+
  scale_shape_manual(values = c(21,24))+
  facet_wrap(.~regn_wb_biosys_lbl, scales = "free_y")

######
## ts vs individual carbon contents: selected taxa ####
n <- 20 # how many rows does the data need to allow plotting
dfall_use %>% 
  dplyr::select(region,wb_id,wb,wb_type,biosys_short,
                sample_date, data_set, lifeform,
                abundance_m3, mn_carb_ugC_per_m3,md_carb_ugC_per_m3) %>% 
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
  dplyr::mutate(
    yday = lubridate::yday(sample_date),
    month = lubridate::month(sample_date)
  ) -> dfall_use_plots

dfall_use_plots %>% 
  dplyr::filter(sample_date > "2022-06-01") %>% #names()
  dplyr::filter(
    lifeform == "Diatom"|lifeform == "Dinoflagellate"|lifeform == "Protozoa"
    ) %>% 
  dplyr::group_by(regn_wb_biosys_lbl) %>% 
  dplyr::filter(n() >= n) %>% ungroup() %>% 
  
  ### SMOOTHED TS #####
  ggplot(., aes(x = yday,
                y = log10(mn_carb_ugC_per_m3+1),
                colour = lifeform,
                shape = lifeform,
                fill = lifeform
                )) +
  geom_point(colour=1) +
  geom_smooth(method="gam",se=FALSE)+
  scale_fill_manual(values=c("green","blue","red"))+
  scale_colour_manual(values=c("green","blue","red"))+
  scale_shape_manual(values = c(21,24, 23))+
  facet_wrap(.~regn_wb_biosys_lbl)

  ### RIDGEPLOT #####
  # ggplot(.,aes(
  #   x = log10(mn_carb_ugC_per_m3+1),
  #   y = as.factor(month),
  #   fill = lifeform
  # ))+
  # ggridges::geom_density_ridges() +
  # scale_fill_manual(values=c("palegreen","lightblue","hotpink1"))+
  # facet_wrap(.~lifeform, ncol=3)
  # # ###

  # ### RADAR PLOT ####
  # ggplot(.,
  #        aes(
  #          x = month,
  #          y = log10(mn_carb_ugC_per_m3+1),
  #          shape = lifeform,
  #          fill = lifeform
  #          )) +
  # geom_jitter(width = 0.2, show.legend = FALSE) +
  # geom_smooth(method = "gam", show.legend = FALSE) +
  # facet_wrap(~ lifeform) +
  # coord_polar(theta = "x", start = 0, direction = 1) +
  # scale_x_continuous(
  #   breaks = 1:12,
  #   labels = 1:12,
  #   expand = c(0, 0)
  # ) +
  # theme_light() +
  # theme(
  #   panel.grid.major.x = element_line(color = "grey70"),
  #   panel.grid.minor.x = element_blank(),
  #   axis.text = element_text(face=2)
  # )

## Time plots #####
dfall_use_plots_2 <- dfall_use_plots %>% 
  dplyr::filter(data_set=="Phytoplankton") %>% 
  # calculate values for Following lifeforms:
  # Diatom, Dinoflagellate, Protozoa, Other
  mutate(
    lifeform_plot = case_when(
      lifeform == "Diatom" ~ "Diatom",
      lifeform == "Dinoflagellate" ~ "Dinoflagellate",
      lifeform == "Protozoa" ~ "Protozoa",
      TRUE ~ "Other"
    )
  )  %>% 
  mutate(lifeform_plot = factor(lifeform_plot,
                                levels = c(
                                  "Diatom",
                                  "Dinoflagellate",
                                  "Protozoa",
                                  "Other"
                                  ))) %>% 
  select(-c(lifeform)) %>% 
  dplyr::group_by(across(c(
    -abundance_m3,
    -mn_carb_ugC_per_m3,
    -md_carb_ugC_per_m3))) %>% 
    summarise(
      abundance_m3 = sum(abundance_m3, na.rm = TRUE),
      mn_carb_ugC_per_m3 = sum(mn_carb_ugC_per_m3,
                                       na.rm = TRUE),
      md_carb_ugC_per_m3=sum(md_carb_ugC_per_m3, na.rm=TRUE),
      .groups = "drop") %>% ungroup()


min_date <- "2007-01-01"
min(dfall_use_plots_2$sample_date)

#####  TOTAL #####
dfall_use_plots_2 %>% 
  mutate(lab="TOTAL") %>% 
  dplyr::filter(
    sample_date >= min_date
  ) %>% 
  group_by(across(c(-lifeform_plot))) %>% 
  summarise(
    abundance_m3 = sum(abundance_m3, na.rm = TRUE),
    mn_carb_ugC_per_m3 = sum(mn_carb_ugC_per_m3,
                             na.rm = TRUE),
    md_carb_ugC_per_m3=sum(md_carb_ugC_per_m3, na.rm=TRUE),
    .groups = "drop") %>% 
  dplyr::filter(
    mn_carb_ugC_per_m3 > 0
  ) %>% 
  ggplot(.,
         aes(
           x = log10(mn_carb_ugC_per_m3+1),
           y = as.factor(month)
         )) +
  facet_wrap(.~lab)+
  ggridges::geom_density_ridges(
    fill = "darkgrey",
    show.legend = FALSE) +
  labs(y = "Month",
       x = bquote(bold(Log[10][(n+1)]~Total~carbon~content~(ug~C~m^-3)))
  )+
  theme(
    strip.text = element_text(face = 2),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = 2),
    panel.background = element_rect(fill = "#EEEEEE")
  ) -> pl1_tot

#####  by LF #####
dfall_use_plots_2 %>% 
  dplyr::filter(
    sample_date >= min_date
  ) %>% 
  dplyr::filter(
    mn_carb_ugC_per_m3 > 0
  ) %>% 
  ggplot(.,
         aes(
           x = log10(md_carb_ugC_per_m3+1),
           y = as.factor(month),
           fill = lifeform_plot
         )) +
  ggridges::geom_density_ridges(show.legend = FALSE)+
  facet_wrap(.~lifeform_plot, ncol = 4)+
  scale_fill_manual(values=c("palegreen","lightblue",
                             "hotpink1","azure2"
                             )) +
  labs(caption = paste0("Phytoplankton data gathered since ",min_date),
       x = bquote(bold(Log[10][(n+1)]~Total~carbon~content~(ug~C~m^-3)))
       )+
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(face = 2)
  ) -> pl1_lfs

pout0 <- pl1_tot+pl1_lfs+plot_layout(widths = c(1, 4),axes = "collect")

ggsave(plot = pout0,
       filename = paste0("outputs/figs/AllRegions.tiff"),
       width = 18, height = 8, units = "in"
)

## By Region ####

filter_to_filename <- function(filter_expr) {
  # Turn expression into a character string, remove spaces & special chars
  expr_text <- rlang::expr_text(filter_expr)
  expr_text |>
    gsub(" ", "", x = _) |>
    gsub("[^A-Za-z0-9_]", "_", x = _)
}

range(dfall_use_plots_2$sample_date)
# dat_use <- list(
#   expr = expr(sample_date >= as.Date("2007-01-01")),
#   label = "From 01/01/2007 onwards",
#   filename = "_since_2007_"
#   )

# dat_use <- list(
#   expr = expr(sample_date >= as.Date("2007-01-01") & sample_date <= as.Date("2012-12-31")),
#   label = "01/01/2007 to 31/12/2012",
#   filename = "_2007_to_2012"
#   )
# 
dat_use <- list(
  expr = expr(sample_date >= as.Date("2020-01-01") & sample_date <= as.Date("2025-03-15")),
  label = "01/01/2020 to 15/03/2025",
  filename = "_2020_to_present_"
)

unique(dfall_use_plots_2$region)
# regname <- "Northumbria"
# regname <- "Humber"
# regname <- "Anglian"
# regname <- "Thames"
# regname <- "South East"
# regname <- "South West"
regname <- "North West"

##### TOTAL ####
dfall_use_plots_2 %>% 
  mutate(lab="TOTAL") %>% 
  dplyr::filter(
    !!dat_use$expr
  ) %>% 
  group_by(across(c(-lifeform_plot))) %>% 
  summarise(
    abundance_m3 = sum(abundance_m3, na.rm = TRUE),
    mn_carb_ugC_per_m3 = sum(mn_carb_ugC_per_m3,
                             na.rm = TRUE),
    md_carb_ugC_per_m3=sum(md_carb_ugC_per_m3, na.rm=TRUE),
    .groups = "drop") %>% 
  dplyr::filter(
    mn_carb_ugC_per_m3 > 0
  ) %>% 
  dplyr::filter(region == regname) %>%
  ggplot(.,
         aes(
           x = log10(mn_carb_ugC_per_m3+1),
           y = as.factor(month)
         )) +
  facet_wrap(.~lab)+
  ggridges::geom_density_ridges(
    fill = "darkgrey",
    show.legend = FALSE) +
  labs(y = "Month",
       x = bquote(bold(Log[10][(n+1)]~Total~carbon~content~(ug~C~m^-3)))
  )+
  ggtitle(regname)+
  theme(
    strip.text = element_text(face = 2),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = 2),
    panel.background = element_rect(fill = "#EEEEEE")
  ) -> pl2_tot_Northumb

#####  by LF #####
dfall_use_plots_2 %>% 
  dplyr::filter(
    !!dat_use$expr
    ) %>% 
  dplyr::filter(
    mn_carb_ugC_per_m3 > 0
  ) %>% 
  dplyr::filter(
    region == regname
  ) %>% 
  ggplot(.,
         aes(
           x = log10(md_carb_ugC_per_m3+1),
           y = as.factor(month),
           fill = lifeform_plot
         )) +
  ggridges::geom_density_ridges(show.legend = FALSE)+
  facet_wrap(.~lifeform_plot, ncol = 4)+
  scale_fill_manual(values=c("palegreen","lightblue",
                             "hotpink1","azure2"
  )) +
  labs(
    caption = paste0("Data window: ", dat_use$label),
    x = bquote(bold(Log[10][(n+1)]~Total~carbon~content~(ug~C~m^-3)))
  )+
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(face = 2)
  ) -> pl2_lfs_Northumb

pout <- (pl2_tot_Northumb+pl2_lfs_Northumb+
           plot_layout(widths = c(1, 4),axes = "collect"))

ggsave(plot = pout,
       filename = paste0("outputs/figs/",regname,"_",dat_use$filename,".tiff"),
       width = 18, height = 8, units = "in"
       )

## Consolidate into single plots ####
## for loop to generate images across all regions ####

# Get all unique region values
regions <- unique(dfall_use_plots_2$region)

# df_region <- dfall_use_plots_2 %>%
#   filter(region == "Anglian")
tic("Generate time shift figures")
for (reg in regions) {
  
  # Filter data for the current region
  df_region <- dfall_use_plots_2 %>%
    filter(region == reg) %>% 
    filter(mn_carb_ugC_per_m3 >0)
  
  # First plot (total)
  pl_tot <- df_region %>%
    mutate(lab = "TOTAL") %>%
    mutate(
      BLOCK = case_when(
        sample_date < ymd("2013-01-01") ~ "Old",
        sample_date > ymd("2020-01-01") ~ "New",
        TRUE ~ "Mid"
      )
    ) %>%
    filter(BLOCK != "Mid") %>%
    ggplot(aes(
      x = log10(mn_carb_ugC_per_m3 + 1),
      y = as.factor(month),
      fill = BLOCK
    )) +
    facet_wrap(. ~ lab) +
    ggridges::geom_density_ridges(
      alpha = 0.5,
      show.legend = FALSE
      ) +
    labs(
      y = "Month",
      x = bquote(bold(Log[10][(n+1)]~Total~carbon~content~(ug~C~m^-3)))
    )+
    ggtitle(reg)+
    theme(
      strip.text = element_text(face = 2),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = 2),
      panel.background = element_rect(fill = "#EEEEEE"),
      legend.title = element_blank()
    )
  
  # Second plot (lifeform)
  pl_lf <- df_region %>%
    mutate(
      BLOCK = case_when(
        sample_date < ymd("2013-01-01") ~ "Old",
        sample_date > ymd("2020-01-01") ~ "New",
        TRUE ~ "Mid"
      )
    ) %>%
    filter(BLOCK != "Mid") %>%
    ggplot(aes(
      x = log10(md_carb_ugC_per_m3 + 1),
      y = as.factor(month),
      fill = BLOCK
    )) +
    ggridges::geom_density_ridges(
      alpha = 0.5,
      # jittered_points = TRUE,
      # position = ggridges::position_points_jitter(width = 0.1, height = -.05),
      # point_shape = '|', point_size = 2, point_alpha = 1, #alpha = 0.7,
      # aes(colour = BLOCK)
      ) +
    facet_wrap(. ~ lifeform_plot, ncol = 4) +
    labs(
      x = bquote(bold(Log[10][(n+1)]~Total~carbon~content~(ug~C~m^-3))),
      caption = "'New' values from phytoplankton samples gathered since 01/01/2020\n'Old' values gathered between 01/01/2007 and 31/12/2022"
    )+
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text = element_text(face = 2),
      legend.title = element_blank()
    )
  
  # Combine plots with patchwork
  pl <- pl_tot + pl_lf +
    patchwork::plot_layout(
      widths = c(1, 4),
      axes = "collect",
      guides = "collect"
    )
  
  # Save to file with region in filename
  out_file <- paste0("outputs/figs/Timeshift_", reg, ".tiff")
  ggsave(out_file, plot = pl, width = 18, height = 8, dpi = 300)
  
  message("Saved plot for region: ", reg)
  }

toc(log=TRUE)
