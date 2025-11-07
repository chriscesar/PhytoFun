# phytoWordcloud.R ####
#### Generate wordclouds of life forms and taxon data

ld_pkgs <- c("tidyverse","ggwordcloud", "tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

ppi <- 300
tictoc::tic.clearlog() ##clear log

#### set universals ####
tic("set meta")
#source("R/folder.links.R") ## data folders
theme_set(ggthemes::theme_few())###set theme for all ggplot objects
ppi <- 300 #image resolution
# colourblind friendly colour palette (RGB values also commented)
cbPalette <- c("#999999", #153/153/153
               "#E69F00",#230/159/000
               "#56B4E9",#086/180/233
               "#CC79A7", #204/121/167
               "#009E73",#000/158/115
               "#F0E442",#240/228/066
               "#0072B2",#000/114/178
               "#D55E00",#213/094/000
               
               "#444444", 
               "#C34D55",
               "#33A2C4",
               "#554C31",
               "#C5C221",
               "#5531A1",
               "#B32C55",
               "#BB3593" 
               
)

cbPalette2 <- c("#646464", #100/100/100
                "#B46D00",#180/109/0
                "#2482BA",#036/130/186
                "#006C41",#000/108/065
                "#BEB210",#190/178/016
                "#004080",#000/064/128
                "#A32C00",#163/044/000
                "#9A4775"#154/071/117
                )
toc(log=TRUE)

tic("load & format data")
# load data ####

df0 <- readRDS("outputs/Phyto_2000_2025_with_Lifeforms_USE_updateCalcs.Rdat")

df <- df0 %>%
  select(river_basi,wb_name,biosys_code,site_id,sample_date,
         # name_use,
         taxa_reported,
         md_carbon_tot_ug_C_per_litre,
         plankton_type:phyto_lf
         ) %>%
  dplyr::mutate(smpref = paste0(site_id,sample_date)) %>% 
  ### drop pre 2011 data
  dplyr::filter(sample_date > "2009-12-31") %>%
  # define colours
  dplyr::mutate(phyto_lf_col = case_when(
    phyto_lf == "Diatom" ~ "chartreuse",
    phyto_lf == "Dinoflagellate" ~ "dodgerblue",
    phyto_lf == "Protozoa" ~ "blueviolet",
    phyto_lf == "Haptophyte" ~ "mediumblue",
    phyto_lf == "Cyanobacteria" ~ "cyan3",
    phyto_lf == "Chrysophyte" ~ "sienna",
    phyto_lf == "Chlorophyte" ~ "aquamarine2",
    phyto_lf == "Silicoflagellate" ~ "deeppink2",
    phyto_lf == "Charophyte" ~ "darkgreen",
    phyto_lf == "Ciliate" ~ "darkorange",
    phyto_lf == "Raphidophyte" ~ "darkred",
    TRUE ~ "darkgrey"
  )) %>% 
  mutate(phyto_lf_col = factor(phyto_lf_col,
                               levels = c("chartreuse","dodgerblue","blueviolet",
                                          "mediumblue","cyan3","sienna","aquamarine2",
                                          "deeppink2","darkgreen","darkorange",
                                          "darkred","darkgrey")
    
  )) %>% 
  dplyr::mutate(phyto_lf = ifelse(is.na(phyto_lf), "Unknown", phyto_lf))
  
          
pal <- c("chartreuse","dodgerblue","blueviolet","mediumblue","cyan3","sienna",
         "aquamarine2","deeppink2","darkgreen","darkorange","darkred","darkgrey")


phyto_col <- c("Diatom" = "chartreuse",
               "Dinoflagellate" = "dodgerblue",
               "Protozoa" = "blueviolet",
               "Haptophyte" = "mediumblue",
               "Cyanobacteria" = "cyan3",
               "Chrysophyte" = "sienna",
               "Chlorophyte" = "aquamarine2",
               "Silicoflagellate" = "deeppink2",
               "Charophyte" = "darkgreen",
               "Ciliate" = "darkorange",
               "Raphidophyte" = "darkred",
               "Unknown" = "darkgrey")
toc(log=TRUE)

##########################
### Mean carbon v2 ####
tic("Generate wordclouds: Mean carbon per taxon v.2")
# all_years <- 2011:2025  # ✅ Updated to full range
all_years <- 2010:2025  # ✅ Updated to full range

set.seed(3123)
for(name in unique(df$river_basi)) {
  
  message(paste("Generating wordcloud for:", name))
  
  # Sanitize filename
  safe_name <- make.names(name)
  
  # Create full grid of year × taxa × phyto_lf combinations
  taxa_phyto <- df %>%
    dplyr::filter(river_basi == name) %>%
    dplyr::distinct(taxa_reported, phyto_lf)
  
  full_grid <- tidyr::expand_grid(
    year = all_years,
    taxa_phyto
  )
  
  site_counts <- df %>%
    dplyr::filter(river_basi == name) %>%
    dplyr::mutate(year = lubridate::year(sample_date)) %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(n_sites = n_distinct(smpref), .groups = "drop")
  
  # Replace original summarisation and complete() with join to full grid
  df %>%
    dplyr::filter(river_basi == name) %>%
    dplyr::mutate(year = lubridate::year(sample_date)) %>%
    dplyr::group_by(year) %>% 
    dplyr::mutate(num = length(unique(smpref))) %>% ungroup() %>% 
    dplyr::filter(!is.na(phyto_lf)) %>%
    dplyr::group_by(year, taxa_reported, phyto_lf) %>%
    dplyr::summarise(mean = mean(md_carbon_tot_ug_C_per_litre, na.rm = TRUE), .groups = "drop") %>%
    dplyr::ungroup() %>%
    
    dplyr::filter(!is.nan(mean)) %>%
    # tidyr::complete(year = all_years, taxa_reported, phyto_lf, fill = list(mean = 0)) %>%  # Commented out
    # Join with full grid to ensure all year-taxonomy combinations are present
    dplyr::right_join(full_grid, by = c("year", "taxa_reported", "phyto_lf")) %>%
    dplyr::mutate(mean = tidyr::replace_na(mean, 0)) %>%
    dplyr::group_by(year) %>%
    dplyr::mutate(tot = sum(mean, na.rm = TRUE)) %>%
    dplyr::mutate(prop = ifelse(tot > 0, mean / tot, 0)) %>%
    dplyr::ungroup() %>%
    # dplyr::filter(prop > 0.001) %>%  # Filter out very small proportions
    
    # Create custom facet labels
    # dplyr::group_by(year) %>%
    # dplyr::mutate(
    #   facet_label = paste0(
    #     year, "\nTotal C: ", round(tot[1], 1), "ug C per litre\nn = ", distin
    #   )
    # ) %>%
    # ungroup() ->x#%>%
    
    dplyr::left_join(site_counts, by = "year") %>%
    dplyr::group_by(year) %>%
    # dplyr::mutate(
    #   facet_label = paste0(
    #     year, "\nTotal C: ", round(tot[1], 1), " ug C per litre\nn = ", n_sites[1]
    #   )
    # ) %>%
    dplyr::mutate(
      facet_label = paste0(
        year,"n = ",n_sites[1],")", "\nMean: ", round(tot[1], 1), " ug C/l"
        )
    ) %>%
    ungroup() %>% 
  
    ggplot(aes(
      label = taxa_reported,
      size = prop,
      colour = phyto_lf
    )) +
    # ggwordcloud::geom_text_wordcloud_area(rm_outside = TRUE, show.legend = TRUE) +
    ggwordcloud::geom_text_wordcloud(rm_outside = TRUE,
                                     show.legend = TRUE,
                                     scale.size = TRUE) +  # Improved layout
    ggplot2::scale_size_area(
      max_size = 50,
      trans = power_trans(1 / .7)
    ) +
    theme_minimal() +
    labs(
      title = name,
      caption = "\nText size reflects proportion of mean carbon content provided by that taxon"
    ) +
    # facet_wrap(. ~ year, scales = "free") +
    facet_wrap(. ~ facet_label, scales = "free") +
    scale_colour_manual(values = phyto_col) +
    guides(
      size = "none",
      colour = guide_legend(nrow = 2, byrow = TRUE)
    ) +
    theme(
      strip.text = element_text(face = "bold", size = 14, colour = "red"),
      title = element_text(face = "bold", size = 20),
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.title = element_blank()
    ) -> pl
  
  # Save PDF
  tic(paste0("Writing ",name," image to disk"))
  ggsave(pl,
         filename = paste0("outputs/figs/ts_wordcloud_mnC_", safe_name, "_v2.pdf"),
         device = "pdf",
         width = 18.12, height = 11.56, units = "in");rm(name, pl)
  message(paste("Saved pdf to: outputs/figs/ts_wordcloud_", safe_name, "_v2.pdf"))
  toc(log=TRUE)
  rm(safe_name)
  }

toc(log = TRUE)

