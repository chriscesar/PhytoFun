# phytoWordcloud.R ####
#### Generate wordclouds of life forms and taxon data

ld_pkgs <- c("tidyverse","ggwordcloud", "tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

ppi <- 300
tictoc::tic.clearlog() ##clear log

#### set universals ####
tic()
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
# load data ####

df0 <- readRDS("outputs/Phyto_2000_2025_with_Lifeforms_USE_updateCalcs.Rdat")

df <- df0 %>%
  select(river_basi,wb_name,biosys_code,sample_date,
         # name_use,
         taxa_reported,
         md_carbon_tot_ug_C_per_litre,
         plankton_type:phyto_lf
         ) %>%
  ### drop pre 2011 data
  dplyr::filter(sample_date > "2010-12-31") %>% 
  # define colours
  dplyr::mutate(phyto_lf_col = case_when(
    phyto_lf == "Diatom" ~ "chartreuse",
    phyto_lf == "Dinoflagellate" ~ "dodgerblue",
    phyto_lf == "Protozoa" ~ "blueviolet",
    phyto_lf == "Haptophyte" ~ "burlywood4",
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
                                          "burlywood4","cyan3","sienna","aquamarine2",
                                          "deeppink2","darkgreen","darkorange",
                                          "darkred","darkgrey")
    
  )) %>% 
  dplyr::mutate(phyto_lf = ifelse(is.na(phyto_lf), "Unknown", phyto_lf))
  
          
pal <- c("chartreuse","dodgerblue","blueviolet","burlywood4","cyan3","sienna",
         "aquamarine2","deeppink2","darkgreen","darkorange","darkred","darkgrey")


phyto_col <- c("Diatom" = "chartreuse",
               "Dinoflagellate" = "dodgerblue",
               "Protozoa" = "blueviolet",
               "Haptophyte" = "burlywood4",
               "Cyanobacteria" = "cyan3",
               "Chrysophyte" = "sienna",
               "Chlorophyte" = "aquamarine2",
               "Silicoflagellate" = "deeppink2",
               "Charophyte" = "darkgreen",
               "Ciliate" = "darkorange",
               "Raphidophyte" = "darkred",
               "Unknown" = "darkgrey")

# tic("Generate wordcloud");set.seed(3123);df %>% 
#   # dplyr::filter(river_basi == "South East") %>% 
#   dplyr::filter(!is.na(phyto_lf)) %>% 
#   dplyr::group_by(river_basi,
#                   taxa_reported,phyto_lf) %>% 
#   dplyr::summarise(sum = sum(md_carbon_tot_ug_C_per_litre,
#                              na.rm=TRUE),
#                    .groups = "drop") %>% 
#   dplyr::ungroup() %>% 
#   dplyr::mutate(tot = sum(sum)) %>% 
#   dplyr::mutate(prop = sum/tot) %>% 
#   ggplot(.,aes(
#     label = taxa_reported,
#     size = prop,
#     colour = phyto_lf
#   ))+
#   geom_text_wordcloud_area()+
#   scale_size_area(max_size = 85,#originally 150
#                   trans = power_trans(1/.7)) + #trans = power_trans(1/.7) "better fit human area perception"
#   theme_minimal()+
#   facet_wrap(. ~ river_basi, scales = "free")+
#   theme(strip.text = element_text(face = "bold",
#                                   size = 14,
#                                   colour = "red")) -> pl1
# 
# toc(log=TRUE)
# 
# tic("Write image")
# ggsave(pl1,filename = "outputs/figs/wordcloud01.png")
# toc(log = TRUE)
# 
# ## by year ####
# tic("Generate wordcloud");set.seed(3123);df %>% 
#   dplyr::filter(river_basi == "South East") %>% 
#   dplyr::mutate(year = lubridate::year(sample_date)) %>% 
#   dplyr::filter(!is.na(phyto_lf)) %>% #View()
#   dplyr::group_by(year,
#                   taxa_reported,phyto_lf) %>% 
#   dplyr::summarise(sum = sum(md_carbon_tot_ug_C_per_litre,
#                              na.rm=TRUE),
#                    .groups = "drop") %>% 
#   dplyr::ungroup() %>% 
#   dplyr::mutate(tot = sum(sum)) %>% ### proportions & totals by year? (add group_by?)
#   dplyr::mutate(prop = sum/tot) %>% 
#   ggplot(.,aes(
#     label = taxa_reported,
#     size = prop,
#     colour = phyto_lf
#   ))+
#   ggwordcloud::geom_text_wordcloud_area(rm_outside = TRUE)+
#   ggplot2::scale_size_area(
#     max_size = 50,#originally 150
#     # max_size = 85,#originally 150
#     trans = power_trans(1/.7)) + #trans = power_trans(1/.7) "better fit human area perception"
#   theme_minimal()+
#   labs(title = "South East")+
#   facet_wrap(. ~ year, scales = "free")+
#   theme(strip.text = element_text(face = "bold",
#                                   size = 14,
#                                   colour = "red")) -> pl2
# 
# toc(log=TRUE)
# 
# tic("Write image")
# ggsave(pl2,filename = "outputs/figs/wordcloudYear.png")
# toc(log = TRUE)

## Regional plots faceted by year ####
### Total SUMMED carbon ####
# tic("Generate wordcloud");set.seed(3123);for(name in unique(df$river_basi)){
#   
#   message(paste("Generating wordcloud for:", name))
#   
#   # Sanitize filename
#   safe_name <- make.names(name)
#   
#   df %>%
#     dplyr::filter(river_basi == name) %>%
#     dplyr::mutate(year = lubridate::year(sample_date)) %>%
#     dplyr::filter(!is.na(phyto_lf)) %>% #View()
#     dplyr::group_by(year,
#                     taxa_reported,phyto_lf) %>%
#     dplyr::summarise(sum = sum(md_carbon_tot_ug_C_per_litre,
#                                na.rm=TRUE),
#                      .groups = "drop") %>%
#     dplyr::ungroup() %>% 
#     group_by(year) %>%
#     dplyr::mutate(tot = sum(sum)) %>%
#     dplyr::mutate(prop = sum/tot) %>% #View()### proportions & totals by year? (add group_by?)
#     ungroup() %>% 
#     ggplot(.,aes(
#       label = taxa_reported,
#       size = prop,
#       colour = phyto_lf
#       ))+
#     ggwordcloud::geom_text_wordcloud_area(rm_outside = TRUE, show.legend = TRUE)+
#     ggplot2::scale_size_area(
#       max_size = 50,#originally 150
#       # max_size = 85,#originally 150
#       trans = power_trans(1/.7)) + #trans = power_trans(1/.7) "better fit human area perception"
#     theme_minimal()+
#     labs(title = name,
#          caption = "Text size reflects proportion of annual total carbon content provided by that taxon")+
#     facet_wrap(. ~ year, scales = "free")+
#     guides(
#       size = "none",  # removes size legend
#       colour = guide_legend(nrow = 2, byrow = TRUE)  # keeps colour legend in 2 rows
#     )+
#     theme(strip.text = element_text(face = "bold",
#                                     size = 14,
#                                     colour = "red"),
#           title = element_text(face = "bold", size=16),
#           legend.position = "bottom",
#           legend.text = element_text(size = 10),
#           legend.title = element_blank()
#           ) -> pl
#   
#   ggsave(pl,
#          filename = paste0("outputs/figs/ts_wordcloud_", safe_name, ".png"),
#          width = 18.12, height = 11.56, units = "in");rm(pl)
#   
#   message(paste("Saved plot to: outputs/figs/ts_wordcloud_totC_", safe_name, ".png"))
#   
#   }
# 
# toc(log=TRUE)

### Mean carbon ####
all_years <- 2011:2015
tic("Generate wordcloud");set.seed(3123);for(name in unique(df$river_basi)){
  
  message(paste("Generating wordcloud for:", name))
  
  # Sanitize filename
  safe_name <- make.names(name)
  
  df %>%
    dplyr::filter(river_basi == name) %>%
    dplyr::mutate(year = lubridate::year(sample_date)) %>%
    dplyr::filter(!is.na(phyto_lf)) %>% #View()
    dplyr::group_by(year,
                    taxa_reported,phyto_lf) %>%
    dplyr::summarise(mean = mean(md_carbon_tot_ug_C_per_litre,
                                 na.rm=TRUE),
                     .groups = "drop") %>%
    dplyr::ungroup() %>% 
    # removes NaN values
    dplyr::filter(!is.nan(mean)) %>% 
    ## ensure that all years are included
    tidyr::complete(year = all_years, taxa_reported, phyto_lf, fill = list(mean = 0)) %>%
    # sum mean carbon contents by year & calculate proportion
    dplyr::group_by(year) %>%
    dplyr::mutate(tot = sum(mean,na.rm = TRUE)) %>%
    dplyr::mutate(prop = mean/tot) %>% #View()### proportions & totals by year? (add group_by?)
    ungroup() %>% 
    ggplot(.,aes(
      label = taxa_reported,
      size = prop,
      colour = phyto_lf
    ))+
    ggwordcloud::geom_text_wordcloud_area(rm_outside = TRUE, show.legend = TRUE)+
    ggplot2::scale_size_area(
      max_size = 50,#originally 150
      # max_size = 85,#originally 150
      trans = power_trans(1/.7)) + #trans = power_trans(1/.7) "better fit human area perception"
    theme_minimal()+
    labs(title = name,
         caption = "\n
         Text size reflects proportion of mean carbon content provided by that taxon")+
    facet_wrap(. ~ year, scales = "free")+
    # scale_color_discrete(pal)+
    scale_colour_manual(values = phyto_col)+
    guides(
      size = "none",  # removes size legend
      colour = guide_legend(nrow = 2, byrow = TRUE)  # keeps colour legend in 2 rows
    )+
    theme(strip.text = element_text(face = "bold",
                                    size = 14,
                                    colour = "red"),
          title = element_text(face = "bold", size=18),
          legend.position = "bottom",
          legend.text = element_text(size = 10),
          legend.title = element_blank()
    ) -> pl

  # Save PDF
  ggsave(pl,
         filename = paste0("outputs/figs/ts_wordcloud_mnC_", safe_name, ".pdf"),
         device = "pdf",
         width = 18.12, height = 11.56, units = "in")
  
  message(paste("Saved pdf to: outputs/figs/ts_wordcloud_", safe_name, ".pdf"))
  }

toc(log=TRUE)

#############################################
################## To Do: ###################
######### MEAN PLANKTON ABUNDANCES ##########
#############################################
#############################################
#############################################

