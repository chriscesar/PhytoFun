# Solent_investigate.R ####

# Investigate Solent data

# load packages ####
ld_pkgs <- c("tidyverse","tictoc","stringr")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

source("R/000setup.R")

# load data ####
tic("Load data")
# df0_raw <- readRDS(file="outputs/Phyto_2000_2025_with_Lifeforms_USE_updateCalcs.Rdat") %>% 
#   janitor::clean_names() #
# df0 <- df0_raw %>% 
#   # dplyr::filter(str_detect(water_body, "^SOLENT")) %>% 
#   dplyr::select(-c(
#     region, area,sample_time,sample_date_time,
#     sample_reason_primary,analysis_id,
#     replicate_code,settling_time,taxon_qualifier,
#     colonies_per_litre_millilitre,sample_ngr,
#     sample_method,wfd_waterbody_id,catchment,
#     version,analysis_type,analysis_method,date_of_analysis,
#     distance_to_wb,eastings,northings,aphia_id,
#     aphia_id_2, scientificname,status, unacceptreason,
#     taxon_rank_id,rank,parent_name_usage_id_x,
#     is_marine,is_brackish,is_freshwater,is_terrestrial,is_extinct,
#     x2,parent_name_usage_id_y,taxa_reported,valid_name,size_class,
#     qa_flag))
df0_raw <- readRDS("outputs/Phyto_2000_2025_with_Lifeforms_USE_updateCalcs_updateWB.Rdat")
toc(log=TRUE)

tic("Sum carbon across taxa")
df0_raw %>%
  dplyr::select(
    rbd,
    wbid, 
    wb_name,
    wb_type,
    hm_desig,
    area_name,
    man_catch,
    op_catch,
    typology,
    site_id, biosys_code, biosys_short, site_station_name,
    eastings,northings,
    sample_date,
    replicate_code, settling_time, sample_method, 
    colonies_per_litre_millilitre, analysis_type, analysis_method,
    date_of_analysis, surveillan,
    cells_per_litre,
    cells_per_m3,
    #mean_vol_per_cell_um3
    #median_vol_per_cell_um3
    #mean_c_per_cell_pg_c
    #mean_c_per_cell_ug_c
    #mean_c_per_cell_mg_c
    #median_c_per_cell_pg_c
    #median_c_per_cell_ug_c
    #median_c_per_cell_mg_c
    #tot_mn_vol_um3
    #tot_md_vol_um3
    mn_carbon_tot_pg_c_per_litre,
    mn_carbon_tot_ug_c_per_litre,
    mn_carbon_tot_mg_c_per_litre,
    mn_carbon_tot_pg_c_per_m3,
    mn_carbon_tot_ug_c_per_m3,
    mn_carbon_tot_mg_c_per_m3,
    md_carbon_tot_pg_c_per_litre,
    md_carbon_tot_ug_c_per_litre,
    md_carbon_tot_mg_c_per_litre,
    md_carbon_tot_pg_c_per_m3,
    md_carbon_tot_ug_c_per_m3,
    md_carbon_tot_mg_c_per_m3) %>% 
  dplyr::group_by(dplyr::across(-c(
    cells_per_litre,
    cells_per_m3,
    mn_carbon_tot_pg_c_per_litre,
    mn_carbon_tot_ug_c_per_litre,
    mn_carbon_tot_mg_c_per_litre,
    mn_carbon_tot_pg_c_per_m3,
    mn_carbon_tot_ug_c_per_m3,
    mn_carbon_tot_mg_c_per_m3,
    md_carbon_tot_pg_c_per_litre,
    md_carbon_tot_ug_c_per_litre,
    md_carbon_tot_mg_c_per_litre,
    md_carbon_tot_pg_c_per_m3,
    md_carbon_tot_ug_c_per_m3,
    md_carbon_tot_mg_c_per_m3
    ))) %>% 
  dplyr::summarise(
    cells_per_litre = sum(cells_per_litre, na.rm = FALSE),
    cells_per_m3 = sum(cells_per_m3, na.rm = FALSE),
    mn_carbon_tot_pg_c_per_litre = sum(mn_carbon_tot_pg_c_per_litre, na.rm=TRUE),
    mn_carbon_tot_ug_c_per_litre = sum(mn_carbon_tot_ug_c_per_litre, na.rm=TRUE),
    mn_carbon_tot_mg_c_per_litre = sum(mn_carbon_tot_mg_c_per_litre, na.rm=TRUE),
    mn_carbon_tot_pg_c_per_m3 = sum(mn_carbon_tot_pg_c_per_m3, na.rm=TRUE),
    mn_carbon_tot_ug_c_per_m3 = sum(mn_carbon_tot_ug_c_per_m3, na.rm=TRUE),
    mn_carbon_tot_mg_c_per_m3 = sum(mn_carbon_tot_mg_c_per_m3, na.rm=TRUE),
    md_carbon_tot_pg_c_per_litre = sum(md_carbon_tot_pg_c_per_litre, na.rm=TRUE),
    md_carbon_tot_ug_c_per_litre = sum(md_carbon_tot_ug_c_per_litre, na.rm=TRUE),
    md_carbon_tot_mg_c_per_litre = sum(md_carbon_tot_mg_c_per_litre, na.rm=TRUE),
    md_carbon_tot_pg_c_per_m3 = sum(md_carbon_tot_pg_c_per_m3, na.rm=TRUE),
    md_carbon_tot_ug_c_per_m3 = sum(md_carbon_tot_ug_c_per_m3, na.rm=TRUE),
    md_carbon_tot_mg_c_per_m3 = sum(md_carbon_tot_mg_c_per_m3, na.rm=TRUE),
    .groups = "drop") %>%# names()
  dplyr::mutate(rbd = factor(rbd,
                             levels = c(
                               "Northumbria",
                               "Humber",
                               "Anglian",
                               "Thames",
                               "South East",
                               "South West",
                               "Severn",
                               "North West",
                               "Solway Tweed"
                               )
                             )) %>% 
  dplyr::mutate(wb_name = factor(wb_name,
                                 levels = unique(wb_name[order(rbd)])
                                 )
                ) -> dfsummary

write.csv(
  dfsummary,
  file = "outputs/Phyto_by_sample.csv",
  row.names = FALSE
  )
toc(log=TRUE)

### SOLENT
dfsummary %>% 
  dplyr::filter(sample_date > "2008-01-01") %>% 
  dplyr::filter(wb_name == "Solent") %>% 
  # dplyr::filter(!is.na(biosys_code)) %>% 
  ggplot(., aes(
    x = sample_date,
    y = log10(md_carbon_tot_mg_c_per_m3)
    )) +
  geom_point()+
  geom_hline(yintercept = log10(7.5),lty=2)+
  #facet_wrap(. ~ wb_name) +
  geom_smooth(method = "gam")+
  labs(title = "Log10 carbon content per phytoplankton sample over time in the Solent",
       caption = "Blue line indicates GAM smoother\nDashed line represents log10(7.5)",
       y = bquote(Log[10]~Total~carbon~content~(as~mg~carbon~m^-3)))+
  theme(
    axis.title.x = element_blank(),
    axis.text.x =  element_text(face=2),
    axis.title.y = element_text(face=2)
  )
  
### Thames
dfsummary %>% 
  dplyr::filter(sample_date > "2008-01-01") %>% 
  dplyr::filter(stringr::str_detect(wb_name, stringr::regex("Thames", ignore_case = TRUE))) %>% 
  # dplyr::filter(!is.na(biosys_code)) %>% 
  ggplot(., aes(
    x = sample_date,
    y = log10(md_carbon_tot_mg_c_per_m3)
  )) +
  geom_point()+
  geom_hline(yintercept = log10(7.5),lty=2)+
  #facet_wrap(. ~ wb_name) +
  geom_smooth(method = "gam")+
  facet_wrap(.~wb_name)+
  labs(title = "Log10 carbon content per phytoplankton sample over time in Thames water bodies",
       caption = "Blue line indicates GAM smoother\nDashed line represents log10(7.5)",
       y = bquote(Log[10]~Total~carbon~content~(as~mg~carbon~m^-3)))+
  theme(
    axis.title.x = element_blank(),
    strip.text = element_text(face=2),
    axis.text.x =  element_text(face=2),
    axis.title.y = element_text(face=2)
  )

### Essex & Kent North
dfsummary %>% 
  # dplyr::filter(sample_date > "2008-01-01") %>% 
  dplyr::filter(stringr::str_detect(
    wb_name,
    stringr::regex("^(Thames|Essex|Kent North)", ignore_case = TRUE)
  )) %>% 
  # dplyr::filter(!is.na(biosys_code)) %>% 
  ggplot(., aes(
    x = sample_date,
    y = log10(md_carbon_tot_mg_c_per_m3)
  )) +
  geom_point(
    #aes(colour = wb_name)
  )+
  geom_hline(yintercept = log10(7.5),lty=2)+
  #facet_wrap(. ~ wb_name) +
  geom_smooth(
    # aes(colour = wb_name),
    method = "gam",
    se=FALSE
    )+
  #facet_wrap(.~wb_name)+
  labs(title = "Log10 carbon content per phytoplankton sample over time in Thames water bodies",
       caption = "Blue line indicates GAM smoother\nDashed line represents log10(7.5)",
       y = bquote(Log[10]~Total~carbon~content~(as~mg~carbon~m^-3)))+
  theme(
    axis.title.x = element_blank(),
    strip.text = element_text(face=2),
    axis.text.x =  element_text(face=2),
    axis.title.y = element_text(face=2)
  )

# Total C all WBs ####
# throw away WBs with fewer than this number of samples:
number <- 500

dfsummary %>% 
  dplyr::group_by(wb_name) %>%
  dplyr::filter(n() >= number) %>%
  dplyr::filter(sample_date > "2008-01-01") %>% 
  ggplot(.,
         aes(
           x = sample_date,
           y = log10(md_carbon_tot_mg_c_per_m3)
         )) +
  geom_point(
    # aes(colour = wb_name),show.legend = FALSE,
    alpha=0.25
    )+
  # facet_wrap(.~wb_name) +
  facet_wrap(.~rbd) +
  labs(title = bquote(bold(Log[10]~carbon~content~per~phytoplankton~sample~over~time~'in'~selected~EA~water~bodies)),
       caption = "Blue line indicates GAM smoother\nDashed line represents log10(7.5)",
       y = bquote(bold(Log[10]~Total~carbon~content~(as~mg~carbon~m^-3))))+
  geom_hline(yintercept = log10(7.5), lty = 2)+
  geom_smooth(
    aes(colour = wb_name),show.legend = FALSE,
    method = "gam",
    se = FALSE
  ) +
  geom_smooth(
    linewidth = 2,
    method = "gam",
    se = FALSE
  ) +
  theme(
    axis.title.x = element_blank(),
    strip.text = element_text(face=2)
  )
         

## diatom - dinof ratio
df0_raw %>% 
  # remove non-diatom/non-dinoflagellate
  dplyr::filter(phyto_lf == "Diatom"|phyto_lf == "Dinoflagellate") %>% 
  dplyr::select(
    rbd,
    wbid, 
    wb_name,
    wb_type,
    hm_desig,
    area_name,
    man_catch,
    op_catch,
    typology,
    site_id, biosys_code, biosys_short, site_station_name,
    sample_date,
    phyto_lf,
    replicate_code, settling_time, sample_method, 
    colonies_per_litre_millilitre, analysis_type, analysis_method,
    date_of_analysis, surveillan,
    cells_per_litre,
    cells_per_m3,
    mn_carbon_tot_pg_c_per_litre,
    mn_carbon_tot_ug_c_per_litre,
    mn_carbon_tot_mg_c_per_litre,
    mn_carbon_tot_pg_c_per_m3,
    mn_carbon_tot_ug_c_per_m3,
    mn_carbon_tot_mg_c_per_m3,
    md_carbon_tot_pg_c_per_litre,
    md_carbon_tot_ug_c_per_litre,
    md_carbon_tot_mg_c_per_litre,
    md_carbon_tot_pg_c_per_m3,
    md_carbon_tot_ug_c_per_m3,
    md_carbon_tot_mg_c_per_m3) %>% 
  dplyr::group_by(dplyr::across(-c(
    cells_per_litre,
    cells_per_m3,
    mn_carbon_tot_pg_c_per_litre,
    mn_carbon_tot_ug_c_per_litre,
    mn_carbon_tot_mg_c_per_litre,
    mn_carbon_tot_pg_c_per_m3,
    mn_carbon_tot_ug_c_per_m3,
    mn_carbon_tot_mg_c_per_m3,
    md_carbon_tot_pg_c_per_litre,
    md_carbon_tot_ug_c_per_litre,
    md_carbon_tot_mg_c_per_litre,
    md_carbon_tot_pg_c_per_m3,
    md_carbon_tot_ug_c_per_m3,
    md_carbon_tot_mg_c_per_m3
  ))) %>% 
  dplyr::summarise(
    cells_per_litre = sum(cells_per_litre, na.rm = FALSE),
    cells_per_m3 = sum(cells_per_m3, na.rm = FALSE),
    mn_carbon_tot_pg_c_per_litre = sum(mn_carbon_tot_pg_c_per_litre, na.rm=TRUE),
    mn_carbon_tot_ug_c_per_litre = sum(mn_carbon_tot_ug_c_per_litre, na.rm=TRUE),
    mn_carbon_tot_mg_c_per_litre = sum(mn_carbon_tot_mg_c_per_litre, na.rm=TRUE),
    mn_carbon_tot_pg_c_per_m3 = sum(mn_carbon_tot_pg_c_per_m3, na.rm=TRUE),
    mn_carbon_tot_ug_c_per_m3 = sum(mn_carbon_tot_ug_c_per_m3, na.rm=TRUE),
    mn_carbon_tot_mg_c_per_m3 = sum(mn_carbon_tot_mg_c_per_m3, na.rm=TRUE),
    md_carbon_tot_pg_c_per_litre = sum(md_carbon_tot_pg_c_per_litre, na.rm=TRUE),
    md_carbon_tot_ug_c_per_litre = sum(md_carbon_tot_ug_c_per_litre, na.rm=TRUE),
    md_carbon_tot_mg_c_per_litre = sum(md_carbon_tot_mg_c_per_litre, na.rm=TRUE),
    md_carbon_tot_pg_c_per_m3 = sum(md_carbon_tot_pg_c_per_m3, na.rm=TRUE),
    md_carbon_tot_ug_c_per_m3 = sum(md_carbon_tot_ug_c_per_m3, na.rm=TRUE),
    md_carbon_tot_mg_c_per_m3 = sum(md_carbon_tot_mg_c_per_m3, na.rm=TRUE),
    .groups = "drop") %>%# names()
  dplyr::mutate(rbd = factor(rbd,
                             levels = c(
                               "Northumbria",
                               "Humber",
                               "Anglian",
                               "Thames",
                               "South East",
                               "South West",
                               "Severn",
                               "North West",
                               "Solway Tweed"
                             )
  )) %>% #names()
  dplyr::mutate(wb_name = factor(wb_name,
                                 levels = unique(wb_name[order(rbd)])
                                 )
                ) %>% 
  # widen for calculation of ratio
  tidyr::pivot_wider(
    names_from = phyto_lf,
    values_from = c(cells_per_litre,
    cells_per_m3,
    mn_carbon_tot_pg_c_per_litre,
    mn_carbon_tot_ug_c_per_litre,
    mn_carbon_tot_mg_c_per_litre,
    mn_carbon_tot_pg_c_per_m3,
    mn_carbon_tot_ug_c_per_m3,
    mn_carbon_tot_mg_c_per_m3,
    md_carbon_tot_pg_c_per_litre,
    md_carbon_tot_ug_c_per_litre,
    md_carbon_tot_mg_c_per_litre,
    md_carbon_tot_pg_c_per_m3,
    md_carbon_tot_ug_c_per_m3,
    md_carbon_tot_mg_c_per_m3)
    ) %>% 
  dplyr::mutate(
    ratio_md_carbon_tot_mg_c_per_m3 = (md_carbon_tot_mg_c_per_m3_Diatom/md_carbon_tot_mg_c_per_m3_Dinoflagellate)) %>% 
  # remove NA and Inf values
  dplyr::filter(!is.na(ratio_md_carbon_tot_mg_c_per_m3)) %>% 
  dplyr::filter(!is.infinite(ratio_md_carbon_tot_mg_c_per_m3)) %>% 
  mutate(
    gtr0 = case_when(
      log10(ratio_md_carbon_tot_mg_c_per_m3) > 0  ~ "Above",
      log10(ratio_md_carbon_tot_mg_c_per_m3) < 0  ~ "Below",
      log10(ratio_md_carbon_tot_mg_c_per_m3) == 0 ~ "Zero")) %>% 
  group_by(wb_name) %>% dplyr::filter(n() >= number) %>% ungroup() %>% 
  dplyr::filter(sample_date > "2008-01-01") %>% 
  ggplot(aes(
    x = sample_date,
    y = log10(ratio_md_carbon_tot_mg_c_per_m3)
  )
  )+
  geom_point(
    aes(colour = gtr0), show.legend = FALSE,
    alpha = 0.25
  )+
  facet_wrap(.~rbd)+
  geom_smooth(
    method="gam",
    se=TRUE
    )+
  geom_hline(yintercept = log10(1))+
  labs(title =
         bquote(bold(Ratio~of~log[10]~carbon~content~'in'~diatoms~taxa~to~that~'in'~donoflagellates)),
       y= bquote(bold(Log[10]~ratio~of~diatom~carbon~to~dinoflagellate~carbon)
                 ),
       caption="Blue line indicates GAM smoother\nOnly samples where both diatoms and dinoflagellates were recorded are displayed",
       # subtitle = "Just for fun"
       )+
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(face=2),
    strip.text = element_text(face=2)
  )
###########

# by lifeform ####
tic("Sum carbon across taxa")
df0_raw %>%
  dplyr::select(
    rbd,
    wbid, 
    wb_name,
    wb_type,
    hm_desig,
    area_name,
    man_catch,
    op_catch,
    typology,
    site_id, biosys_code, biosys_short, site_station_name,
    eastings,northings,
    sample_date,
    replicate_code, 
    phyto_lf,
    colonies_per_litre_millilitre, analysis_type, analysis_method,
    date_of_analysis, surveillan,
    cells_per_litre,
    cells_per_m3,
    #mean_vol_per_cell_um3
    #median_vol_per_cell_um3
    #mean_c_per_cell_pg_c
    #mean_c_per_cell_ug_c
    #mean_c_per_cell_mg_c
    #median_c_per_cell_pg_c
    #median_c_per_cell_ug_c
    #median_c_per_cell_mg_c
    #tot_mn_vol_um3
    #tot_md_vol_um3
    mn_carbon_tot_pg_c_per_litre,
    mn_carbon_tot_ug_c_per_litre,
    mn_carbon_tot_mg_c_per_litre,
    mn_carbon_tot_pg_c_per_m3,
    mn_carbon_tot_ug_c_per_m3,
    mn_carbon_tot_mg_c_per_m3,
    md_carbon_tot_pg_c_per_litre,
    md_carbon_tot_ug_c_per_litre,
    md_carbon_tot_mg_c_per_litre,
    md_carbon_tot_pg_c_per_m3,
    md_carbon_tot_ug_c_per_m3,
    md_carbon_tot_mg_c_per_m3) %>% 
  dplyr::group_by(dplyr::across(-c(
    cells_per_litre,
    cells_per_m3,
    mn_carbon_tot_pg_c_per_litre,
    mn_carbon_tot_ug_c_per_litre,
    mn_carbon_tot_mg_c_per_litre,
    mn_carbon_tot_pg_c_per_m3,
    mn_carbon_tot_ug_c_per_m3,
    mn_carbon_tot_mg_c_per_m3,
    md_carbon_tot_pg_c_per_litre,
    md_carbon_tot_ug_c_per_litre,
    md_carbon_tot_mg_c_per_litre,
    md_carbon_tot_pg_c_per_m3,
    md_carbon_tot_ug_c_per_m3,
    md_carbon_tot_mg_c_per_m3
  ))) %>% 
  dplyr::summarise(
    cells_per_litre = sum(cells_per_litre, na.rm = FALSE),
    cells_per_m3 = sum(cells_per_m3, na.rm = FALSE),
    mn_carbon_tot_pg_c_per_litre = sum(mn_carbon_tot_pg_c_per_litre, na.rm=TRUE),
    mn_carbon_tot_ug_c_per_litre = sum(mn_carbon_tot_ug_c_per_litre, na.rm=TRUE),
    mn_carbon_tot_mg_c_per_litre = sum(mn_carbon_tot_mg_c_per_litre, na.rm=TRUE),
    mn_carbon_tot_pg_c_per_m3 = sum(mn_carbon_tot_pg_c_per_m3, na.rm=TRUE),
    mn_carbon_tot_ug_c_per_m3 = sum(mn_carbon_tot_ug_c_per_m3, na.rm=TRUE),
    mn_carbon_tot_mg_c_per_m3 = sum(mn_carbon_tot_mg_c_per_m3, na.rm=TRUE),
    md_carbon_tot_pg_c_per_litre = sum(md_carbon_tot_pg_c_per_litre, na.rm=TRUE),
    md_carbon_tot_ug_c_per_litre = sum(md_carbon_tot_ug_c_per_litre, na.rm=TRUE),
    md_carbon_tot_mg_c_per_litre = sum(md_carbon_tot_mg_c_per_litre, na.rm=TRUE),
    md_carbon_tot_pg_c_per_m3 = sum(md_carbon_tot_pg_c_per_m3, na.rm=TRUE),
    md_carbon_tot_ug_c_per_m3 = sum(md_carbon_tot_ug_c_per_m3, na.rm=TRUE),
    md_carbon_tot_mg_c_per_m3 = sum(md_carbon_tot_mg_c_per_m3, na.rm=TRUE),
    .groups = "drop") %>%# names()
  dplyr::mutate(rbd = factor(rbd,
                             levels = c(
                               "Northumbria",
                               "Humber",
                               "Anglian",
                               "Thames",
                               "South East",
                               "South West",
                               "Severn",
                               "North West",
                               "Solway Tweed"
                             )
  )) %>% 
  dplyr::mutate(wb_name = factor(wb_name,
                                 levels = unique(wb_name[order(rbd)])
  )
  ) -> dfsummary_lf

dfsummary_lf %>%
  dplyr::filter(!is.na(phyto_lf)) %>% 
  #dplyr::filter(rbd == "South East") %>% 
  dplyr::filter(phyto_lf %in% c("Diatom","Dinoflagellate")) %>% 
  dplyr::group_by(across(-c(cells_per_litre,
                            cells_per_m3,
                            mn_carbon_tot_pg_c_per_litre,
                            mn_carbon_tot_ug_c_per_litre,
                            mn_carbon_tot_mg_c_per_litre,
                            mn_carbon_tot_pg_c_per_m3,
                            mn_carbon_tot_ug_c_per_m3,
                            mn_carbon_tot_mg_c_per_m3,
                            md_carbon_tot_pg_c_per_litre,
                            md_carbon_tot_ug_c_per_litre,
                            md_carbon_tot_mg_c_per_litre,
                            md_carbon_tot_pg_c_per_m3,
                            md_carbon_tot_ug_c_per_m3,
                            md_carbon_tot_mg_c_per_m3))) %>% 
  dplyr::summarise(
    cells_per_litre = mean(cells_per_litre, na.rm = TRUE),
    cells_per_m3 = mean(cells_per_m3, na.rm = TRUE),
    mn_carbon_tot_pg_c_per_litre = mean(mn_carbon_tot_pg_c_per_litre, na.rm = TRUE),
    mn_carbon_tot_ug_c_per_litre = mean(mn_carbon_tot_ug_c_per_litre, na.rm = TRUE),
    mn_carbon_tot_mg_c_per_litre = mean(mn_carbon_tot_mg_c_per_litre, na.rm = TRUE),
    mn_carbon_tot_pg_c_per_m3 = mean(mn_carbon_tot_pg_c_per_m3, na.rm = TRUE),
    mn_carbon_tot_ug_c_per_m3 = mean(mn_carbon_tot_ug_c_per_m3, na.rm = TRUE),
    mn_carbon_tot_mg_c_per_m3 = mean(mn_carbon_tot_mg_c_per_m3, na.rm = TRUE),
    md_carbon_tot_pg_c_per_litre = mean(md_carbon_tot_pg_c_per_litre, na.rm = TRUE),
    md_carbon_tot_ug_c_per_litre = mean(md_carbon_tot_ug_c_per_litre, na.rm = TRUE),
    md_carbon_tot_mg_c_per_litre = mean(md_carbon_tot_mg_c_per_litre, na.rm = TRUE),
    md_carbon_tot_pg_c_per_m3 = mean(md_carbon_tot_pg_c_per_m3, na.rm = TRUE),
    md_carbon_tot_ug_c_per_m3 = mean(md_carbon_tot_ug_c_per_m3, na.rm = TRUE),
    md_carbon_tot_mg_c_per_m3 = mean(md_carbon_tot_mg_c_per_m3, na.rm = TRUE),
    .groups = "drop"
    ) %>% 
  dplyr::mutate(
    month = lubridate::month(sample_date),
    year = lubridate::year(sample_date)) %>% 
  dplyr::filter(year>2007) %>% 
  dplyr::filter(rbd != "Severn") %>% 
  dplyr::filter(rbd != "Solway Tweed") %>% 
  ggplot(aes(
    x = year,
    y = month
    # z = log10(md_carbon_tot_mg_c_per_m3+1)
    )) +
  # geom_point()
  # stat_contour(geom="raster", interpolate=TRUE)+
  # scale_fill_gradientn(colours=rev(rainbow(10)))+
  # geom_tile(aes(fill=md_carbon_tot_mg_c_per_m3))+
  # facet_wrap(~phyto_lf)
  geom_hline(yintercept = c(1:12),lty = 2, colour="grey")+
  geom_vline(xintercept = seq(from=2007.5, to=2025.5, by=1), lty = 2, colour="grey") +
  geom_tile(aes(fill = log10(cells_per_m3)),show.legend = TRUE)+
  # facet_wrap(.~phyto_lf,ncol = 1)+
  facet_grid(rbd~phyto_lf)+
  scale_fill_gradient(low = "lightskyblue2", high = "red")+
  # scale_fill_gradientn(colours = c("blue", "white", "red"),
  #                      values = scales::rescale(c(0,5,10)))+
  labs(title = "Heatmaps of plankton abundances",
       subtitle = "Data expressed as log10(cells per m3)",
       x="",
       y="Month")+
  guides(fill=guide_legend(title=bquote(bold(Log[10]~Cells~m^-3)),
                           reverse = TRUE))+
  scale_y_continuous(breaks= scales::pretty_breaks())+
  coord_cartesian(xlim = c(2008,NA),
                  ylim = c(1,12))+
theme(
  axis.title = element_text(face=2),
  axis.text = element_text(face=2),
  strip.text = element_text(face=2),
  legend.title = element_text(face=2)
  )

# by taxon ####

plnk <- "Other phytoplankton"
mxCells <- log10(max(df0_raw$cells_per_m3[df0_raw$name_use == plnk],
                     na.rm = TRUE))

df0_raw %>%
  # dplyr::filter(grepl(plnk,name_use)) %>%
  dplyr::filter(grepl(paste0("^",plnk),name_use)) %>%
  dplyr::mutate(rbd_wb_biosys = paste0(rbd,"_",wb_name,"_",biosys_code)) %>% 
  dplyr::mutate(rbd_wb_biosys = factor(rbd_wb_biosys,
                                       levels = unique(rbd_wb_biosys[order(rbd)]))) %>%
  dplyr::select(
    rbd,
    wbid, 
    wb_name,
    rbd_wb_biosys,
    wb_type,
    hm_desig,
    area_name,
    man_catch,
    op_catch,
    typology,
    site_id, biosys_code, biosys_short, site_station_name,
    eastings,northings,
    sample_date,
    replicate_code, 
    colonies_per_litre_millilitre, analysis_type, analysis_method,
    date_of_analysis, surveillan,
    cells_per_litre,
    cells_per_m3,
    mn_carbon_tot_pg_c_per_litre,
    mn_carbon_tot_ug_c_per_litre,
    mn_carbon_tot_mg_c_per_litre,
    mn_carbon_tot_pg_c_per_m3,
    mn_carbon_tot_ug_c_per_m3,
    mn_carbon_tot_mg_c_per_m3,
    md_carbon_tot_pg_c_per_litre,
    md_carbon_tot_ug_c_per_litre,
    md_carbon_tot_mg_c_per_litre,
    md_carbon_tot_pg_c_per_m3,
    md_carbon_tot_ug_c_per_m3,
    md_carbon_tot_mg_c_per_m3) %>% 
  dplyr::group_by(dplyr::across(-c(
    cells_per_litre,
    cells_per_m3,
    mn_carbon_tot_pg_c_per_litre,
    mn_carbon_tot_ug_c_per_litre,
    mn_carbon_tot_mg_c_per_litre,
    mn_carbon_tot_pg_c_per_m3,
    mn_carbon_tot_ug_c_per_m3,
    mn_carbon_tot_mg_c_per_m3,
    md_carbon_tot_pg_c_per_litre,
    md_carbon_tot_ug_c_per_litre,
    md_carbon_tot_mg_c_per_litre,
    md_carbon_tot_pg_c_per_m3,
    md_carbon_tot_ug_c_per_m3,
    md_carbon_tot_mg_c_per_m3
  ))) %>% 
  dplyr::summarise(
    cells_per_litre = sum(cells_per_litre, na.rm = FALSE),
    cells_per_m3 = sum(cells_per_m3, na.rm = FALSE),
    mn_carbon_tot_pg_c_per_litre = sum(mn_carbon_tot_pg_c_per_litre, na.rm=TRUE),
    mn_carbon_tot_ug_c_per_litre = sum(mn_carbon_tot_ug_c_per_litre, na.rm=TRUE),
    mn_carbon_tot_mg_c_per_litre = sum(mn_carbon_tot_mg_c_per_litre, na.rm=TRUE),
    mn_carbon_tot_pg_c_per_m3 = sum(mn_carbon_tot_pg_c_per_m3, na.rm=TRUE),
    mn_carbon_tot_ug_c_per_m3 = sum(mn_carbon_tot_ug_c_per_m3, na.rm=TRUE),
    mn_carbon_tot_mg_c_per_m3 = sum(mn_carbon_tot_mg_c_per_m3, na.rm=TRUE),
    md_carbon_tot_pg_c_per_litre = sum(md_carbon_tot_pg_c_per_litre, na.rm=TRUE),
    md_carbon_tot_ug_c_per_litre = sum(md_carbon_tot_ug_c_per_litre, na.rm=TRUE),
    md_carbon_tot_mg_c_per_litre = sum(md_carbon_tot_mg_c_per_litre, na.rm=TRUE),
    md_carbon_tot_pg_c_per_m3 = sum(md_carbon_tot_pg_c_per_m3, na.rm=TRUE),
    md_carbon_tot_ug_c_per_m3 = sum(md_carbon_tot_ug_c_per_m3, na.rm=TRUE),
    md_carbon_tot_mg_c_per_m3 = sum(md_carbon_tot_mg_c_per_m3, na.rm=TRUE),
    .groups = "drop") %>%# names()
  dplyr::mutate(rbd = factor(rbd,
                             levels = c(
                               "Northumbria",
                               "Humber",
                               "Anglian",
                               "Thames",
                               "South East",
                               "South West",
                               "Severn",
                               "North West",
                               "Solway Tweed"
                             )
  )) %>% 
  dplyr::mutate(wb_name = factor(wb_name,
                                 levels = unique(wb_name[order(rbd)])
  )
  ) %>% 
  dplyr::mutate(
    month = lubridate::month(sample_date),
    year = lubridate::year(sample_date)) %>% 
  # dplyr::filter(year>2007) %>% # Keep only most recent data
  # dplyr::filter(rbd != "Severn") %>% # throw out data-scarce sites
  # dplyr::filter(rbd != "Solway Tweed") %>% # throw out data-scarce sites
  ggplot(aes(
    x = year,
    y = month
    )) +
  geom_rect(xmin = -Inf,
            xmax = 2030,
            ymin = 0,
            ymax = Inf,
            # fill="lightskyblue2",
            fill="paleturquoise",
            alpha=0.1
            )+
  geom_tile(aes(fill = log10(cells_per_m3)),show.legend = TRUE)+
  geom_hline(yintercept = seq(0.5,12.5,by=1),lty = 2, colour="grey")+
  geom_vline(xintercept = seq(from=1999.5, to=2025.5, by=1), lty = 2, colour="grey") +
  # facet_wrap(.~wb_name)+ #By water body
  facet_wrap(.~rbd)+ #By river basin
  scale_fill_gradient(limits = c(0,mxCells),
                      low = "paleturquoise", high = "orange3")+
  labs(title = paste0("'",plnk,"' abundances"),
       subtitle = bquote(Values~expressed~as~log[10]~cells~~m^-3),
       x="",
       y="Month")+
  guides(fill=guide_legend(title=bquote(bold(Log[10]~Cells~m^-3)),
                           reverse = TRUE))+
  scale_y_continuous(breaks= scales::pretty_breaks())+
  coord_cartesian(
    xlim = c(2007,NA),
    ylim = c(1,12))+
  theme(
    axis.title = element_text(face=2),
    axis.text = element_text(face=2),
    strip.text = element_text(face=2),
    legend.title = element_text(face=2)
  )

#### examine acf plots: Critical slowing down????
df0_raw %>%
  filter(wb_name=="Mersey Mouth") %>%
  filter(biosys_short=="LIV001") %>% 
  filter(name_use == "Pseudo-nitzschia") %>% 
  bayesforecast::ggacf(md_carbon_tot_mg_c_per_m3)
  
