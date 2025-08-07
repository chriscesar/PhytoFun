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
dfsol <- readRDS(file="outputs/Phyto_2000_2025_with_Lifeforms_USE_updateCalcs.Rdat") %>% 
  janitor::clean_names() %>% 
  dplyr::filter(str_detect(water_body, "^SOLENT")) %>% 
  dplyr::select(-c(
    region, area,sample_time,sample_date_time,
    sample_reason_primary,analysis_id,
    replicate_code,settling_time,taxon_qualifier,
    colonies_per_litre_millilitre,sample_ngr,
    sample_method,wfd_waterbody_id,catchment,
    version,analysis_type,analysis_method,date_of_analysis,
    distance_to_wb,eastings,northings,aphia_id,
    aphia_id_2, scientificname,status, unacceptreason,
    taxon_rank_id,rank,parent_name_usage_id_x,
    is_marine,is_brackish,is_freshwater,is_terrestrial,is_extinct,
    x2,parent_name_usage_id_y,taxa_reported,valid_name,size_class,
    qa_flag))
    
toc(log=TRUE)

dfsol %>% 
  dplyr::select(
    river_basi:surveillan,
    cells_per_litre,
    cells_per_m3,
    #mean_vol_per_cell_um3,
    #median_vol_per_cell_um3,
    #mean_c_per_cell_pg_c,
    #mean_c_per_cell_ug_c,
    #mean_c_per_cell_mg_c,
    #median_c_per_cell_pg_c,
    #median_c_per_cell_ug_c,
    #median_c_per_cell_mg_c,
    #tot_mn_vol_um3,
    #tot_md_vol_um3,
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
    ) %>% 
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
    .groups = "drop") -> dfsummary

write.csv(
  dfsummary,
  file = "outputs/Phyto_Solent.csv",
  row.names = FALSE
  )
toc(log=TRUE)

dfsummary %>% 
  dplyr::filter(sample_date > "2008-01-01") %>% 
  # dplyr::filter(!is.na(biosys_code)) %>% 
  ggplot(., aes(
    x = sample_date,
    y = log10(mn_carbon_tot_mg_c_per_m3)
    )) +
  geom_point()+
  geom_hline(yintercept = log10(7.5),lty=2)+
  #facet_wrap(. ~ wb_name) +
  geom_smooth(method = "gam")+
  labs(title = "Log10 carbon content per sample over time in the Solent",
       caption = "Blue line indicates GAM smoother\nDashed line represents log10(7.5)")+
  theme(
    axis.title.x = element_blank(),
    axis.text.x =  element_text(face=2),
    axis.title.y = element_text(face=2)
  )
  
