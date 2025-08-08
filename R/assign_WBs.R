# assign_WBs.R ####
### assigns sites to water bodies

# load packages ####
ld_pkgs <- c("tidyverse","tictoc","stringr","sf","maps","nngeo")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)

rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

source("R/000setup.R")

# load data & convert to spatial ####
tic("Load data")
df0_raw <- readRDS(file="outputs/Phyto_2000_2025_with_Lifeforms_USE_updateCalcs.Rdat") %>% 
  janitor::clean_names() #
toc(log=TRUE)

tic("extract unique Site Names, convert to spatial & assign to WB")
## extract unique Site Names & assign to WB
df0_raw %>% 
  dplyr::select(site_id,ngr,eastings,northings) %>% 
  dplyr::distinct() -> sites

# convert to spatial df
sites_sf <- sf::st_as_sf(sites,
                         coords = c("eastings","northings"),
                         crs = 27700)

# WB shapefile
base_WBs <- sf::read_sf(paste0(GISfol,"C3_WFD_Waterbody/","EnglandTRAC_C3_ReducedFields.shp"))

# Function to find nearest polygon and distance
nearest_polygon <- function(point, polygons) {
  distances <- st_distance(point, polygons)
  nearest_index <- which.min(distances)
  nearest_poly <- polygons[nearest_index, ]
  distance <- as.numeric(distances[nearest_index])
  return(data.frame(nearest_poly, distance = distance))
}


tic("Join")
# Apply the function to each point
joined_data <- do.call(rbind, lapply(1:nrow(sites_sf), function(i) {
  nearest_polygon(sites_sf[i, ], base_WBs)
}))
toc(log=TRUE)

toc(log = TRUE)

tic("Export")
# Join point data with the joined data
joined_data_out <- cbind(sites_sf, joined_data) %>% as.data.frame(.) %>% 
  dplyr::filter(!is.na(site_id)) %>% 
  dplyr::select(site_id,WBID,
                WB_Name,WB_Type,HM_Desig,Area_Name,Man_Catch,
                Op_Catch,Typology,RBD,distance) %>% 
  dplyr::distinct() %>% 
  janitor::clean_names() 

saveRDS(joined_data_out, file = "outputs/phytoSitesWBs.Rdat")
write.csv(joined_data_out, file = "outputs/phytoSitesWBs.csv",row.names = FALSE)
toc(log = TRUE)

tic("Join phyto data to 'correct' site-WB info & export")
# Join phyto data to 'correct' site-WB info & export####
df0_raw %>% 
  dplyr::select(
    site_id,
    biosys_code,
    biosys_short,
    sample_date,
    site_station_name,
    replicate_code,
    settling_time,
    sample_method,
    taxon_qualifier,
    colonies_per_litre_millilitre,
    analysis_type,
    analysis_method,
    date_of_analysis,
    surveillan,
    ngr,
    eastings,
    northings,
    taxon_rank_id,
    rank,
    parent_name_usage_id_x,
    parent_name_usage_id_y,
    kingdom,
    phylum,
    class,
    order,
    family,
    genus,
    is_marine,
    is_brackish,
    is_freshwater,
    is_terrestrial,
    is_extinct,
    carbon_value_match,
    valid_aphia_id,
    name_use,
    cells_per_litre,
    cells_per_m3,
    mean_vol_per_cell_um3,
    median_vol_per_cell_um3,
    mean_c_per_cell_pg_c,
    mean_c_per_cell_ug_c,
    mean_c_per_cell_mg_c,
    median_c_per_cell_pg_c,
    median_c_per_cell_ug_c,
    median_c_per_cell_mg_c,
    tot_mn_vol_um3,
    tot_md_vol_um3,
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
    md_carbon_tot_mg_c_per_m3,
    plankton_type,
    phytoplankton_type,
    phytoplankton_size,
    phyto_depth,
    phyto_feeding_mech,
    toxic_nuisance,
    phyto_habitat,
    protozoa_type,
    protozoa_size,
    protozoa_habitat,
    protozoa_feeding,
    phyto_lf
  ) %>% 
  dplyr::left_join(., joined_data_out, by = "site_id") -> df0_out

saveRDS(df0_out, file = "outputs/Phyto_2000_2025_with_Lifeforms_USE_updateCalcs_updateWB.Rdat")
write.csv(df0_out, file = "outputs/Phyto_2000_2025_with_Lifeforms_USE_updateCalcs_updateWB.csv",
          row.names = FALSE)
toc(log=TRUE)

unlist(tictoc::tic.log())

# Tidy up ####
rm(list=ls(pattern = "^cb"))
rm(list=ls(pattern = "^site"))
rm(df0_out,
   base_WBs,df0_raw,joined_data,joined_data_out, theme_use,
   GISfol,ppi,zoopfol, nearest_polygon)
