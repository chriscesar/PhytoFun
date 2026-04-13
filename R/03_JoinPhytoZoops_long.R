# 03_JoinPhytoZoops_long.R ####
# import the current zooplankton database and append for analysis
# This appends the 'shorter' zooplankton data (since 2022) to the complete phyto data
# (since 2000)

# load packages ####
ld_pkgs <- c("tidyverse","tictoc","stringr","sf","maps","nngeo")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

source("R/000setup.R")

#load data ####
tic("load data")
## ZOOPS
dfzoop <- read.csv((paste0(zoopfol,"processedData/zoopsAll.csv"))) %>% 
  # dfzoop <- read.csv("data/zoopsAll.csv") %>% 
  mutate(Biosys_short = if_else(
    is.na(BIOSYS.Code),
    NA_character_,
    substr(BIOSYS.Code, 1, nchar(BIOSYS.Code) - 1)
  )
  )
dfzoop$DataSet <- "Zooplankton"

zoop_meta <- readxl::read_xlsx(paste0(zoopfol,
                                      # "processedData/MBA_Returns_Amalgamated_USE.xlsx"),
                                      "processedData/260114_MBA_Returns_Amalgamated_USE.xlsx"),
                               sheet = "SiteMeta")

## phyto
dfphyto <- readRDS(file="outputs/Phyto_2000_2025_with_Lifeforms_USE.Rdat") %>% 
  # created in script '02_AppendLifeforms.R'
  mutate(Biosys_short = if_else(
    is.na(biosys_code),
    NA_character_,
    substr(biosys_code, 1, nchar(biosys_code) - 1)
  )
  )
dfphyto$DataSet <- "Phytoplankton"
toc(log=TRUE)

tic("Trim phyto by sites within zoop data & homogenise variable names")
# Homogenise variable names ####
dfphyto_trm <- dfphyto

## Prep zoop data ####
# retain 'sensible' variables
zooptrm <- dfzoop %>%
  dplyr::select(
    Pot.Number,
    sample.date,
    Biosys_short,
    DataSet,
    Aphia.ID,
    Taxa,
    Region,
    WBID,
    WB,
    Eastings,
    Northings,
    LF02,
    Kingdom, Phylum,Class,Order,Family,Genus,
    Abund_m3,
    mn_carbTot_m3_ug,
    md_carbTot_m3_ug,
    mnCPerIndiv_ug,
    mdCPerIndiv_ug,
    mnlongMaxAxis_mm,
    mdlongMaxAxis_mm
  ) %>% 
  ## create 'units' columns
  dplyr::mutate(
    Abund_m3_units = "Count_per_m3",
    mn_carbTot_units = "ugC_per_m3",
    md_carbTot_units = "ugC_per_m3",
    mnCPerIndiv_ug_units = "ugC_per_inividual",
    mdCPerIndiv_ug_units = "ugC_per_inividual",
    mnlongMaxAxis_mm_units = "mm",
    mdlongMaxAxis_mm_units = "mm"
  ) %>% 
  #rename for consistency
  dplyr::rename(
    sample_id = Pot.Number,
    region = Region,
    wb_id = WBID,
    wb = WB,
    eastings=Eastings,
    northings=Northings,
    sample_date = sample.date,
    aphia_id = Aphia.ID,
    taxon = Taxa,
    lifeform = LF02,
    abundance = Abund_m3,
    abundance_units = Abund_m3_units,
    mn_carbTot = mn_carbTot_m3_ug,
    md_carbTot = md_carbTot_m3_ug,
    mn_carbTot_units = mn_carbTot_units,
    md_carbTot_units = md_carbTot_units,
    mn_carbInd = mnCPerIndiv_ug,
    mn_carbInd_units = mnCPerIndiv_ug_units,
    md_carbInd = mdCPerIndiv_ug,
    md_carbInd_units = mdCPerIndiv_ug_units,
    mn_size = mnlongMaxAxis_mm,
    md_size = mdlongMaxAxis_mm,
    mn_size_units = mnlongMaxAxis_mm_units,
    md_size_units = mdlongMaxAxis_mm_units
  ) %>% #names()
  janitor::clean_names() %>% 
  ## reorder columns to match
  dplyr::select(
    sample_id,
    region,
    wb_id,
    wb,
    eastings,
    northings,
    biosys_short,
    sample_date,
    data_set,
    aphia_id,
    taxon,
    lifeform,
    kingdom,phylum,class,order,family,genus,
    abundance,abundance_units,
    mn_carb_tot,mn_carb_tot_units,
    md_carb_tot,md_carb_tot_units,
    mn_carb_ind,mn_carb_ind_units,
    md_carb_ind,md_carb_ind_units,
    mn_size,mn_size_units,
    md_size,md_size_units
  ) %>% 
  # convert all variables to character to allow easier joining
  mutate(across(everything(), as.character))

## Prep phyto data ####
# retain 'sensible' variables
phytotrm <- dfphyto_trm %>% dplyr::select(
  sample_id,
  sample_date,
  Biosys_short,
  DataSet,
  valid_aphia_id,
  name_use,
  river_basi,
  wb_id,
  wb_name,
  eastings,
  northings,
  phyto_lf,
  kingdom,phylum,class,order,family,genus,
  cells_per_litre_millilitre,
  tot_mn_c_pg_c_per_l,
  tot_md_c_pg_c_per_l,
  mean_c_per_cell_pg_c,
  median_c_per_cell_pg_c,
  mean_vol_per_cell_um3,
  median_vol_per_cell_um3,
  mean_c_per_cell_pg_c,
  median_c_per_cell_pg_c
) %>% 
  ## create 'units' columns
  dplyr::mutate(
    cells_per_litre_millilitre_units = "cells_per_litre",
    tot_mn_c_pg_c_per_l_units = "pg_C_per_litre",
    tot_md_c_pg_c_per_l_units = "pg_C_per_litre",
    mean_c_per_cell_pg_c_units = "pg_C_per_cell",
    median_c_per_cell_pg_c_units = "pg_C_per_cell",
    mean_vol_per_cell_um3_units = "individual_cell_volume_um3",
    median_vol_per_cell_um3_units = "individual_cell_volume_um3"
  ) %>% 
  # rename for consistency
  dplyr::rename(
    region = river_basi,
    wb_id = wb_id,
    wb = wb_name,
    biosys_short = Biosys_short,
    sample_date=sample_date,
    data_set = DataSet,
    aphia_id = valid_aphia_id,
    taxon = name_use,
    lifeform = phyto_lf,
    kingdom = kingdom, phylum = phylum, class = class,
    order = order, family = family, genus = genus,
    abundance = cells_per_litre_millilitre,
    abundance_units = cells_per_litre_millilitre_units,
    mn_carb_tot = tot_mn_c_pg_c_per_l,
    mn_carb_tot_units = tot_mn_c_pg_c_per_l_units,
    md_carb_tot = tot_md_c_pg_c_per_l,
    md_carb_tot_units = tot_md_c_pg_c_per_l_units,
    mn_carb_ind = mean_c_per_cell_pg_c,
    mn_carb_ind_units = mean_c_per_cell_pg_c_units,
    md_carb_ind = median_c_per_cell_pg_c,
    md_carb_ind_units = median_c_per_cell_pg_c_units,
    mn_size = mean_vol_per_cell_um3,
    mn_size_units = mean_vol_per_cell_um3_units,
    md_size = median_vol_per_cell_um3,
    md_size_units = median_vol_per_cell_um3_units
  ) %>% 
  janitor::clean_names() %>% 
  ## reorder columns to match
  dplyr::select(
    sample_id,region,wb_id,wb,
    eastings,
    northings,biosys_short,
    sample_date,
    data_set,
    aphia_id,
    taxon,
    lifeform,
    kingdom,phylum,class,order,family,genus,
    abundance,abundance_units,
    mn_carb_tot,mn_carb_tot_units,
    md_carb_tot,md_carb_tot_units,
    mn_carb_ind,mn_carb_ind_units,
    md_carb_ind,md_carb_ind_units,
    mn_size,mn_size_units,
    md_size,md_size_units
  ) %>% 
  # convert all variables to character to allow easier joining
  mutate(across(everything(), as.character))

## check names match
table(names(zooptrm) == names(phytotrm))
toc(log=TRUE)

tic("Bind zoops and phyto to single df")
# Bind zoops and phyto to single df and export ####
dfall <- bind_rows(zooptrm,phytotrm) %>% 
  mutate(sample_date = as.Date(sample_date),
         abundance = as.numeric(abundance),
         mn_carb_tot = as.numeric(mn_carb_tot),
         md_carb_tot = as.numeric(md_carb_tot),
         mn_size = as.numeric(mn_size),
         md_size = as.numeric(md_size)
  )

# homogenise region names and wb names ####
unique(dfall$region)

zoop_meta %>%
  mutate(biosys_short = substr(BIOSYS_ID, 1, nchar(BIOSYS_ID) - 1)) -> zoop_meta

dfall0 <- dfall

dfall %>% left_join(.,zoop_meta, by = "biosys_short") %>% 
  dplyr::select(-c(region:wb)) %>% #->d#%>% 
  dplyr::rename(
    region = RBID,
    wb_id = WBID,
    wb = WBName,
    wb_type = Type
  ) %>% 
  dplyr::relocate(wb_type) %>% 
  dplyr::relocate(wb) %>% 
  dplyr::relocate(wb_id) %>% 
  dplyr::relocate(region) %>% 
  dplyr::select(.,-c(WIMS:BIOSYS_ID)) %>% 
  mutate(abundance = as.numeric(abundance),
         mn_carb_tot = as.numeric(mn_carb_tot),
         md_carb_tot = as.numeric(md_carb_tot),
         mn_carb_ind = as.numeric(mn_carb_ind),
         md_carb_ind = as.numeric(md_carb_ind),
         mn_size = as.numeric(mn_size),
         md_size = as.numeric(md_size)
  ) %>% 
  dplyr::mutate(region = factor(region,
                                levels = c(
                                  "Northumbria","Humber",
                                  "Anglian","Thames",
                                  "South East","South West",
                                  "North West"
                                ))) -> dfall

toc(log = TRUE)  

tic("calculate values in consistent units")
# calculate values in consistent units ####
dfall %>% 
  ## create variable for micrograms of Carbon per INDIVIDUAL
  mutate(
    mn_carb_ind_as_ugC = if_else(mn_carb_ind_units == "ugC_per_inividual",
                                 as.numeric(mn_carb_ind),
                                 as.numeric(mn_carb_ind)/1000),
    md_carb_ind_as_ugC = if_else(md_carb_ind_units == "ugC_per_inividual",
                                 as.numeric(md_carb_ind),
                                 as.numeric(md_carb_ind)/1000)
  ) %>%
  ## create variable to show abundances per m3
  mutate(abundance_m3 = if_else(abundance_units == "Count_per_m3",
                                abundance,
                                abundance*1000)) %>% 
  ## multiply carbon per individual by abundance/m3
  ##to show carbon per m3
  mutate(mn_carb_ugC_per_m3 = abundance_m3*mn_carb_ind_as_ugC,
         md_carb_ugC_per_m3 = abundance_m3*md_carb_ind_as_ugC) %>% 
  mutate(eastings = as.numeric(eastings),
         northings = as.numeric(northings)) -> dfall
toc(log=TRUE)

## Ensure consistency in site WB names & regions ####
## extract unique grid refs ####
dfall %>%
  dplyr::select(eastings,northings) %>%
  mutate(EN=paste0(str_pad(eastings,width = 6, pad = "0"),
                   str_pad(northings,width = 6, pad = "0")
                   )) %>% distinct() %>% 
  as_tibble() -> sites

## join to WB ####

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
  result <- cbind(
    point,
    nearest_poly,
    distance = distance
    )
  return(result)
  }

tic("Join updated georef")
# Apply the function to each point
joined_data <- do.call(rbind, lapply(1:nrow(sites_sf), function(i) {
  nearest_polygon(sites_sf[i, ], base_WBs)
  }))

## join data ####
joined_data %>% 
  as_tibble() %>% 
  dplyr::select(EN,Country,WBID,WB_Name,WB_Type,Area_Name,
                Op_Catch,Typology,RBD,distance) %>% 
  mutate(RBD_lb = case_when(
    RBD == "South East" ~ "SE",
    RBD == "Northumbria" ~ "Nth",
    RBD == "Anglian" ~ "Ang",
    RBD == "Thames" ~ "Tha",
    RBD == "Humber" ~ "Hum",
    RBD == "North West" ~ "NW",
    RBD == "South West" ~ "SW",
    RBD == "Severn" ~ "Sev",
    RBD == "Solway Tweed" ~ "NW",
    TRUE ~ NA
    ),
    WB_lb = vegan::make.cepnames(WB_Name),
    RBD_WB_lb = paste0(RBD_lb,"_",WB_lb)
    ) -> sitedata

dfall %>% 
  mutate(EN=paste0(str_pad(eastings,width = 6, pad = "0"),
                   str_pad(northings,width = 6, pad = "0")
  )) %>% 
  left_join(., sitedata, by = "EN")%>% 
  dplyr::select(-c(region,wb_id,wb,wb_type,EN)) %>% 
  relocate(Country,Op_Catch,RBD,Area_Name,WBID,WB_Name,WB_Type,
           Typology,distance,RBD_lb,WB_lb,RBD_WB_lb) %>% 
  rename(dist_to_WB_m = distance) %>% 
  janitor::clean_names(.) -> dfall_join

toc(log=TRUE)

## investigate 'big' distance values
dfall_join %>% 
  filter(dist_to_wb_m >1000) %>% #View()
  ggplot(.,
         aes(y=dist_to_wb_m,
             x = rbd_wb_lb))+
  facet_wrap(.~rbd,scales = "free_x")+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(alpha = 0.05)+
  theme(
    strip.text = element_text(face=2)
  )

tic("Assign RBD factors")
dfall_join$rbd <- factor(dfall_join$rbd,levels = c(
  "Solway Tweed","Northumbria","Humber","Anglian",
  "Thames", "South East","South West","Severn",
  "North West"
))
toc(log=TRUE)

tic("Export data")
# Export data ####
saveRDS(dfall_join, file = "outputs/ZoopPhytoMatchingBIOSYS_long.Rdat")
write.csv(dfall_join, file = "outputs/ZoopPhytoMatchingBIOSYS_long.csv", row.names = FALSE)
toc(log=TRUE)

unlist(tictoc::tic.log())

# Tidy up ####
rm(list = ls(pattern = "^df"))
rm(list = ls(pattern = "^cb"))
rm(phytotrm, zooptrm, theme_use,ppi, zoopfol, zoop_meta,GISfol,
   base_WBs,joined_data,sitedata,sites,sites_sf,nearest_polygon)

detach("package:tidyverse", unload=TRUE)
detach("package:tidyr", unload=TRUE)
detach("package:tictoc", unload=TRUE)
detach("package:janitor", unload=TRUE)
detach("package:stringr", unload=TRUE)
