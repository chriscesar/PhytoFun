# 02b_ConvertUnits_Phyto.R ####
# Load Phyto data and convert carbon values

# load packages ####
ld_pkgs <- c("tidyverse","tictoc","stringr")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

source("R/000setup.R")

# load data ####
tic("Load data")
dfphyto <- readRDS(file="outputs/Phyto_2000_2025_with_Lifeforms_USE.Rdat") %>% 
  # created in script '02_AppendLifeforms.R'
  mutate(Biosys_short = if_else(
    is.na(biosys_code),
    NA_character_,
    substr(biosys_code, 1, nchar(biosys_code) - 1)
  )
  )

toc(log=TRUE)

tic("Convert units")
# Convert units ####
# convert 'picograms carbon per litre' to 'micrograms carbon per litre'
# convert 'micrograms carbon per litre' to 'micrograms carbon per m3'
# Convert 'micrograms carbon per m3' to 'milligrams carbon per m3'

dfphyto %>% 
  ### mean_c_per_cell_pg_c
  # convert pg to ug (multiply by 1e-6)
  dplyr::mutate(mean_c_per_cell_ug_c = mean_c_per_cell_pg_c*1e-6) %>% 
  dplyr::relocate(mean_c_per_cell_ug_c, .after = mean_c_per_cell_pg_c) %>% 
  
  # convert um to mg (multiply by 0.001)
  dplyr::mutate(mean_c_per_cell_mg_c = mean_c_per_cell_ug_c*0.001) %>% 
  dplyr::relocate(mean_c_per_cell_mg_c, .after = mean_c_per_cell_ug_c) %>%
  
  ### median_c_per_cell_pg_c
  # convert pg to ug (multiply by 1e-6)
  dplyr::mutate(median_c_per_cell_ug_c = median_c_per_cell_pg_c*1e-6) %>% 
  dplyr::relocate(median_c_per_cell_ug_c, .after = median_c_per_cell_pg_c) %>% 
  
  # convert um to mg (multiply by 0.001)
  dplyr::mutate(median_c_per_cell_mg_c = median_c_per_cell_ug_c*0.001) %>% 
  dplyr::relocate(median_c_per_cell_mg_c, .after = median_c_per_cell_ug_c) %>%
  
  ### Convert cells_per_litre_millilitre to cells per m3
  # convert cells per litre to cells per m3 (multiply value/litre by 1000)
  dplyr::mutate(cells_per_m3 = cells_per_litre_millilitre*1000) %>% 
  dplyr::relocate(cells_per_m3, .after = cells_per_litre_millilitre) %>% 
  
  ### Convert carbon per cell to carbon per litre
  ## MEAN values
  # Picograms carbon per litre (multiply cells per litre by pg carbon per cell)
  dplyr::mutate(mn_carbon_tot_pg_C_per_litre = cells_per_litre_millilitre * mean_c_per_cell_pg_c) %>% 
  dplyr::relocate(mn_carbon_tot_pg_C_per_litre, .after = tot_mn_c_pg_c_per_l) %>% 
  
  # Micrograms carbon per litre (multiply cells per litre by ug carbon per cell)
  dplyr::mutate(mn_carbon_tot_ug_C_per_litre = cells_per_litre_millilitre * mean_c_per_cell_ug_c) %>% 
  dplyr::relocate(mn_carbon_tot_ug_C_per_litre, .after = mn_carbon_tot_pg_C_per_litre) %>%
  
  # Milligrams carbon per litre (multiply cells per litre by mg carbon per cell)
  dplyr::mutate(mn_carbon_tot_mg_C_per_litre = cells_per_litre_millilitre*mean_c_per_cell_mg_c) %>% 
  dplyr::relocate(mn_carbon_tot_mg_C_per_litre, .after = mn_carbon_tot_ug_C_per_litre) %>% 
  
  ### Convert carbon per litre to carbon per m3
  # picograms per litre to picograms per m3 (multiply value/litre by 1000)
  dplyr::mutate(mn_carbon_tot_pg_C_per_m3 = mn_carbon_tot_pg_C_per_litre*1000) %>% 
  dplyr::relocate(mn_carbon_tot_pg_C_per_m3, .after = mn_carbon_tot_mg_C_per_litre) %>% 
  
  # micrograms per litre to micrograms per m3 (multiply value/litre by 1000)
  dplyr::mutate(mn_carbon_tot_ug_C_per_m3 = mn_carbon_tot_ug_C_per_litre*1000) %>%
  dplyr::relocate(mn_carbon_tot_ug_C_per_m3,.after = mn_carbon_tot_pg_C_per_m3) %>% 
  
  # milligrams per litre to milligrams per m3 (multiply value/litre by 1000)
  dplyr::mutate(mn_carbon_tot_mg_C_per_m3 = mn_carbon_tot_mg_C_per_litre*1000) %>% 
  dplyr::relocate(mn_carbon_tot_mg_C_per_m3, .after = mn_carbon_tot_ug_C_per_m3) %>% 
  
  ## MEDIAN values
  # Picograms carbon per litre (multiply cells per litre by pg carbon per cell)
  dplyr::mutate(md_carbon_tot_pg_C_per_litre = cells_per_litre_millilitre * median_c_per_cell_pg_c) %>% 
  dplyr::relocate(md_carbon_tot_pg_C_per_litre, .after = tot_md_c_pg_c_per_l) %>% 
  
  # Micrograms carbon per litre (multiply cells per litre by ug carbon per cell)
  dplyr::mutate(md_carbon_tot_ug_C_per_litre = cells_per_litre_millilitre * median_c_per_cell_ug_c) %>% 
  dplyr::relocate(md_carbon_tot_ug_C_per_litre, .after = md_carbon_tot_pg_C_per_litre) %>%
  
  # Milligrams carbon per litre (multiply cells per litre by mg carbon per cell)
  dplyr::mutate(md_carbon_tot_mg_C_per_litre = cells_per_litre_millilitre*median_c_per_cell_mg_c) %>% 
  dplyr::relocate(md_carbon_tot_mg_C_per_litre, .after = md_carbon_tot_ug_C_per_litre) %>% 
  
  ### Convert carbon per litre to carbon per m3
  # picograms per litre to picograms per m3 (multiply value/litre by 1000)
  dplyr::mutate(md_carbon_tot_pg_C_per_m3 = md_carbon_tot_pg_C_per_litre*1000) %>% 
  dplyr::relocate(md_carbon_tot_pg_C_per_m3, .after = md_carbon_tot_mg_C_per_litre) %>% 
  # micrograms per litre to micrograms per m3 (multiply value/litre by 1000)
  dplyr::mutate(md_carbon_tot_ug_C_per_m3 = md_carbon_tot_ug_C_per_litre*1000) %>%
  dplyr::relocate(md_carbon_tot_ug_C_per_m3,.after = md_carbon_tot_pg_C_per_m3) %>% 
  # milligrams per litre to milligrams per m3 (multiply value/litre by 1000)
  dplyr::mutate(md_carbon_tot_mg_C_per_m3 = md_carbon_tot_mg_C_per_litre*1000) %>% 
  dplyr::relocate(md_carbon_tot_mg_C_per_m3, .after = md_carbon_tot_ug_C_per_m3) %>% 
  dplyr::select(-c(tot_mn_c_pg_c_per_l,tot_md_c_pg_c_per_l)) %>% 
  dplyr::rename(cells_per_litre = cells_per_litre_millilitre) -> df0

write.csv(df0,
          file = "outputs/Phyto_2000_2025_with_Lifeforms_USE_updateCalcs.csv",
          row.names = FALSE)
saveRDS(df0,
        file = "outputs/Phyto_2000_2025_with_Lifeforms_USE_updateCalcs.Rdat"
        )
toc(log = TRUE)

unlist(tictoc::tic.log())

# Tidy up ####
rm(list = ls(pattern = "^df"))
rm(list = ls(pattern = "^cb"))
rm(theme_use, ppi,zoopfol)

detach("package:tidyverse", unload=TRUE)
detach("package:janitor", unload=TRUE)
detach("package:stringr", unload=TRUE)
detach("package:tictoc", unload=TRUE)
