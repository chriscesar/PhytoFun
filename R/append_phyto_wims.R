# append_phyto_wims.R ####
## append wims and phyto taxon/lifeform data

# load packages ####
ld_pkgs <- c("tidyverse","tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

tic("load phyto data")
# load phyto data ####
dftx <- readRDS(file="outputs/Phyto_2000_2025_with_Lifeforms_USE.Rdat")
toc(log=TRUE)

tic("load WIMS data")
source("R/loadWIMS.R")
toc(log=TRUE)

# create 'code' variable for phyto data ####
dftx %>% 
  dplyr::mutate(code = paste0(wb_id,"_",sample_date)) %>% 
  dplyr::mutate(code = str_remove_all(code, "-")) -> dftx

# flag which phyto samples have linked wims data & retain only those
x <- dftx$code %in% wimsdat_wb$code
dftxtrm <- dftx[x,]
rm(x)

# flag which WIMS samples have linked phyto data & retain only those
x <-  wimsdat_wb$code %in% dftx$code
wimsdat_wb_trm <- wimsdat_wb[x,]
rm(x)

# trim out unneccessary columns ####
## phyto ####
#names(dftxtrm)
dftxtrm %>% 
  dplyr::select(
    river_basi,
    wb_name,
    site_id,
    biosys_code,
    sample_date,
    # site_station_name,
    # ngr,
    # region,
    # area,
    sample_id,
    # sample_time,
    # sample_date_time,
    # sample_reason_primary,
    # analysis_id,
    # water_body,
    # water_body_type,
    # replicate_code,
    # settling_time,
    # taxon_qualifier,
    # colonies_per_litre_millilitre,
    # sample_ngr,
    # sample_method,
    # wfd_waterbody_id,
    # catchment,
    # version,
    # analysis_type,
    # analysis_method,
    # date_of_analysis,
    wb_id,
    wb_type,
    wb_hm_designation,
    wb_typology,
    surveillan,
    # distance_to_wb,
    # eastings,
    # northings,
    # aphia_id,
    # aphia_id_2,
    # scientificname,
    # status,
    # unacceptreason,
    # taxon_rank_id,                
    # rank,
    # parent_name_usage_id_x,
    # kingdom,
    # phylum,
    # class,
    # order,
    # family,
    # genus,
    # is_marine,
    # is_brackish,
    # is_freshwater,
    # is_terrestrial,
    # is_extinct,
    # x2,
    # parent_name_usage_id_y,
    carbon_value_match,
    # taxa_reported,
    valid_aphia_id,
    # valid_name,
    name_use,
    cells_per_litre_millilitre,   
    mean_vol_per_cell_um3,
    median_vol_per_cell_um3,
    mean_c_per_cell_pg_c,
    median_c_per_cell_pg_c,
    tot_mn_vol_um3,
    tot_md_vol_um3,
    tot_mn_c_pg_c_per_l,
    tot_md_c_pg_c_per_l,
    # size_class,
    # qa_flag,
    # plankton_type,
    # phytoplankton_type,
    phytoplankton_size,
    # phyto_depth,
    # phyto_feeding_mech,
    # toxic_nuisance,
    # phyto_habitat,
    # protozoa_type,
    protozoa_size,
    # protozoa_habitat,
    # protozoa_feeding,
    phyto_lf,
    code
  ) %>% 
  # filter animals from phyto data
  dplyr::filter(name_use != "Lyonsia") %>% 
  dplyr::filter(name_use != "Strombus") -> xx

tmp_phyto <- xx %>% dplyr::filter(phyto_lf %in% c(
  "Diatom","Dinoflagellate","Haptophyte","Cyanobacteria","Chrysophyte",
  "Chlorophyte","Silicoflagellate","Charophyte","Dictyochophyte","Raphidophyte",
  "Cryptophyte","Phytoplankton","Xanthophyte"
  )) %>% 
  dplyr::mutate(phyto_lf2 = paste0(phyto_lf,"_",phytoplankton_size))

tmp_prot <- xx %>% dplyr::filter(phyto_lf %in% c(
  "Protozoa","Ciliate"
  )) %>% 
  dplyr::mutate(phyto_lf2 = paste0(phyto_lf,"_",protozoa_size))  

xx_na <- xx %>% filter(is.na(phyto_lf)) %>% 
  dplyr::mutate(
    phyto_lf = "Other phytoplankton",
    phyto_lf2 = "Other phytoplankton")

dftxtrm <- rbind(tmp_phyto, tmp_prot, xx_na)

rm(xx,xx_na,tmp_phyto,tmp_prot)
dftxtrm <- dftxtrm %>% 
  dplyr::select(
    -c(phytoplankton_size,protozoa_size)
    )

## pull out metadata, abundance data and carbon data into separate list elements
## and widen
## do this for both taxa and lifeforms

