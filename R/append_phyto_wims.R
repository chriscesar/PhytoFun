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

########################
## Metadata ####

dfmeta <- dftxtrm %>% 
  dplyr::select(
    code,
    river_basi,
    wb_id,
    wb_name,
    site_id,
    biosys_code,
    sample_date,
    sample_id,
    wb_type,
    wb_hm_designation,
    wb_typology,
    surveillan
  ) %>% 
  dplyr::rename(sample_id_biosys = sample_id) %>% 
  distinct()

########################
## Abundance data versions ####

df_abund_taxa_per_l <- dftxtrm %>% 
  dplyr::rename(sample_id_biosys = sample_id) %>% 
  dplyr::select(
    code,
    sample_id_biosys,
    name_use,
    cells_per_litre_millilitre
  ) %>% 
  dplyr::group_by(across(-cells_per_litre_millilitre)) %>% 
  dplyr::summarise(cells_per_litre_millilitre=sum(cells_per_litre_millilitre),
                   .groups = "drop") %>% 
  tidyr::pivot_wider(
    names_from = name_use,
    values_from = cells_per_litre_millilitre
  )

df_abund_lf <- dftxtrm %>% 
  dplyr::rename(sample_id_biosys = sample_id) %>% 
  dplyr::select(
    code,
    sample_id_biosys,
    phyto_lf,
    cells_per_litre_millilitre
  ) %>% 
  dplyr::group_by(across(-cells_per_litre_millilitre)) %>% 
  dplyr::summarise(cells_per_litre_millilitre=sum(cells_per_litre_millilitre),
                   .groups = "drop") %>% 
  tidyr::pivot_wider(
    names_from = phyto_lf,
    values_from = cells_per_litre_millilitre
  )

df_abund_lf2 <- dftxtrm %>% 
  dplyr::rename(sample_id_biosys = sample_id) %>% 
  dplyr::select(
    code,
    sample_id_biosys,
    phyto_lf2,
    cells_per_litre_millilitre
  ) %>% 
  dplyr::group_by(across(-cells_per_litre_millilitre)) %>% 
  dplyr::summarise(cells_per_litre_millilitre=sum(cells_per_litre_millilitre),
                   .groups = "drop") %>% 
  tidyr::pivot_wider(
    names_from = phyto_lf2,
    values_from = cells_per_litre_millilitre
  )
########################
## Carbon data versions ####
df_carb_taxa_md_pg_per_l <- dftxtrm %>% 
  dplyr::rename(sample_id_biosys = sample_id) %>% 
  dplyr::select(
    code,
    sample_id_biosys,
    name_use,
    tot_md_c_pg_c_per_l
  ) %>% 
  dplyr::group_by(across(-tot_md_c_pg_c_per_l)) %>% 
  dplyr::summarise(tot_md_c_pg_c_per_l=sum(tot_md_c_pg_c_per_l),
                   .groups = "drop") %>% 
  tidyr::pivot_wider(
    names_from = name_use,
    values_from = tot_md_c_pg_c_per_l
  )

df_carb_lf_md_pg_per_l <- dftxtrm %>% 
  dplyr::rename(sample_id_biosys = sample_id) %>% 
  dplyr::select(
    code,
    sample_id_biosys,
    phyto_lf,
    tot_md_c_pg_c_per_l
  ) %>% 
  dplyr::group_by(across(-tot_md_c_pg_c_per_l)) %>% 
  dplyr::summarise(tot_md_c_pg_c_per_l=sum(tot_md_c_pg_c_per_l),
                   .groups = "drop") %>% 
  tidyr::pivot_wider(
    names_from = phyto_lf,
    values_from = tot_md_c_pg_c_per_l
  )

df_carb_lf2_md_pg_per_l <- dftxtrm %>% 
  dplyr::rename(sample_id_biosys = sample_id) %>% 
  dplyr::select(
    code,
    sample_id_biosys,
    phyto_lf2,
    tot_md_c_pg_c_per_l
  ) %>% 
  dplyr::group_by(across(-tot_md_c_pg_c_per_l)) %>% 
  dplyr::summarise(tot_md_c_pg_c_per_l=sum(tot_md_c_pg_c_per_l),
                   .groups = "drop") %>% 
  tidyr::pivot_wider(
    names_from = phyto_lf2,
    values_from = tot_md_c_pg_c_per_l
  )

################
## WIMS data ####
tictoc::tic("Reformat WIMS data")
wimsdat_wb_trm %>% 
  dplyr::mutate(
    det_units = paste0(DETE_SHORT_DESC,
                       "_",
                       UNIT_SHORT_DESC)) %>%
  dplyr::mutate(
    det_units = snakecase::to_snake_case(det_units)
  ) %>% 
  ## remove ODD samples which flag as duplicates
  dplyr::filter(SAMP_ID != "1568725" & det_units != "123789_hx_cdd_pg_l") %>% 
  dplyr::filter(SAMP_ID != "1563520" & det_units != "123789_hx_cdd_pg_l") %>% 
  dplyr::filter(SAMP_ID != "1564003" & det_units != "123789_hx_cdd_pg_l") %>% 
  dplyr::filter(SAMP_ID != "1568839" & det_units != "123789_hx_cdd_pg_l") %>% 
  dplyr::filter(SAMP_ID != "5196058" & det_units != "123789_hx_cdd_pg_l") %>% 
  dplyr::filter(SAMP_ID != "5198151" & det_units != "123789_hx_cdd_pg_l") %>% 
  dplyr::filter(SAMP_ID != "926030" & det_units != "orthophs_filt_mg_l") %>% 
  dplyr::filter(SAMP_ID != "926033" & det_units != "orthophs_filt_mg_l") %>% 
  dplyr::filter(SAMP_ID != "926038" & det_units != "orthophs_filt_mg_l") %>% 
  dplyr::filter(SAMP_ID != "938194" & det_units != "1_3_dichlor_ug_l") %>% 
  dplyr::filter(SAMP_ID != "979790" & det_units != "orthophs_filt_mg_l") %>% 
  dplyr::filter(SAMP_ID != "983856" & det_units != "demeton_s_ug_l") %>% 
  dplyr::filter(SAMP_ID != "979785" & det_units != "orthophs_filt_mg_l") %>% 
  dplyr::filter(SAMP_ID != "3295414" & det_units != "123789_hx_cdd_pg_l") %>% 
  dplyr::filter(SAMP_ID != "3295127" & det_units != "123789_hx_cdd_pg_l") %>% 
  dplyr::filter(SAMP_ID != "1563813" & det_units != "123789_hx_cdd_pg_l") %>% 
  dplyr::filter(SAMP_ID != "1568958" & det_units != "123789_hx_cdd_pg_l") %>% 
  dplyr::filter(SAMP_ID != "1054911" & det_units != "demeton_s_ug_l") %>% 
  
  ## sometimes get multiple samples within a WB in a day.  Calculate mean values
  ## across wbs for each date
  dplyr::group_by(across(!MEAS_RESULT)) %>% 
  
  
  
  dplyr::mutate(
    result = if_else(
      is.na(MEAS_SIGN),
      as.character(MEAS_RESULT),
      paste(MEAS_SIGN, MEAS_RESULT)
      )
    ) %>%
  # dplyr::select(-c(DETE_SHORT_DESC,
  #                  UNIT_SHORT_DESC,
  #                  SAMP_CONFIDENTIAL,
  #                  SAMP_MATERIAL,
  #                  DATE_TIME,
  #                  RIVER_BASI,
  #                  notation,
  #                  MEAS_SIGN,
  #                  MEAS_RESULT,
  #                  MEAS_ANAL_METH_CODE,
  #                  MEAS_DETERMINAND_CODE,
  #                  MEAS_LIMITS,
  #                  res,
  #                  wims_region,
  #                  time_period,
  #                  material_type,
  #                  sample_date,
  #                  WBID,
  #                  WBName,
  #                  Type,
  #                  EA_AREA_NA
  #                  )
  #               ) %>%
  dplyr::select(
    SAMP_SMPT_USER_REFERENCE,
    SAMP_ID,
    #SAMP_PURPOSE_CODE,
    code,
    det_units,
    result
  ) %>% 
  # ## remove ODD samples which flag as duplicates
  # dplyr::filter(SAMP_ID != "1568725" & det_units != "123789_hx_cdd_pg_l") %>% 
  # dplyr::filter(SAMP_ID != "1563520" & det_units != "123789_hx_cdd_pg_l") %>% 
  # dplyr::filter(SAMP_ID != "1564003" & det_units != "123789_hx_cdd_pg_l") %>% 
  # dplyr::filter(SAMP_ID != "1568839" & det_units != "123789_hx_cdd_pg_l") %>% 
  # dplyr::filter(SAMP_ID != "5196058" & det_units != "123789_hx_cdd_pg_l") %>% 
  # dplyr::filter(SAMP_ID != "5198151" & det_units != "123789_hx_cdd_pg_l") %>% 
  # dplyr::filter(SAMP_ID != "926030" & det_units != "orthophs_filt_mg_l") %>% 
  # dplyr::filter(SAMP_ID != "926033" & det_units != "orthophs_filt_mg_l") %>% 
  # dplyr::filter(SAMP_ID != "926038" & det_units != "orthophs_filt_mg_l") %>% 
  # dplyr::filter(SAMP_ID != "938194" & det_units != "1_3_dichlor_ug_l") %>% 
  # dplyr::filter(SAMP_ID != "979790" & det_units != "orthophs_filt_mg_l") %>% 
  # dplyr::filter(SAMP_ID != "983856" & det_units != "demeton_s_ug_l") %>% 
  # dplyr::filter(SAMP_ID != "979785" & det_units != "orthophs_filt_mg_l") %>% 
  # dplyr::filter(SAMP_ID != "3295414" & det_units != "123789_hx_cdd_pg_l") %>% 
  # dplyr::filter(SAMP_ID != "3295127" & det_units != "123789_hx_cdd_pg_l") %>% 
  # dplyr::filter(SAMP_ID != "1563813" & det_units != "123789_hx_cdd_pg_l") %>% 
  # dplyr::filter(SAMP_ID != "1568958" & det_units != "123789_hx_cdd_pg_l") %>% 
  # dplyr::filter(SAMP_ID != "1054911" & det_units != "demeton_s_ug_l") %>% 
  
  
  tidyr::pivot_wider(
    names_from = det_units,
    values_from = result
  ) -> dfwims

toc(log=TRUE)

# dfwims_tmp|>
#   dplyr::summarise(n = dplyr::n(), .by = c(SAMP_SMPT_USER_REFERENCE, SAMP_ID, SAMP_PURPOSE_CODE, code,
#                                            det_units)) |>
#   dplyr::filter(n > 1L) 
