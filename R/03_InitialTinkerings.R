# 03_InitialTinkerings.R ####

# load packages ####
ld_pkgs <- c("tidyverse","tictoc","tsibble","feasts","zoo")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

source("R/000setup.R")

# Load data ####
tic("Load data")
# Phyto Time Series, generated in 02_AppendLifeforms.R
df <- readRDS(file="outputs/Phyto_2000_2025_with_Lifeforms_USE.Rdat")
toc(log=TRUE)

# Summarise total carbon/Abundance per sample ####
tic("Summarise total carbon/Abundance per sample")
tm <- df %>% dplyr::select(wb_name, sample_date) %>% distinct() %>% 
  arrange(wb_name,sample_date)

table(tm$wb_name) %>% as.data.frame() %>% arrange(-Freq)
# > 300 records:
# Mersey Mouth  647
# MERSEY  475
# Solent  441
# THAMES LOWER  420
# GREAT OUSE  418
# CARRICK ROADS INNER  394
# Wash Outer  394
# SOUTHAMPTON WATER  389
# POOLE HARBOUR  367
# THAMES MIDDLE  364
# TAW / TORRIDGE  362
# Northumberland North  345
# RIBBLE  339
# Morecambe Bay  315
# Whitstable Bay  311
# CAMEL  309
# LUNE  308
# TEES  302
rm(tm)

wb_dig <- c("TEES","LUNE","CAMEL","Whitstable Bay","Morecambe Bay",
            "RIBBLE","Mersey Mouth","MERSEY","Solent","THAMES LOWER",
            "GREAT OUSE","CARRICK ROADS INNER","Wash Outer","SOUTHAMPTON WATER",
            "POOLE HARBOUR","THAMES MIDDLE","TAW / TORRIDGE",
            "Northumberland North")

df %>% 
  # Keep only WBs with most data (n=300)
  #dplyr::filter(wb_name %in% wb_dig) %>% 
  # Remove NA values for tot_mn_c_pg_c
  dplyr::filter(!is.na(tot_mn_c_pg_c_per_l)) %>% 
  # Remove unwanted variables unique to taxa identity
  # aphia ID, name, lifeforms, taxonomy
  dplyr::select(-c(taxon_qualifier, colonies_per_litre_millilitre,
                   aphia_id, aphia_id_2,
                   scientificname, status, unacceptreason, taxon_rank_id,
                   rank, parent_name_usage_id, kingdom, phylum, class,
                   order, family, genus, is_marine, is_brackish,
                   is_freshwater, is_terrestrial, is_extinct,
                   size_class, qa_flag, plankton_type, phytoplankton_type,
                   phytoplankton_size, phyto_depth, phyto_feeding_mech,
                   toxic_nuisance, phyto_habitat, protozoa_type,
                   protozoa_size, protozoa_habitat, protozoa_feeding,
                   taxa_reported, valid_aphia_id, valid_name, name_use)) %>% 
  # Remove numerical variables not currently interested in
  # (retaining them will interfere with grouping)
  # dplyr::select(-c(mean_vol_per_cell_um3, median_vol_per_cell_um3,
  # mean_c_per_cell_pg_c,median_c_per_cell_pg_c)) %>% 
  group_by(across(!c(cells_per_litre_millilitre, tot_mn_vol_um3,
  tot_md_vol_um3, tot_mn_c_pg_c_per_l, tot_md_c_pg_c_per_l))) %>% 
  summarise(cells_per_litre_millilitre = sum(cells_per_litre_millilitre),
            tot_mn_vol_um3 = sum(tot_mn_vol_um3),
            tot_md_vol_um3 = sum(tot_md_vol_um3),
            tot_mn_c_pg_c_per_l = sum(tot_mn_c_pg_c_per_l),
            tot_md_c_pg_c_per_l = sum(tot_md_c_pg_c_per_l),
            mean_vol_per_cell_um3 = mean(mean_vol_per_cell_um3),
            median_vol_per_cell_um3=median(median_vol_per_cell_um3),
            mean_c_per_cell_pg_c=mean(mean_c_per_cell_pg_c),
            median_c_per_cell_pg_c=median(mean_c_per_cell_pg_c),
            .groups = "drop") %>% 
  ungroup() %>%
  arrange(river_basi,wb_name,site_id,sample_date) %>% 
  mutate(site_id = as.character(site_id)) -> df_carbsum

toc(log=TRUE)

tic("Initial plot of trends over time")

# df_carbsum %>% 
#   dplyr::filter(!is.na(biosys_code)) %>% 
#   #dplyr::filter(river_basi == "Northumbria") %>% 
#   #dplyr::filter(wb_name == "Northumberland North") %>% 
#   ggplot(aes(x = sample_date,
#              y = log(cells_per_litre_millilitre)
#                )
#          ) +
#   # facet_wrap(.~wb_name)+
#   facet_wrap(. ~ site_id)+
#   geom_point(aes(colour=site_id),show.legend = FALSE)+
#   theme_use

### only retain sites with N or more data points
num <- 50

# Step 1: Get sites with >num samples
sites_with_many_samples <- df_carbsum |>
  count(biosys_code) |>
  filter(n > num) %>% arrange(-n)

# Step 2: Filter and plot
df_carbsum |>
  semi_join(sites_with_many_samples, by = "biosys_code") |> #View()
  filter(!is.na(biosys_code)) |>
  filter(river_basi == "Anglian"|river_basi=="Humber") |>
  ggplot(aes(
    x = sample_date,
    # y = log(tot_mn_vol_um3)))+ #summed cell volume
    # y = log(cells_per_litre_millilitre))) + #summed cells per ml
    # y = log(tot_mn_c_pg_c)))+ #summed (log) carbon per sample
    # y = mean_c_per_cell_pg_c))+ #mean carbon per cell
    y = log(mean_vol_per_cell_um3)))+ #(log) mean cell volume
  facet_wrap(
    .~ water_body,
    # .~ biosys_code,
    scale = "free_y") +
  geom_point(aes(
    colour = site_id
    # colour = water_body
    ),
    alpha=0.05,
    show.legend = FALSE
    ) +
  # geom_smooth(method="gam",
              # aes(group = site_id, col = site_id),se=FALSE,
              # show.legend = FALSE, linewidth=2)+
  geom_smooth(method="gam", colour="red")+
  # geom_textsmooth(aes(label = biosys_code, group = biosys_code),
  #                 lty = 2,
  #                 hjust = 1,
  #                 method="gam"
  #                 )+
  theme_use

#############################################################################
#############################################################################
#############################################################################
#### TO SORT:
#### Samples with Biosys code as part of their name only go back as far as ~2007
#### This truncates the 

#### =================== plot vs Julian Day ===================
