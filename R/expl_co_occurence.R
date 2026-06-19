# expl_co_occurence.R ----

## Exploration of how  selected phytoplankters co-occur
## DATE: 18/06/2026
## DATE Last Updated: 
## AUTHOR: Dr Christopher Cesar
## EDITORS:

# load packages & set meta ----
ld_pkgs <- c("tidyverse","tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

tictoc::tic("Load data and set metadata")
## set meta inputs
source("R/000setup.R")
# set ouput folder to save files
output.location <- "//prodds.ntnl/Shared/AN/KFH/Groups/N_Marine/02 Projects_Tasks/04 E&B/Ecology and Ecosystems/Planktonic_Food_Webs/Co_occurrence/"

## load data
df0 <- readRDS("outputs/ZoopPhyto_long_USE.Rdat") %>% 
  # convert date from character to date
  dplyr::mutate(sample_date = lubridate::date(sample_date)) %>% 
  # create year-day variable
  dplyr::mutate(yday = lubridate::yday(sample_date)) %>% 
  # create season variable
  dplyr::mutate(seas = seas::mkseas(sample_date, width = "DJF"))
tictoc::toc(log = TRUE)

#check sample_id 607949
# generate output of unique taxa in data ####
# names(df0)
# df0 %>%
#   dplyr::select(
#     data_set,
#     taxon,
#     aphia_id,
#     lifeform,
#     kingdom,
#     subkingdom,
#     infrakingdom,
#     phylum,
#     phylum_division,
#     subphylum,
#     subphylum_subdivision,
#     infraphylum,
#     parvphylum,
#     gigaclass,
#     superclass,
#     class,
#     subclass,
#     infraclass,
#     subterclass,
#     superorder,
#     order,
#     suborder,
#     infraorder,
#     parvorder,
#     superfamily,
#     family,
#     subfamily,
#     tribe,
#     genus,
#     subgenus,
#     species,
#     forma,
#     variety,
#     section,
#     subsection,
#     subspecies
#     ) %>%
#   dplyr::count(dplyr::across(dplyr::everything()),
#                name = "n_occurrences") %>% 
#   # distinct() #%>%
#   write.csv(.,
#             file = "outputs/unique_taxa.csv",
#             row.names = FALSE)

# Flag taxa of interest ----
tictoc::tic("Flag & generate data for taxa of interest")

# create new version of file
df_all <- df0

# Phaeocystis: Genus Phaeocystis
df_all %>% 
  dplyr::filter(!is.na(genus) & genus == "Phaeocystis") %>% 
  pull(.,aphia_id) %>% unique() -> aphia_phaeocystis

df_all %>% dplyr::mutate(
  phaeocystis_abund_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_phaeocystis ~ abundance_m3,
    TRUE ~ NA_real_
    ),
  phaeocystis_ugC_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_phaeocystis ~ md_carb_ugC_per_m3,
    TRUE ~ NA_real_
    )
  )-> df_all

## Coccolithophores: Class Coccolithophyceae
### (this also includes Genus Phaeocystis)
df_all %>% 
  dplyr::filter(!is.na(class) & class == "Coccolithophyceae") %>% 
  pull(.,aphia_id) %>% unique() -> aphia_coccolith

df_all %>% dplyr::mutate(
  coccolithophyceae_abund_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_coccolith ~ abundance_m3,
    TRUE ~ NA_real_
  ),
  coccolithophyceae_ugC_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_coccolith ~ md_carb_ugC_per_m3,
    TRUE ~ NA_real_
    )
  ) -> df_all

## Chaetoceros: Genus Chaetoceros
df_all %>% 
  dplyr::filter(!is.na(genus) & genus == "Chaetoceros") %>% 
  pull(.,aphia_id) %>% unique() -> aphia_chaetoceros

df_all %>% dplyr::mutate(
  chaetoceros_abund_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_chaetoceros ~ abundance_m3,
    TRUE ~ NA_real_
  ),
  chaetoceros_ugC_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_chaetoceros ~ md_carb_ugC_per_m3,
    TRUE ~ NA_real_
    )
  ) -> df_all

## Skeletonema: Genus Skeletonema
df_all %>% 
  dplyr::filter(!is.na(genus) & genus == "Skeletonema") %>% 
  pull(.,aphia_id) %>% unique() -> aphia_skeletonema

df_all %>% dplyr::mutate(
  skeletonema_abund_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_skeletonema ~ abundance_m3,
    TRUE ~ NA_real_
  ),
  skeletonema_ugC_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_skeletonema ~ md_carb_ugC_per_m3,
    TRUE ~ NA_real_
    )
  ) -> df_all

## Thalassiosira: Genus Thalassiosira
df_all %>% 
  dplyr::filter(!is.na(genus) & genus == "Thalassiosira") %>% 
  pull(.,aphia_id) %>% unique() -> aphia_thalassiosira

df_all %>% dplyr::mutate(
  thalassiosira_abund_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_thalassiosira ~ abundance_m3,
    TRUE ~ NA_real_
  ),
  thalassiosira_ugC_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_thalassiosira ~ md_carb_ugC_per_m3,
    TRUE ~ NA_real_
    )
  ) -> df_all

## Cerataulina: Genus Cerataulina
df_all %>% 
  dplyr::filter(!is.na(genus) & genus == "Cerataulina") %>% 
  pull(.,aphia_id) %>% unique() -> aphia_cerataulina

df_all %>% dplyr::mutate(
  cerataulina_abund_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_cerataulina ~ abundance_m3,
    TRUE ~ NA_real_
  ),
  cerataulina_ugC_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_cerataulina ~ md_carb_ugC_per_m3,
    TRUE ~ NA_real_
    )
  ) -> df_all

## Diatoms: Class Bacillariophyceae
### Contains ~200 taxa
df_all %>% 
  dplyr::filter(!is.na(class) & class == "Bacillariophyceae") %>% 
  pull(.,aphia_id) %>% unique() -> aphia_diatoms

df_all %>% dplyr::mutate(
  diatoms_abund_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_diatoms ~ abundance_m3,
    TRUE ~ NA_real_
  ),
  diatoms_ugC_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_diatoms ~ md_carb_ugC_per_m3,
    TRUE ~ NA_real_
    )
  ) -> df_all

## Odontella: Genus Odontella
df_all %>% 
  dplyr::filter(!is.na(genus) & genus == "Odontella") %>% 
  pull(.,aphia_id) %>% unique() -> aphia_odontella

df_all %>% dplyr::mutate(
  odontella_abund_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_odontella ~ abundance_m3,
    TRUE ~ NA_real_
  ),
  odontella_ugC_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_odontella ~ md_carb_ugC_per_m3,
    TRUE ~ NA_real_
    )
  ) -> df_all

## Coscinodiscus: Genus Coscinodiscus
df_all %>% 
  dplyr::filter(!is.na(genus) & genus == "Coscinodiscus") %>% 
  pull(.,aphia_id) %>% unique() -> aphia_coscinodiscus

df_all %>% dplyr::mutate(
  coscinodiscus_abund_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_coscinodiscus ~ abundance_m3,
    TRUE ~ NA_real_
  ),
  coscinodiscus_ugC_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_coscinodiscus ~ md_carb_ugC_per_m3,
    TRUE ~ NA_real_
    )
  ) -> df_all

## Actinocyclus: Genus Actinocyclus OR Genus Actinoptychus
df_all %>% 
  dplyr::filter(!is.na(genus) & genus %in% c("Actinocyclus", "Actinoptychus")) %>% 
  pull(.,aphia_id) %>% unique() -> aphia_actinocyclus_actinoptychus

df_all %>% dplyr::mutate(
  actinocyclus_actinoptychus_abund_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_actinocyclus_actinoptychus ~ abundance_m3,
    TRUE ~ NA_real_
  ),
  actinocyclus_actinoptychus_ugC_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_actinocyclus_actinoptychus ~ md_carb_ugC_per_m3,
    TRUE ~ NA_real_
    )
  ) -> df_all

## Dinoflagellates: Infraphylum Dinoflagellata
df_all %>% 
  dplyr::filter(!is.na(infraphylum) & infraphylum == "Dinoflagellata") %>% 
  pull(.,aphia_id) %>% unique() -> aphia_dinoflagellata

df_all %>% dplyr::mutate(
  dinoflagellata_abund_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_dinoflagellata ~ abundance_m3,
    TRUE ~ NA_real_
  ),
  dinoflagellata_ugC_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_dinoflagellata ~ md_carb_ugC_per_m3,
    TRUE ~ NA_real_
    )
  ) -> df_all

## Pseudo-nitzschia: Genus Pseudo-nitzschia
df_all %>% 
  dplyr::filter(!is.na(genus) & genus == "Pseudo-nitzschia") %>% 
  pull(.,aphia_id) %>% unique() -> aphia_pseudo_nitzschia

df_all %>% dplyr::mutate(
  pseudo_nitzschia_abund_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_pseudo_nitzschia ~ abundance_m3,
    TRUE ~ NA_real_
  ),
  pseudo_nitzschia_ugC_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_pseudo_nitzschia ~ md_carb_ugC_per_m3,
    TRUE ~ NA_real_
    )
  ) -> df_all

## Karenia: Genus Karenia
df_all %>% 
  dplyr::filter(!is.na(genus) & genus == "Karenia") %>% 
  pull(.,aphia_id) %>% unique() -> aphia_karenia

df_all %>% dplyr::mutate(
  karenia_abund_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_karenia ~ abundance_m3,
    TRUE ~ NA_real_
  ),
  karenia_ugC_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_karenia ~ md_carb_ugC_per_m3,
    TRUE ~ NA_real_
    )
  ) -> df_all

## Dinophysis: Genus Dinophysis
df_all %>% 
  dplyr::filter(!is.na(genus) & genus == "Dinophysis") %>% 
  pull(.,aphia_id) %>% unique() -> aphia_dinophysis

df_all %>% dplyr::mutate(
  dinophysis_abund_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_dinophysis ~ abundance_m3,
    TRUE ~ NA_real_
  ),
  dinophysis_ugC_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_dinophysis ~ md_carb_ugC_per_m3,
    TRUE ~ NA_real_
    )
  ) -> df_all

## Noctiluca: Order Noctilucales
df_all %>% 
  dplyr::filter(!is.na(order) & order == "Noctilucales") %>% 
  pull(.,aphia_id) %>% unique() -> aphia_noctilucales

df_all %>% dplyr::mutate(
  noctilucales_abund_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_noctilucales ~ abundance_m3,
    TRUE ~ NA_real_
  ),
  noctilucales_ugC_m3 = case_when(
    !is.na(aphia_id) & aphia_id %in% aphia_noctilucales ~ md_carb_ugC_per_m3,
    TRUE ~ NA_real_
  )
) -> df_all

# Tidy up
rm(aphia_phaeocystis,aphia_coccolith,aphia_chaetoceros,
   aphia_actinocyclus_actinoptychus, aphia_skeletonema,aphia_thalassiosira,
   aphia_cerataulina,aphia_diatoms,aphia_odontella,aphia_coscinodiscus,
   aphia_dinoflagellata,aphia_pseudo_nitzschia,
   aphia_karenia,aphia_dinophysis,aphia_noctilucales)

tictoc::toc(log = TRUE)

# tidy up data & convert to long ----
tictoc::tic("Convert to long")
df_all %>% #names() %>% 
  ## trim data into variables of interest
  dplyr::select(
    sample_id,
    region,
    region_lab,
    wb_id, wb,
    sample_date,
    yday,
    seas,
    # abundance_m3,
    # md_carb_ugC_per_m3,
    62:dplyr::last_col()
    ) %>%
  dplyr::filter(
    dplyr::if_any(
      8:dplyr::last_col(),
      ~ !is.na(.x)
    )) %>%
  pivot_longer(
    cols = -c(sample_id, region, region_lab, wb_id, wb,sample_date, yday, seas),
    names_to = c("taxon", ".value"),
    names_pattern = "(.*)_(abund_m3|ugC_m3)"
  ) %>%
  rename(
    abundance_m3 = abund_m3,
    ug_c_m3 = ugC_m3
  ) %>% 
  # retain rows the have either abundance or carbon values
  dplyr::filter(!is.na(abundance_m3) | !is.na(ug_c_m3)
                ) %>% 
  ## arrange by sample_id
  dplyr::arrange(.,sample_id) -> df_trim_l

tictoc::toc(log = TRUE)

tictoc::tic("Generate summed taxa by sample dataframe")
# Generate summed taxa by sample dataframe ----

df_trim_l %>% 
  dplyr::group_by(dplyr::across(-c(abundance_m3, ug_c_m3))) %>% 
  dplyr::summarise(abundance_m3 = sum(abundance_m3, na.rm = TRUE),
                   ug_c_m3 = sum(ug_c_m3, na.rm = TRUE),
                   .groups = "drop") %>% dplyr::ungroup() %>% 
  ## arrange by sample_id
  dplyr::arrange(.,sample_id) -> df_trim_l_summary

# write data

write.csv(df_trim_l_summary,
          file = paste0(output.location,
                        "trimmed_phyto_taxa_of_interest_long.csv"),
          row.names = FALSE
          )

tictoc::toc(log = TRUE)

tictoc::tic("Produce correlation heatmap")
# correlation heatmap ----
taxa_mat <-
  df_trim_l_summary %>%
  select(sample_id, taxon, abundance_m3) %>%
  group_by(sample_id, taxon) %>%
  summarise(
    abundance_m3 = sum(abundance_m3),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = taxon,
    values_from = abundance_m3,
    values_fill = 0
  )

# Load correlation package
library(correlation);library(reshape2)

cor_mat <-
  taxa_mat %>%
  select(-sample_id) %>%
  cor(
    method = "spearman",
    use = "pairwise.complete.obs"
  )

# Keep lower triangle only
cor_df <- as_tibble(melt(cor_mat))

cor_df <- cor_df %>%
  dplyr::mutate(Var1 = stringr::str_to_camel(as.character(Var1),first_upper = TRUE),
                Var2 = stringr::str_to_camel(as.character(Var2),first_upper = TRUE),
                Var1 = as.factor(Var1),
                Var2 = as.factor(Var2),
  ) %>% 
  dplyr::filter(
    as.numeric(Var1) > as.numeric(Var2)
  )

ggplot(
  cor_df,
  aes(
    x = stringr::str_to_camel(as.character(Var1),first_upper = TRUE),
    y = stringr::str_to_camel(as.character(Var2),first_upper = TRUE),
    fill = value,
  )
) +
  geom_hline(yintercept = seq(from = 0.5,
                              to = length(unique(cor_df$Var1)),
                              by = 1),
             colour = "lightgrey",
             lty = 2)+
  geom_vline(xintercept = seq(from = 0.5,
                              to = length(unique(cor_df$Var1)),
                              by = 1),
             colour = "lightgrey",
             lty = 2)+
  geom_tile(colour = "grey") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Spearman\nrho"
  ) +
  coord_equal() +
  # theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )
  ) +
  labs(
    title = "Correlation of abundances of selected phytoplankton taxa",
    caption = paste0("Phytoplankton data gathered between ",min(year(df_trim_l_summary$sample_date))," & ", min(year(df_trim_l_summary$sample_date))),
    x = NULL,
    y = NULL
  )+
  theme(
    plot.title = element_text(face = 2),
    axis.text = element_text(face = 2),
    )

rm(cor_df, cor_mat, taxa_mat)
tictoc::toc(log=TRUE)  

# Plot
df_trim_l_summary %>%
  dplyr::mutate(taxon = stringr::str_to_camel(as.character(taxon),
                                              first_upper = TRUE)) %>% 
  # remove older data which is potentially more clunky(?)
  # dplyr::filter(sample_date > "2009-12-31") %>% 
  ggplot(.,
         aes(
           x = yday,
           y = log10(abundance_m3+1),
           group = taxon,
         )
  )+
  geom_point(
    # aes(colour = taxon),
    alpha = 0.05)+
  geom_smooth(
    # aes(colour = taxon),
    method = "gam",
    se=FALSE)+
  facet_wrap(.~taxon, scales = "free_y")+
  labs(
    title = "Abundance of selected phytoplankton taxa throughout the year",
    caption = paste0("Phytoplankton data gathered between ",
                     min(year(df_trim_l_summary$sample_date)),
                     " & ", min(year(df_trim_l_summary$sample_date)),".","\n",
                     "Note: scales on y axis are independent."),
    x = "Day of year",
    y = "log10(abundance per m3 +1)"
    )+
  theme(
    plot.title = element_text(face = 2),
    plot.caption = element_text(face = 2),
    strip.text = element_text(face = 2),
    axis.title = element_text(face = 2),
    )
