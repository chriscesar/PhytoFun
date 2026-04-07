# Solway_dig.R ####

# load packages ####
ld_pkgs <- c("tidyverse","tictoc","sf","units")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log
tictoc::tic("TOTAL TIME")
source("R/helperFunctions.R")

# tictoc::tic("Load_data_exported_from_BOXI_Query")
# df0 <- readxl::read_xlsx("data/Phyto_2000_2025.xlsx",sheet = "PhyoData",
#                                guess_max = 2000) %>% 
#   janitor::clean_names()
# saveRDS(janitor::clean_names(df_phyto0), file = "outputs/ospar_Phyto_2000_2025.raw.Rdat")
# tictic::toc(log = TRUE)

df0 <- readRDS("outputs/ospar_Phyto_2000_2025.raw.Rdat")

# dfnw <- df0 %>% 
#   filter(river_basi == "North West"|river_basi =="Solway Tweed")
# 
# unique(dfnw$taxa_name)

# Convert taxon names refs and update to 'correct' Aphia ID ####
# taxa <- unique(dfnw$taxa_name)
taxa <- unique(df0$taxa_name)

# extract Aphia IDs  ####
tic("Query AphiaID for each name")
aphia_ids <- map_dfr(taxa, function(taxon) {
  result <- tryCatch(
    {
      # Safe querying
      res <- worrms::wm_name2id(name = taxon)
      
      # Ensure a tibble/data frame with expected fields
      if (length(res) == 0) {
        tibble(name = taxon, aphia_id = NA_integer_)
      } else {
        tibble(name = taxon, aphia_id = res)
      }
    },
    error = function(e) {
      tibble(name = taxon, aphia_id = NA_integer_)
    }
  )
})
toc(log=TRUE)

### Import missing Aphia values ####
blankAPHIA <- read.csv("outputs/MissingAphiaLookupImport.csv") #created manually
#

# tictoc::tic("retrieve record details: non-na")
# # Filter out rows with NA aphia_id
# aphia_ids_non_na <- aphia_ids %>% 
#   dplyr::filter(!is.na(aphia_id)) %>% 
#   dplyr::mutate(
#     RLIST = "ERID"
#   )
# 
# # For each aphia_id, retrieve record details
# record_info_non_na <- aphia_ids_non_na |> 
#   dplyr::mutate(
#     wm_data = purrr::map(aphia_id, ~ tryCatch(worrms::wm_record(.x),
#                                               error = function(e) NULL))
#   ) |>
#   tidyr::unnest_wider(wm_data)
# tictoc::toc(log=TRUE)
# 
# 
# tictoc::tic("retrieve record details: na")
# # For each aphia_id, retrieve record details
# record_info_missing <- blankAPHIA |> 
#   mutate(
#     wm_data = purrr::map(APHIA_ID, ~ tryCatch(worrms::wm_record(.x),
#                                               error = function(e) NULL))
#   ) |>
#   tidyr::unnest_wider(wm_data)
# record_info_missing |> rename("taxa_name" = "Taxon") -> record_info_missing
# tictoc::toc(log=TRUE)
# 
# # Combine WORMS records ####
# tic("Combine WORMS records")
# 
# # tidy names
# record_info_missing0 <- janitor::clean_names(record_info_missing) %>% 
#   dplyr::select(.,-notes, -flag)
# record_info_non_na0 <- janitor::clean_names(record_info_non_na) %>% dplyr::rename(taxa_name=name)
# 
# records <- rbind(record_info_non_na0,record_info_missing0) %>%
#   dplyr::rename(reported_taxon = taxa_name)

tictoc::tic("Get taxon IDs")
## =========================~
## Record details: non-NA IDs ####
## =========================~
tictoc::tic("retrieve record details: non-na")

aphia_ids_non_na <- aphia_ids %>%
  dplyr::filter(!is.na(aphia_id)) %>%
  dplyr::mutate(RLIST = "ERID")

# Core record info
record_info_non_na <- aphia_ids_non_na |>
  dplyr::mutate(
    wm_data = purrr::map(aphia_id, ~ tryCatch(worrms::wm_record(.x),
                                              error = function(e) NULL))
  ) |>
  tidyr::unnest_wider(wm_data)

tictoc::toc(log = TRUE)

## =========================~
## Record details: NA branch (blankAPHIA) ####
## =========================~
tictoc::tic("retrieve record details: na")

record_info_missing <- blankAPHIA |>
  dplyr::mutate(
    wm_data = purrr::map(APHIA_ID, ~ tryCatch(worrms::wm_record(.x),
                                              error = function(e) NULL))
  ) |>
  tidyr::unnest_wider(wm_data)

record_info_missing |> dplyr::rename("taxa_name" = "Taxon") -> record_info_missing

tictoc::toc(log = TRUE)

## =========================~
## Tidy names (as in your code) ####
## =========================~
tic("Combine WORMS records")

record_info_missing0 <- janitor::clean_names(record_info_missing) %>%
  dplyr::select(-notes, -flag)

record_info_non_na0 <- janitor::clean_names(record_info_non_na) %>%
  dplyr::rename(taxa_name = name)

records <- dplyr::bind_rows(record_info_non_na0, record_info_missing0) %>%
  dplyr::rename(reported_taxon = taxa_name)

## =========================~
## NEW: Full classification per Aphia ID ####
## =========================~
# Helper: build classification wide table for a given vector of Aphia IDs
build_classification_wide <- function(aphia_vec, id_col_name = "aphia_id", names_prefix = "") {
  # Fetch classification list for each Aphia ID
  class_list <- tibble::tibble(!!id_col_name := aphia_vec) %>%
    dplyr::mutate(
      wm_class = purrr::map(
        !!rlang::sym(id_col_name),
        ~ tryCatch(worrms::wm_classification(.x), error = function(e) NULL)
      )
    ) %>%
    tidyr::unnest(wm_class) %>%
    # Some calls may return NULL; drop those
    dplyr::filter(!is.na(AphiaID) | !is.na(scientificname) | !is.na(rank)) %>%
    dplyr::mutate(rank = tolower(rank)) %>%
    # avoid rare duplicated ranks; keep first occurrence
    dplyr::group_by(!!rlang::sym(id_col_name), rank) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    # Pivot rank -> columns, values = scientificname
    tidyr::pivot_wider(
      id_cols = !!rlang::sym(id_col_name),
      names_from = rank,
      values_from = scientificname,
      names_prefix = names_prefix
    )
  
  class_list
}

## =========================~
## Classification for NON-NA records (using aphia_id) ####
## =========================~
class_non_na_wide <- build_classification_wide(aphia_ids_non_na$aphia_id, id_col_name = "aphia_id")

# Join into your existing non-NA table
record_info_non_na_enh <- record_info_non_na0 %>%
  dplyr::left_join(class_non_na_wide, by = "aphia_id")

## =========================~
## Classification for NA branch ####
## =========================~
# For blankAPHIA we’ll assume classification should use APHIA_ID if present post-unnest.
# After clean_names, the column is likely 'aphia_id' already, but we guard against naming.
missing_id_col <- if ("aphia_id" %in% names(record_info_missing0)) "aphia_id" else "APHIA_ID"

class_missing_wide <- build_classification_wide(
  aphia_vec = record_info_missing0[[missing_id_col]],
  id_col_name = missing_id_col
)

record_info_missing_enh <- record_info_missing0 %>%
  dplyr::left_join(class_missing_wide, by = setNames("aphia_id", missing_id_col))

## =========================~
## Final combined with classification ####
## =========================~
records_with_classification <- dplyr::bind_rows(
  record_info_non_na_enh,
  record_info_missing_enh
) %>%
  dplyr::rename(reported_taxon = taxa_name)

tictoc::toc(log = TRUE)

# =========================~
## OPTIONAL: Prefer valid AphiaID classification instead of original ####
# =========================~
# If you want ranks based on valid_AphiaID (more stable), build a second wide table and join.
if ("valid_aphia_id" %in% names(records_with_classification)) {
  class_valid_wide <- build_classification_wide(
    aphia_vec = records_with_classification$valid_aphia_id,
    id_col_name = "valid_aphia_id",
    names_prefix = "valid_"
  )
  
  records_with_classification <- records_with_classification %>%
    dplyr::left_join(class_valid_wide, by = "valid_aphia_id")
}
tictoc::toc(log=TRUE)

tictoc::tic("Pull out required taxon variables & append to data")
# dfout <- dfnw
dfout <- df0

dfout %>% 
  tibble::as_tibble() %>%  
  # dplyr::rename(speci = value) %>% 
  dplyr::left_join(
    .,records_with_classification,
    # records_with_classification, %>% 
    #   dplyr::select(
    #     valid_kingdom,valid_infrakingdom,valid_phylum,valid_subkingdom,
    #     valid_class,"valid_phylum (division)", "valid_subphylum (subdivision)",
    #     "valid_family","valid_genus","valid_infraphylum","valid_order",
    #     "valid_subclass","valid_subfamily","valid_subphylum","valid_species",
    #     "valid_gigaclass","valid_parvphylum","valid_superclass",
    #     "valid_infraclass","valid_subterclass","valid_superfamily",
    #     "valid_superorder","valid_forma","valid_variety","valid_suborder",
    #     "valid_subgenus") %>% 
    #     
    #     
    #  )
      # dplyr::select(reported_taxon, rlist,aphia_id,aphia_id_2,valid_aphia_id,
      #               scientificname, 
      #               kingdom,phylum,class,order,
      #               family,genus,
      #               
      #               valid_name) %>% 
      # dplyr::rename(speci = reported_taxon,
      #               rlist = rlist),
    by = c("taxa_name"= "reported_taxon")
  ) ->dfout

# dfout %>% 
#   dplyr::filter(river_basi=="North West"|river_basi=="Solway Tweed") %>% 
#   write.csv(., file = "outputs/Solway_region_phyto.csv",row.names = FALSE)

dino <- dfout %>%
  dplyr::filter(!is.na(valid_infraphylum)) %>% 
  dplyr::filter(infraphylum=="Dinoflagellata")

dino %>% dplyr::filter(river_basi=="North West"|river_basi=="Solway Tweed") -> nwdino


# Create a sequence of January 1st for each year spanned by your data
year_lines <- tibble(
  year_start = seq(
    floor_date(min(nwdino$sample_date, na.rm = TRUE), unit = "year"),
    floor_date(max(nwdino$sample_date, na.rm = TRUE), unit = "year"),
    by = "1 year"
  )
)

year_lines5<- tibble(
  year_start = seq(
    floor_date(min(nwdino$sample_date, na.rm = TRUE), unit = "year"),
    floor_date(max(nwdino$sample_date, na.rm = TRUE), unit = "year"),
    by = "5 years"
  )
)

# DINOFLAG PLOT ####
nwdino %>% 
  dplyr::select(.,
                sample_date, cells_per_litre_millilitre, wb_name,sample_id) %>% 
  group_by(across(!cells_per_litre_millilitre)) %>% 
  summarise(cells_per_litre_millilitre = sum(cells_per_litre_millilitre,na.rm = TRUE),
            .groups = "drop") %>% 
  ungroup() %>% 
  ggplot(.,
         aes(
           x = sample_date,
           y = log10(cells_per_litre_millilitre+1),
           ))+
  geom_vline(data = year_lines,
             aes(xintercept = as.numeric(year_start)),
             color = "grey70", linewidth = 0.3,lty=3)+
  geom_vline(data = year_lines5,
             aes(xintercept = as.numeric(year_start)),
             color = "grey50", linewidth = 0.3,lty=2)+
  geom_point()+
  facet_wrap(.~wb_name)+
  ggthemes::theme_few()+
  geom_smooth(method="gam")+
  labs(
    title = "Dinoflagellate abundances recorded in selected English water bodies",
    subtitle = "Data faceted by Environment Agency water body",
    caption = "EA data filtered to retain taxa belonging to Infraphylum Dinoflagellata",
    y="Log10 cells per litre (n+1)"
    )+
  theme(
    strip.text = element_text(face=2,size = 14),
    axis.text = element_text(face=2,size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2,size = 14),
    plot.title = element_text(face=2,size = 16),
    plot.subtitle = element_text(face=2,size = 14),
    plot.caption = element_text(face=2,size = 14),
  )

ggsave(plot = get_last_plot(),
       filename = paste0("outputs/figs/SolwayNW_Dinoflag.png"),
       width = 18, height = 8, units = "in"
       )

# DIATOM PLOT ####
diat <- dfout %>%
  # dplyr::filter(!is.na(valid_infraphylum)) %>% 
  dplyr::filter(subphylum=="Bacillariophytina")

diat %>%
  dplyr::filter(river_basi=="North West"|river_basi=="Solway Tweed") -> nwdiat

nwdiat %>% 
dplyr::select(.,
              sample_date, cells_per_litre_millilitre, wb_name,sample_id) %>% 
  group_by(across(!cells_per_litre_millilitre)) %>% 
  summarise(cells_per_litre_millilitre = sum(cells_per_litre_millilitre,na.rm = TRUE),
            .groups = "drop") %>% 
  ungroup() %>% 
  ggplot(.,
         aes(
           x = sample_date,
           y = log10(cells_per_litre_millilitre+1),
         ))+
  geom_vline(data = year_lines,
             aes(xintercept = as.numeric(year_start)),
             color = "grey70", linewidth = 0.3,lty=3)+
  geom_vline(data = year_lines5,
             aes(xintercept = as.numeric(year_start)),
             color = "grey50", linewidth = 0.3,lty=2)+
  geom_point()+
  facet_wrap(.~wb_name)+
  ggthemes::theme_few()+
  geom_smooth(method="gam")+
  labs(
    title = "Diatom abundances recorded in selected English water bodies",
    subtitle = "Data faceted by Environment Agency water body",
    caption = "EA data filtered to retain taxa belonging to Subphylum Bacillariophytina",
    y="Log10 cells per litre (n+1)"
  )+
  theme(
    strip.text = element_text(face=2,size = 14),
    axis.text = element_text(face=2,size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2,size = 14),
    plot.title = element_text(face=2,size = 16),
    plot.subtitle = element_text(face=2,size = 14),
    plot.caption = element_text(face=2,size = 14),
  )

ggsave(plot = get_last_plot(),
       filename = paste0("outputs/figs/SolwayNW_Diatoms.png"),
       width = 18, height = 8, units = "in"
       )
