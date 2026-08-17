# phytoDoughnuts_spt.R ----
# Create PCA plot of individual sample points ---

# load pacakges ----
ld_pkgs <- c("tidyverse","tictoc","ggfortify", "plotly","ggrepel")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

## load metadate ----
source("R/000setup.R")

## load phyto & wims data ----
df0 <- readRDS("outputs/phyto_wims_list_use.Rdat")

# Data prep ----
df <- df0 # create working version separate from backup
df$meta %>% dplyr::select(wb_id,wb) %>% dplyr::distinct() -> wb_inf

# Prep data ----
n <- 3  # minimum number of samples for taxon to be retained

## Create regional dataset (your existing workflow) ----
rgn_meta <- df$meta %>%
  droplevels()

rgn_lf <- df$abund_lf_per_l

rgn_df <- rgn_meta %>%
  dplyr::select(sample_id_biosys, region, wb, biosys_short,sample_date) %>%
  left_join(rgn_lf, by = "sample_id_biosys") %>%
  dplyr::select(-c(
    # sample_id_biosys,
    code
    )) %>% 
  dplyr::filter(!is.na(biosys_short))

# which sample point has most data?

rgn_df %>%
  filter(
    biosys_short == names(which.max(table(biosys_short)))
  ) -> pca_data0

## calculate annual means ----


#   mutate(across(where(is.numeric), ~ tidyr::replace_na(., 0)))

pca_data0 %>% 
  dplyr::mutate(year = lubridate::year(sample_date)) %>% 
  dplyr::relocate(year) %>% 
  # Replace NAs with 0 in numeric columns
  mutate(across(where(is.numeric), ~ tidyr::replace_na(., 0))) %>% 
  # remove sample ID & meta data
  dplyr::select(-c(sample_id_biosys,region,wb,sample_date)) %>% 
  dplyr::group_by(biosys_short, year) %>%
  dplyr::summarise(
    dplyr::across(
      where(is.numeric),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>% #-> pca_data_mn
  dplyr::select(-dplyr::any_of("NA")) -> pca_data_mn
#   
# vegan::decostand(rgn_mean %>%
#                                #   select(-year),
#                                #   method = "hellinger"
                               #     ) %>%
                               #   # Keep taxa with >= n non-zero values
                               #   select(where(~ sum(.x > 0, na.rm = TRUE) >= n)) %>%
                               #   # Log transform
                               #   mutate(across(everything(), ~ log10(.x + 1)))

n <- 3 # min number samples to retain taxon

# drop meta & under-represented taxa
pca_data_mn %>% dplyr::select(-c(biosys_short,year)) %>% 
    # Keep taxa with >= n non-zero values
    select(where(~ sum(.x > 0, na.rm = TRUE) >= n)) -> pca_data
#--------------------------------------------------
# PCA
#--------------------------------------------------

pca <- 
  stats::prcomp(
    pca_data,
    center = TRUE,
    scale. = TRUE
  )

# scores <- as.data.frame(pca$x) %>%
#   dplyr::mutate(year = years) %>%
#   dplyr::arrange(year) %>%
#   dplyr::mutate(
#     PC1_next = dplyr::lead(PC1),
#     PC2_next = dplyr::lead(PC2)
#   )
# 
# loadings <- as.data.frame(pca$rotation) %>%
#   tibble::rownames_to_column("taxon")
# 
# var_explained <- summary(pca)$importance[2, ] * 100
# 
# mult <- 5
# 
# p <- ggplot2::ggplot(
#   scores,
#   ggplot2::aes(PC1, PC2)
# ) +
#   ggplot2::geom_segment(
#     data = dplyr::filter(scores, !is.na(PC1_next)),
#     ggplot2::aes(
#       x = PC1,
#       y = PC2,
#       xend = PC1_next,
#       yend = PC2_next,
#       colour = year
#     ),
#     linewidth = 1.2,
#     arrow = grid::arrow(
#       length = unit(0.2, "cm"),
#       type = "open"
#     )
#   ) +
#   ggplot2::geom_point(
#     ggplot2::aes(colour = year),
#     size = 4
#   ) +
#   ggplot2::scale_colour_viridis_c(direction = -1) +
#   ggplot2::geom_segment(
#     data = loadings,
#     ggplot2::aes(
#       x = 0,
#       y = 0,
#       xend = PC1 * mult,
#       yend = PC2 * mult
#     ),
#     colour = "red",
#     linetype = 2,
#     linewidth = 0.75,
#     arrow = grid::arrow(length = unit(0.2, "cm"))
#   ) +
#   ggrepel::geom_text_repel(
#     data = loadings,
#     ggplot2::aes(
#       x = PC1 * mult,
#       y = PC2 * mult,
#       label = taxon
#     ),
#     colour = "red",
#     fontface = "bold",
#     size = 4
#   ) +
#   ggplot2::labs(
#     title = region_use,
#     x = paste0(
#       "PC1 (",
#       round(var_explained[1], 1),
#       "% variance explained)"
#     ),
#     y = paste0(
#       "PC2 (",
#       round(var_explained[2], 1),
#       "% variance explained)"
#     ),
#     caption = paste0(
#       "Principal components analysis of Hellinger-transformed annual mean abundances of phytoplankton groups.\n",
#       "Taxa appearing in fewer than ", n, " years excluded.\n",
#       "Dashed lines indicate PCA loadings. Coloured lines and points show trends over time.\n",
#       "Excluded taxa (if any): ",
#       paste(exclude_taxa, collapse = ", "), "."
#     )
#   ) +
#   ggplot2::theme_bw() +
#   ggplot2::theme(
#     plot.title = ggplot2::element_text(face = 2),
#     axis.title = ggplot2::element_text(face = 2)
#   )
# 
# p