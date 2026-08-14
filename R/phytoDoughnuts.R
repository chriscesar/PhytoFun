# phytoDoughnuts.R ----

# Aims: Import water quality & phytoplankton data & explore cycles over time
# Author: Christopher Cesar
# Created: 04/08/2026
# Last edited:
# Notes:

# load packages & data ----
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

## WIMS ----
### define date
df$wims <- df$wims %>%
  dplyr::mutate(
    date = stringr::str_sub(code, start = -8),
    date = lubridate::ymd(date)
  ) %>%
  relocate(date, .after = code) %>% 
  dplyr::mutate(wb_id = stringr::str_sub(code, end = -10,)) %>% 
  relocate(wb_id, .after = code) %>% 
  dplyr::left_join(., wb_inf,by = "wb_id") %>% 
  dplyr::relocate(wb, .after = wb_id)


summary(nchar(df0$wims$code))

names(df$meta)

# df$wims %>%
#   dplyr::filter(wb == "WYRE") %>%
#   dplyr::select_if(~ !all(is.na(.))) %>% #names()
#   {ggplot2::ggplot(
#     .,
#     aes(
#       x = chlorophyll_ug_l,
#       # y = temp_water_cel,
#       y = n_dis_inorg_mg_l,
#     )
#   ) +
#       ggplot2::geom_point() +
#       ggplot2::labs(title = unique(.$wb))
#   }

# # PCA over time ----
# 
# ## Pull data from one region/WB ----
# table(df$meta$region,df$meta$wb)
# rg <- "South West"
# wb_use <- "POOLE HARBOUR"
# n <- 3 # minimum number of samples for taxon to be retained
# 
# ## Data prep ----
# # Which rows correspond to that region?
# df$meta %>% 
#   filter(region == rg) %>% droplevels(.) -> rgn_meta
# 
# df$abund_lf_per_l %>% 
#   filter(sample_id_biosys %in% rgn_meta$sample_id_biosys) %>% 
#   droplevels(.) -> rgn_lf
# 
# rgn_meta %>% 
#   dplyr::select(sample_id_biosys, wb, sample_date) %>% 
#   left_join(., rgn_lf, by = "sample_id_biosys") %>% 
#   # drop unneccessary cols
#   dplyr::select(-c(sample_id_biosys,code))-> rgn_df
# rm(rgn_lf, rgn_meta)
# 
# # fill NAs with zero & generate year variable
# rgn_df_clean <- rgn_df %>%
#   # create year variable and delete date
#   dplyr::mutate(year = lubridate::year(sample_date)) %>% 
#   relocate(year) %>% 
#   dplyr::select(-sample_date) %>% 
#   # Remove columns that are entirely NA
#   select(where(~ !all(is.na(.)))) %>%
#   # Replace remaining NAs with 0 in numeric columns
#   mutate(across(where(is.numeric), ~ tidyr::replace_na(., 0)))
# 
# rgn_mean <- rgn_df_clean %>%
#   group_by(wb, year) %>%
#   summarise(
#     across(
#       where(is.numeric),
#       ~ mean(.x, na.rm = TRUE)
#     ),
#     .groups = "drop"
#   ) %>%
#   # remove variable "NA"
#   dplyr::select(-`NA`) %>%
#   # choose WB of interest
#   dplyr::filter(wb == wb_use) %>%
#   # Remove obsolete variables
#   dplyr::select(-c(wb))
# 
# # Extract year separately
# years <- rgn_mean$year
# 
# # PCA input data (exclude year)
# 
# pca_data <- vegan::decostand(rgn_mean %>%
#   select(-year),
#   method = "hellinger"
#     ) %>%
#   # Keep taxa with >= n non-zero values
#   select(where(~ sum(.x > 0, na.rm = TRUE) >= n)) %>%
#   # Log transform
#   mutate(across(everything(), ~ log10(.x + 1)))
# 
# ## Run PCA ----
# pca <- prcomp(
#   pca_data,
#   center = TRUE,
#   scale. = TRUE
#   )
# 
# # PCA scores (site/year coordinates)
# scores <- as.data.frame(pca$x) %>%
#   mutate(year = years) %>%
#   arrange(year)
# 
# # PCA loadings
# loadings <- as.data.frame(pca$rotation) %>%
#   rownames_to_column("taxon")
# 
# # Variance explained
# var_explained <- summary(pca)$importance[2, ] * 100
# 
# # Scale loadings to fit score space
# mult <- 5
# 
# ## Plot ----
# ggplot(scores, aes(PC1, PC2)) +
#   geom_path(
#   aes(colour = year),
#   linewidth = 1.2,
#   arrow = arrow(length = unit(0.15, "cm"))
#   ) +
#   geom_point(
#     aes(colour = year),
#     size = 4
#   ) +
#   scale_colour_viridis_c() +
#   # loadings
#   geom_segment(
#     data = loadings,
#     aes(
#       x = 0,
#       y = 0,
#       xend = PC1 * mult,
#       yend = PC2 * mult,
#       ),
#     linewidth = .75,
#     linetype = 2,
#     colour = "red",
#     arrow = arrow(length = unit(0.2, "cm"))
#   ) +
#   geom_text_repel(
#     data = loadings,
#     aes(
#       x = PC1 * mult,
#       y = PC2 * mult,
#       label = taxon
#     ),
#     fontface = "bold",
#     colour = "red",
#     size = 5
#   )+
#   labs(
#     title = paste0(rg,": ", wb_use),
#     x = paste0("PC1 (", round(var_explained[1], 1), "% variance explained)"),
#     y = paste0("PC2 (", round(var_explained[2], 1), "% variance explained)")
#     )+
#   theme_bw()+
#   theme(
#     plot.title = element_text(face = 2),
#     axis.title = element_text(face = 2),
#     )

# loop plot ----
n <- 3  # minimum number of samples for taxon to be retained

## Create regional dataset (your existing workflow) ----

rgn_meta <- df$meta %>%
  droplevels()

rgn_lf <- df$abund_lf_per_l

rgn_df <- rgn_meta %>%
  dplyr::select(sample_id_biosys, region, wb, sample_date) %>%
  left_join(rgn_lf, by = "sample_id_biosys") %>%
  dplyr::select(-c(sample_id_biosys, code))

rgn_df_clean <- rgn_df %>%
  dplyr::mutate(year = lubridate::year(sample_date)) %>%
  relocate(year) %>%
  dplyr::select(-sample_date) %>%
  dplyr::select(where(~ !all(is.na(.)))) %>%
  dplyr::mutate(
    dplyr::across(
      where(is.numeric),
      ~ tidyr::replace_na(.x, 0)
    )
  )

## Lookup table of WB-region combinations ----

wb_lookup <- rgn_df_clean %>%
  dplyr::distinct(region, wb)

## Output folder ----

outdir <- "outputs/PCA_plots"
dir.create(outdir, showWarnings = FALSE)

## Function to create PCA plot ----
make_pca_plot <- function(wb_use, region_use, dat, n = 3) {
  
  rgn_mean <- dat %>%
    dplyr::filter(wb == wb_use) %>%
    dplyr::group_by(wb, year) %>%
    dplyr::summarise(
      dplyr::across(
        where(is.numeric),
        ~ mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    ) %>%
    # remove variable "NA"
    dplyr::select(-`NA`) %>%
    dplyr::select(-wb)
  
  # Need enough years for PCA
  if(nrow(rgn_mean) < 3) {
    message("Skipping ", wb_use, ": insufficient years")
    return(NULL)
  }
  
  #--------------------------------------------------
  # Remove years with no taxa present
  #--------------------------------------------------
  
  keep_rows <- rowSums(
    rgn_mean %>% dplyr::select(-year)
  ) > 0
  
  if(sum(keep_rows) < 3) {
    message(
      "Skipping ", wb_use,
      ": fewer than 3 years with data"
    )
    return(NULL)
  }
  
  years <- rgn_mean$year[keep_rows]
  
  #--------------------------------------------------
  # PCA input data
  #--------------------------------------------------
  
  pca_data <- rgn_mean %>%
    dplyr::filter(keep_rows) %>%
    dplyr::select(-year) %>%
    
    # Hellinger transform
    vegan::decostand(method = "hellinger") %>%
    
    # Retain taxa present in >= n years
    dplyr::select(
      where(~ sum(.x > 0, na.rm = TRUE) >= n)
    ) %>%
    
    # Log transform
    dplyr::mutate(
      dplyr::across(
        everything(),
        ~ log10(.x + 1)
      )
    ) %>%
    
    # Remove taxa with zero variance
    dplyr::select(
      where(~ stats::sd(.x, na.rm = TRUE) > 0)
    )
  
  #--------------------------------------------------
  # Sanity checks
  #--------------------------------------------------
  
  if(nrow(pca_data) < 3) {
    message(
      "Skipping ", wb_use,
      ": insufficient observations"
    )
    return(NULL)
  }
  
  if(ncol(pca_data) < 2) {
    message(
      "Skipping ", wb_use,
      ": insufficient taxa"
    )
    return(NULL)
  }
  
  if(any(!is.finite(as.matrix(pca_data)))) {
    message(
      "Skipping ", wb_use,
      ": non-finite values after transformation"
    )
    return(NULL)
  }
  
  #--------------------------------------------------
  # PCA
  #--------------------------------------------------
  
  pca <- tryCatch(
    stats::prcomp(
      pca_data,
      center = TRUE,
      scale. = TRUE
    ),
    error = function(e) {
      message(
        "Skipping ", wb_use,
        ": ", e$message
      )
      return(NULL)
    }
  )
  
  if(is.null(pca)) {
    return(NULL)
  }
  
  scores <- as.data.frame(pca$x) %>%
    dplyr::mutate(year = years) %>%
    dplyr::arrange(year)
  
  loadings <- as.data.frame(pca$rotation) %>%
    tibble::rownames_to_column("taxon")
  
  var_explained <- summary(pca)$importance[2, ] * 100
  
  mult <- 5
  
  p <- ggplot2::ggplot(
    scores,
    ggplot2::aes(PC1, PC2)
  ) +
    ggplot2::geom_path(
      ggplot2::aes(colour = year),
      linewidth = 1.2,
      arrow = grid::arrow(length = unit(0.15, "cm"))
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = year),
      size = 4
    ) +
    ggplot2::scale_colour_viridis_c(direction = -1) +
    ggplot2::geom_segment(
      data = loadings,
      ggplot2::aes(
        x = 0,
        y = 0,
        xend = PC1 * mult,
        yend = PC2 * mult
      ),
      colour = "red",
      linetype = 2,
      linewidth = 0.75,
      arrow = grid::arrow(length = unit(0.2, "cm"))
    ) +
    ggrepel::geom_text_repel(
      data = loadings,
      ggplot2::aes(
        x = PC1 * mult,
        y = PC2 * mult,
        label = taxon
      ),
      colour = "red",
      fontface = "bold",
      size = 4
    ) +
    ggplot2::labs(
      title = paste0(region_use, ": ", wb_use),
      x = paste0(
        "PC1 (",
        round(var_explained[1], 1),
        "% variance explained)"
      ),
      y = paste0(
        "PC2 (",
        round(var_explained[2], 1),
        "% variance explained)"
      ),
      caption = paste0("Principal components analysis of Hellinger-transformed annual mean abundances of phytoplankton groups.\n",
                       "Taxa appearing in ",n," or fewer samples excluded.\n",
                       "Dashed lines indicate PCA loadings. Coloured lines and points show trends over time.")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = 2),
      axis.title = ggplot2::element_text(face = 2)
    )
  
  return(p)
  }

safe_make_pca_plot <- purrr::possibly(
  make_pca_plot,
  otherwise = NULL
)

purrr::pwalk(
  wb_lookup,
  function(region, wb) {
    
    p <- safe_make_pca_plot(
      wb_use = wb,
      region_use = region,
      dat = rgn_df_clean,
      n = n
    )
    
    if(!is.null(p)) {
      
      fname <- paste0(
        outdir, "/",
        gsub("[^A-Za-z0-9_]", "_", region),
        "_",
        gsub("[^A-Za-z0-9_]", "_", wb),
        "_PCA.png"
      )
      
      ggplot2::ggsave(
        filename = fname,
        plot = p,
        width = 10,
        height = 8,
        dpi = 300
      )
      
      message("Saved: ", fname)
    }
  }
)


## Loop through all waterbodies ----
tictoc::tic("Generate individual plots")
purrr::pwalk(
  wb_lookup,
  function(region, wb) {
    
    p <- safe_make_pca_plot(
      wb_use = wb,
      region_use = region,
      dat = rgn_df_clean,
      n = n
    )
    
    if(!is.null(p)) {
      
      fname <- paste0(
        outdir, "/",
        gsub("[^A-Za-z0-9_]", "_", region),
        "_",
        gsub("[^A-Za-z0-9_]", "_", wb),
        "_PCA.png"
      )
      
      ggplot2::ggsave(
        filename = fname,
        plot = p,
        width = 10,
        height = 8,
        dpi = 300
      )
      
      message("Saved: ", fname)
    }
    }
  )
tictoc::toc(log=TRUE)

## save as single pdf ----
wb_lookup_sorted <- wb_lookup %>%
  dplyr::arrange(region, wb)

tictoc::tic("Generate plots in pdf")
pdf("outputs/PCA_plots/All_WB_PCA.pdf", width = 10, height = 8)

for(i in seq_len(nrow(wb_lookup))) {
  
  p <- safe_make_pca_plot(
    wb_lookup_sorted$wb[i],
    wb_lookup_sorted$region[i],
    rgn_df_clean,
    n = n
    )
  if(!is.null(p)) print(p)
  }
dev.off()
tictoc::toc(log=TRUE)

unlist(tictoc::tic.log())
