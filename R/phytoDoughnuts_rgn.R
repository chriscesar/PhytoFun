# phytoDoughnuts_rgn.R ----
# Create PCA plots by region ---

# First run data import from phytoDoughnuts.R, then run code below to plot by region

## Lookup table of regions ----

region_lookup <- rgn_df_clean %>%
  dplyr::distinct(region)

## Output folder ----

outdir <- "outputs/PCA_plots"
dir.create(outdir, showWarnings = FALSE)

## Function to create PCA plot ----

make_pca_plot <- function(region_use, dat, n = 3) {
  
  rgn_mean <- dat %>%
    dplyr::filter(region == region_use) %>%
    dplyr::group_by(region, year) %>%
    dplyr::summarise(
      dplyr::across(
        where(is.numeric),
        ~ mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    ) %>%
    dplyr::select(-dplyr::any_of("NA")) %>%
    dplyr::select(-region)
  
  # Need enough years for PCA
  if (nrow(rgn_mean) < 3) {
    message("Skipping ", region_use, ": insufficient years")
    return(NULL)
  }
  
  #--------------------------------------------------
  # Remove years with no taxa present
  #--------------------------------------------------
  
  keep_rows <- rowSums(
    rgn_mean %>% dplyr::select(-year)
  ) > 0
  
  if (sum(keep_rows) < 3) {
    message(
      "Skipping ", region_use,
      ": fewer than 3 years with data"
    )
    return(NULL)
  }
  
  years <- rgn_mean$year[keep_rows]
  
  #--------------------------------------------------
  # PCA input data
  #--------------------------------------------------
  
  pca_data <- rgn_mean %>%
    dplyr::select(-dplyr::any_of(exclude_taxa)) %>%
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
  
  if (nrow(pca_data) < 3) {
    message(
      "Skipping ", region_use,
      ": insufficient observations"
    )
    return(NULL)
  }
  
  if (ncol(pca_data) < 2) {
    message(
      "Skipping ", region_use,
      ": insufficient taxa"
    )
    return(NULL)
  }
  
  if (any(!is.finite(as.matrix(pca_data)))) {
    message(
      "Skipping ", region_use,
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
        "Skipping ", region_use,
        ": ", e$message
      )
      return(NULL)
    }
  )
  
  if (is.null(pca)) {
    return(NULL)
  }
  
  scores <- as.data.frame(pca$x) %>%
    dplyr::mutate(year = years) %>%
    dplyr::arrange(year) %>%
    dplyr::mutate(
      PC1_next = dplyr::lead(PC1),
      PC2_next = dplyr::lead(PC2)
    )
  
  loadings <- as.data.frame(pca$rotation) %>%
    tibble::rownames_to_column("taxon")
  
  var_explained <- summary(pca)$importance[2, ] * 100
  
  mult <- 5
  
  p <- ggplot2::ggplot(
    scores,
    ggplot2::aes(PC1, PC2)
  ) +
    ggplot2::geom_segment(
      data = dplyr::filter(scores, !is.na(PC1_next)),
      ggplot2::aes(
        x = PC1,
        y = PC2,
        xend = PC1_next,
        yend = PC2_next,
        colour = year
      ),
      linewidth = 1.2,
      arrow = grid::arrow(
        length = unit(0.2, "cm"),
        type = "open"
      )
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
      title = region_use,
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
      caption = paste0(
        "Principal components analysis of Hellinger-transformed annual mean abundances of phytoplankton groups.\n",
        "Taxa appearing in fewer than ", n, " years excluded.\n",
        "Dashed lines indicate PCA loadings. Coloured lines and points show trends over time.\n",
        "Excluded taxa (if any): ",
        paste(exclude_taxa, collapse = ", "), "."
      )
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = 2),
      axis.title = ggplot2::element_text(face = 2)
    )
  
  p
}

safe_make_pca_plot <- purrr::possibly(
  make_pca_plot,
  otherwise = NULL
)

exclude_taxa <- c("NA", "Protozoa")

purrr::walk(
  region_lookup$region,
  function(region_use) {
    
    p <- safe_make_pca_plot(
      region_use = region_use,
      dat = rgn_df_clean,
      n = n
    )
    
    if (!is.null(p)) {
      
      fname <- paste0(
        outdir, "/",
        gsub("[^A-Za-z0-9_]", "_", region_use),
        "_PCA_no_NAProt_rgn.png"
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

# export as a pdf ----
tictoc::tic("Generate region PCA plots in pdf")

pdf(
  "outputs/PCA_plots/All_Region_PCA.pdf",
  width = 10,
  height = 8
)

for(i in seq_len(nrow(region_lookup))) {
  
  p <- safe_make_pca_plot(
    region_use = region_lookup$region[i],
    dat = rgn_df_clean,
    n = n
  )
  
  if(!is.null(p)) {
    print(p)
  }
}

dev.off()

tictoc::toc(log = TRUE)
