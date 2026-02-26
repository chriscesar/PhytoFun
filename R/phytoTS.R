# phytoTS.R ####
## Time series dabblings####
# load packages ####
ld_pkgs <- c("tidyverse","tictoc","mgcv","gratia")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log
tictoc::tic("TOTAL TIME")
theme_set(ggthemes::theme_few())

## load data
tictoc::tic("load data")
df0 <- readRDS("data/Phyto_taxa_long.Rdat")
df <- df0
tictoc::toc(log=TRUE)

# THAMES LOWER
df %>% 
  # remove old data
  dplyr::filter(sample_date > "2009-12-31") %>% 
  # wb of interest
  # dplyr::filter(wb_name == "THAMES LOWER") %>%
  # dplyr::filter(wb_name == "Thames Coastal North") %>%
  # dplyr::filter(wb_name == "Kent South") %>%
  # dplyr::filter(wb_name == "DART") %>%
  # dplyr::filter(str_detect(wb_name, "THAMES")) %>%
  dplyr::filter(str_detect(wb_name, "HUMBER")) %>%
  # dplyr::filter(str_detect(wb_name, "SEVERN")) %>% 
  # dplyr::filter(wb_name == "Mersey Mouth") %>%
  # dplyr::filter(str_detect(wb_name, "SOUTHAMPTON WATER")) %>% 
  # {wb <- unique(pull(., wb_name))[1]; .} %>%
  { wb <<- unique(dplyr::pull(., wb_name))[1]; . } %>% # Assign temporary variable
  # taxa of interest
  dplyr::filter(phylum.y == "Bacillariophyta"|infraphylum == "Dinoflagellata") %>% 
  ## create taxon variable
  dplyr::mutate(taxon = case_when(
    phylum.y == "Bacillariophyta" ~ "Diatom",
    infraphylum == "Dinoflagellata" ~ "Dinoflagellate",
    TRUE ~ "ERROR"
    )) %>% 
  ## sum taxa
  group_by(sample_id, sample_date,taxon) %>% #str()
  summarise(cells_per_litre = sum(cells_per_litre_millilitre,
                                             na.rm = TRUE),
            .groups = "drop") %>% 
    ungroup() %>% 
    dplyr::filter(cells_per_litre != 0) %>% 
  ggplot(., aes(
    x = sample_date,
    y = log10(cells_per_litre),
    colour = taxon,
    )
  ) +
  # geom_vline(xintercept = c(seq(from = min(year(sample_date),
  #                                          to = maz(year(sample_date))),by=1)))+
  geom_point(alpha = .1) +
  geom_smooth(
    method = "gam"
  ) +
  labs(
    title=wb
  ) +
  ggthemes::theme_few()+
  theme(
    plot.title = element_text(face = 2),
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = 2),
    axis.text = element_text(face = 2),
        );rm(wb)

# create data object for model ####
df %>% 
  # remove old data
  dplyr::filter(sample_date > "2009-12-31") %>% 
  # wb of interest
  # dplyr::filter(wb_name == "THAMES LOWER") %>%
  # dplyr::filter(wb_name == "DART") %>%
  # dplyr::filter(str_detect(wb_name, "THAMES")) %>%
  dplyr::filter(str_detect(wb_name, "HUMBER")) %>%
  # dplyr::filter(str_detect(wb_name, "SEVERN")) %>% 
  # dplyr::filter(wb_name == "Mersey Mouth") %>%
  # dplyr::filter(str_detect(wb_name, "SOUTHAMPTON WATER")) %>% 
  # {wb <- unique(pull(., wb_name))[1]; .} %>%
  # taxa of interest
  dplyr::filter(phylum.y == "Bacillariophyta"|infraphylum == "Dinoflagellata") %>% 
  ## create taxon variable
  dplyr::mutate(taxon = case_when(
    phylum.y == "Bacillariophyta" ~ "Diatom",
    infraphylum == "Dinoflagellata" ~ "Dinoflagellate",
    TRUE ~ "ERROR"
  )) %>% 
  ## sum taxa
  group_by(sample_id, sample_date,taxon) %>% #str()
  summarise(cells_per_litre = sum(cells_per_litre_millilitre,
                                  na.rm = TRUE),
            .groups = "drop") %>% 
  ungroup() %>% 
  dplyr::filter(cells_per_litre != 0) %>% 
  as.data.frame(.)-> data_in

origin <- min(data_in$sample_date, na.rm = TRUE)

data_in$days_since <- as.numeric(difftime(data_in$sample_date, origin, units = "days"))
data_in$taxon <- factor(data_in$taxon)
data_in %>% filter(taxon == "Diatom") ->data_in
rm(origin)

# m1 <- mgcv::gam(log10(cells_per_litre) ~ s(days_since, by = taxon, bs = "cr"),

m1 <- mgcv::gam(log10(cells_per_litre) ~ s(days_since, bs = "cr"),
                                family="poisson",
                data = data_in)

### generate model estimates ###
sm <- gratia::smooth_estimates(m1,n = 2000) %>% 
  gratia::add_confint()

terms_date_by <- gratia::smooths(m1)

# Calculate first derivative
gratia::derivatives(
  m1,
  term = terms_date_by,
  order = 1,
  interval = "simultaneous",
  n = 2000,
  type="central",
  unconditional = TRUE
  ) -> deriv1

## bind smooth and derivatives together ####
df_pred <- bind_cols(
  sm %>% rename_with( ~ paste0("sm_", .x)),
  deriv1 %>% rename_with( ~ paste0("d1_", .x))
  ); rm(deriv1,sm)

## calculate changepoints
df_pred <- df_pred %>% 
  arrange(sm_days_since) %>% 
  ## identify direction of change based on derivatives
  mutate(
    sig_pos = d1_.lower_ci >0,
    sig_neg = d1_.upper_ci <0,
    sig_dir = case_when(
      sig_pos ~ "increasing",
      sig_neg ~ "decreasing",
      TRUE ~ "uncertain"
    )
  ) %>% #names()
  ## is this a changepoint?
  mutate(
    direction_prev = lag(sig_dir)) %>% 
  mutate(changepoint = sig_dir != direction_prev &
           sig_dir %in% c("increasing","decreasing"))
rm(terms_date_by)

## extract changepoints
df_pred %>% filter(changepoint == TRUE) %>% 
  arrange(sm_days_since) -> df_pred_ch
