# TaxRichExpl.R ####
## Taxon richness explorations

# load packages ####
ld_pkgs <- c("tidyverse","tictoc","stringr",
             "brms", "tidybayes", "data.table")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

source("R/000setup.R")

# load data ####
tic("Load data")
df0_raw <- readRDS("outputs/Phyto_2000_2025_with_Lifeforms_USE_updateCalcs_updateWB.Rdat")
toc(log=TRUE)

# calculate taxon richness
df0_raw %>% 
  dplyr::select(
    rbd,
    wbid,
    wb_name,
    site_id,
    replicate_code,
    biosys_short,
    sample_date,
    surveillan,
    eastings, northings,
    name_use,
    cells_per_m3
    ) %>% 
  dplyr::mutate(rbd = factor(rbd,
                             levels = c(
                               "Northumbria",
                               "Humber",
                               "Anglian",
                               "Thames",
                               "South East",
                               "South West",
                               "Severn",
                               "North West",
                               "Solway Tweed"
                             )
  )) %>% 
  dplyr::mutate(wb_name = factor(wb_name,
                                 levels = unique(wb_name[order(rbd)])
                                 )
                ) %>% 
  # group and sum
  group_by(across(!cells_per_m3)) %>% 
  summarise(cells_per_m3 = sum(cells_per_m3, na.rm = TRUE),
            .groups = "drop") %>% 
  ### widen chart
  tidyr::pivot_wider(
    names_from = name_use,
    values_from = cells_per_m3,
    values_fill = 0
    ) -> dfw

# calculate taxon richness & Abundance & append to wide data
S <- vegan::specnumber(dfw[,-c(1:10)])  
N <- rowSums(dfw[,-c(1:10)])
dfw$S <- S; rm(S)
dfw$N <- N; rm(N)

# create annual summary for modelling ####
dfw %>% 
  dplyr::filter(sample_date >"2008-01-01") %>% #names() %>% 
  dplyr::select(
    rbd,
    wbid,
    wb_name,
    site_id,
    replicate_code,
    biosys_short,
    sample_date,
    surveillan,
    eastings,
    northings,
    S,
    N
  ) %>% 
  dplyr::mutate(year = lubridate::year(sample_date)) %>% 
  dplyr::select(-c(sample_date,
                   biosys_short,
                   replicate_code,
                   surveillan,
                   eastings,
                   northings,
                   wbid,
                   wb_name,
                   site_id
  )
  ) %>% 
  group_by(across(-c(S, N))) %>% 
  summarise(S = mean(S, na.rm = TRUE),
            N = mean(N, na.rm = TRUE),
            .groups = "drop") %>% ungroup() -> df_S_Annual

# TAXON RICHNESS DATA ####
# throw away WBs with fewer than this number of samples:
number <- 500
dfw %>% 
  dplyr::filter(S > 0) %>% # remove zero S values
  dplyr::filter(sample_date > "2008-01-01") %>% 
  dplyr::group_by(wb_name) %>%
  dplyr::filter(n() >= number) %>% ungroup() %>% 
  ggplot(aes(
    x = sample_date,
    y = S
  )
  ) +
  geom_point()+
  facet_wrap(.~wb_name)+
  geom_smooth(
    method = "gam"
  )

### single WB
# distributional model (models location (mean) and shape (variance))
## fit 1 ####
dfw %>% 
  dplyr::filter(wb_name == "CARRICK ROADS INNER") %>% # remove zero S values
  dplyr::filter(sample_date > "2008-01-01") %>% 
  ggplot(aes(
    x = sample_date,
    y = S
    )
    ) +
  geom_point()+
  geom_smooth(
    method = "gam"
  )

tic("Fit distributional model of S over time")
# fit1 <- brm(bf(S ~ sample_date,
#                sigma ~ sample_date),
#             data = dfw %>% dplyr::filter(wb_name == "CARRICK ROADS INNER") %>% # remove zero S values
#               dplyr::filter(sample_date > "2008-01-01"),
#             family = gaussian())
# saveRDS(fit1, file= "outputs/models/01distrib_S_time.Rdat")
fit1 <- readRDS("outputs/models/01distrib_S_time.Rdat")
toc(log=TRUE)

tic("Model summaries")
summary(fit1)
toc(log=TRUE)
plot(fit1, nvariables = 2, ask = FALSE)

hyp <- c("exp(sigma_Intercept) = 0",
         "exp(sigma_Intercept + sigma_sample_date) = 0")
hypothesis(fit1, hyp)

hyp <- "exp(sigma_Intercept + sigma_sample_date) > exp(sigma_Intercept)"
(hyp <- hypothesis(fit1, hyp))
plot(hyp, chars = NULL)

## fit 2 ####
df_S_Annual %>% 
  ggplot(aes(
    x = year,
    y = S
  )
  ) +
  geom_point() + facet_wrap(. ~ rbd) + geom_smooth(method = "gam")

tic("fit2")
# fit2 <- brm(bf(S ~ s(year) + (1|rbd)),
# fit2 <- brm(bf(S ~ s(year) + rbd),
# fit2 <- brm(bf(S ~ s(year) + rbd + (year|rbd), sigma ~ s(year)),
#             data = df_S_Annual,
#             family = gaussian())
# saveRDS(fit2, file= "outputs/models/02distrib_S_rbd_time.Rdat")
fit2 <- readRDS("outputs/models/02distrib_S_rbd_time.Rdat")
toc(log=TRUE)

summary(fit2)
plot(fit2)

conditional_effects(fit2, surface = TRUE)
conditional_smooths(fit2, surface = TRUE)

####
## fit3 ####
# using info from:
# https://ourcodingclub.github.io/tutorials/brms/#part5

# SEE ALSO:
# https://m-clark.github.io/posts/2021-02-28-practical-bayes-part-ii/#summary-the-practical-approach-to-bayesian-models

tic("fit3")
# fit3 <- brms::brm(S ~ I(year - 2008),
#                            data = df_S_Annual,
#                         family = gaussian(),
#                         chains = 4,
#                         iter = 3000,
#                         warmup = 1000)
# saveRDS(fit3, file = "outputs/models/03brmsSimple.Rdat")
fit3 <- readRDS("outputs/models/03brmsSimple.Rdat")
toc(log=TRUE)

summary(fit3)
pp_check(fit3)

### fit 3.1 ####
## add 'year' as a random effect
tic("fit3.1")
# fit3.1 <- brms::brm(S ~ I(year - 2008) + (1|year),
#                   data = df_S_Annual,
#                   family = gaussian(),
#                   chains = 4,
#                   iter = 3000,
#                   warmup = 1000)
# saveRDS(fit3.1, file = "outputs/models/03.1brmsSimple.Rdat")
fit3.1 <- readRDS("outputs/models/03.1brmsSimple.Rdat")
toc(log=TRUE)

summary(fit3.1);plot(fit3.1)

### fit 3.2 ####
tic("fit3.2") #adding location parameter
# fit3.2 <- brms::brm(S ~ I(year - 2008)+rbd,
#                     data = df_S_Annual,
#                     family = gaussian(),
#                     chains = 4,
#                     iter = 3000,
#                     warmup = 1000)
# saveRDS(fit3.2, file = "outputs/models/03.2brmsSimple.Rdat")
fit3.2 <- readRDS("outputs/models/03.2brmsSimple.Rdat")
toc(log = TRUE)

summary(fit3.2)
plot(fit3.2)
pp_check(fit3.2,ndraws = 100)

### fit 3.3 ####
tic("fit3.3") #adding variable sigma
# fit3.3 <- brms::brm(bf(S ~ I(year - 2008)+rbd,
#                     sigma ~ rbd),
#                     data = df_S_Annual,
#                     family = gaussian(),
#                     chains = 4,
#                     iter = 3000,
#                     warmup = 1000)
# saveRDS(fit3.3, file = "outputs/models/03.3brmsSimple.Rdat")
fit3.3 <- readRDS("outputs/models/03.3brmsSimple.Rdat")
toc(log = TRUE)

tic("fit3.3 summaries")
summary(fit3.3)
plot(fit3.3)
pp_check(fit3.3,ndraws = 100)

df_S_Annual %>% 
  tidybayes::add_predicted_draws(fit3.3) %>% #add posterior distribution
  ggplot(
    aes(x = year,
        y = S)
  ) +
  tidybayes::stat_lineribbon(aes(# regression line and CI)
    y = .prediction),
    .width = c(.95, .80, .50),
    alpha=0.4,
    colour="blue",
    show.legend = FALSE
  )+
  scale_fill_brewer(palette = "Greys")+
  geom_point(data = df_S_Annual, colour = "black", size = 2) +
  facet_wrap(.~rbd) +
  labs(
    title = "Mean taxon richness over time by river basin district",
    y = "Taxon richness"
  )+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2),
    strip.text = element_text(face=2)
  )
toc(log=TRUE)

### fit 3.4 ####
tic("fit3.4")
# fit3.4 <- brms::brm(bf(S ~ s(I(year - 2008))+rbd,
#                        sigma ~ rbd),
#                     data = df_S_Annual,
#                     family = gaussian(),
#                     chains = 4,
#                     iter = 3000,
#                     warmup = 1000)
# saveRDS(fit3.4, file = "outputs/models/03.4brmsSimple.Rdat")
fit3.4 <- readRDS("outputs/models/03.4brmsSimple.Rdat")
toc(log=TRUE)

tic("fit3.4 summaries")
summary(fit3.4)
plot(fit3.4)
pp_check(fit3.4,ndraws = 100)

df_S_Annual %>% 
  dplyr::group_by(rbd) %>% 
  tidybayes::add_predicted_draws(fit3.4) %>% #add posterior distribution
  ggplot(
    aes(x = year,
        y = S)
  ) +
  tidybayes::stat_lineribbon(aes(# regression line and CI)
    y = .prediction),
    .width = c(.95, .80, .50),
    alpha=0.4,
    colour="blue",
    show.legend = FALSE
  )+
  scale_fill_brewer(palette = "Greys")+
  geom_point(data = df_S_Annual, colour = "black", size = 2) +
  facet_wrap(.~rbd) +
  labs(
    title = "Mean taxon richness over time by river basin district",
    y = "Taxon richness",
    caption = "Blue lines indicate model prediction, shaded bands represent 50, 80, and 95% credible intervals"
  )+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2),
    strip.text = element_text(face=2),
    axis.text = element_text(face=2)
  )
toc(log=TRUE)

###
# FUN
df_S_Annual %>%
  group_by(rbd) %>%
  add_predicted_draws(fit3.4) %>%
  ggplot(aes(x = year, y = S,
             # color = ordered(rbd),
             # fill = ordered(rbd),
             # shape = rbd
             )) +
  stat_lineribbon(
    aes(y = .prediction),
    .width = c(.95, .80, .50),
    alpha = 0.3,
    show.legend = FALSE
    ) +
  scale_fill_brewer(palette = "Greys") +
  # scale_color_brewer(palette = "PiYG") +
  # scale_shape_manual(values = c(0:3,0:3,0))+
  geom_point(data = df_S_Annual) +
  labs(
    title = "Mean taxon richness over time by river basin district",
    y = "Taxon richness",
    caption = "Blue lines indicate model prediction, shaded bands represent 50, 80, and 95% credible intervals"
    )+
  facet_wrap(. ~ rbd)+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2),
    strip.text = element_text(face=2),
    axis.text = element_text(face=2)
  )

df_S_Annual %>%
  group_by(rbd) %>%
  add_predicted_draws(fit3.4) %>%
  ggplot(aes(x = year, y = S,
             color = ordered(rbd),
             linetype = ordered(rbd)
             )) +
  stat_lineribbon(
    aes(y = .prediction),
    .width = 0,
    #.width = c(.95, .80, .50),
    # alpha = 0.3,
    #show.legend = FALSE
  ) +
  #scale_fill_brewer(palette="Greys")+
  scale_color_brewer(palette = "Paired")+
  scale_linetype_manual(values = c(
    "solid", "dashed", "dotted", "dotdash",
    "longdash", "twodash", "solid", "dashed", "dotted"
  ))+
  guides(fill = "none")+
  labs(
    title = "Model-predicted taxon richness over time by river basin district",
    y = "Taxon richness",
    caption = "Lines indicate model predictions",
    linetype = "River basin district",
    colour = "River basin district"
  )+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2),
    strip.text = element_text(face=2),
    axis.text = element_text(face=2),
    legend.title = element_text(face=2)
  )

brms::pp_check(fit3.4, ndraws = 200)

### fit 3.5 ####
tic("fit3.5") #truncated gaussian
# fit3.5 <- brms::brm(bf(S|trunc(lb = 0) ~ s(I(year - 2008))+rbd,
#                        sigma ~ rbd),
#                     data = df_S_Annual,
#                     family = gaussian(),
#                     chains = 4,
#                     iter = 10000,
#                     warmup = 3000)
# saveRDS(fit3.5, file = "outputs/models/03.5brmsSimple.Rdat")
fit3.5 <- readRDS("outputs/models/03.5brmsSimple.Rdat")
toc(log=TRUE)

tic("fit3.5 summaries")
(summ_fit3.5 <- summary(fit3.5))
pp_check(fit3.5,
         # ndraws = 800,
         ndraws = 300,
         alpha=0.15
         )
toc(log=TRUE)

tic("fit3.5 plots")
df_S_Annual %>%
  group_by(rbd) %>%
  add_predicted_draws(fit3.5) %>%
  ggplot(aes(x = year, y = S,
             # color = ordered(rbd),
             # fill = ordered(rbd),
             # shape = rbd
  )) +
  stat_lineribbon(
    aes(y = .prediction),
    .width = c(.95, .80, .50),
    alpha = 0.3,
    show.legend = FALSE
  ) +
  scale_fill_brewer(palette = "Greys") +
  geom_point(data = df_S_Annual) +
  labs(
    title = "Mean taxon richness over time by river basin district",
    y = "Taxon richness",
    caption = "Bold lines indicate model prediction, shaded bands represent 50, 80, and 95% credible intervals"
    )+
  facet_wrap(. ~ rbd)+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2),
    strip.text = element_text(face=2),
    axis.text = element_text(face=2)
  )

df_S_Annual %>%
  group_by(rbd) %>%
  add_predicted_draws(fit3.5) %>%
  ggplot(aes(x = year, y = S,
             color = ordered(rbd),
             linetype = ordered(rbd)
  )) +
  stat_lineribbon(
    aes(y = .prediction),
    .width = 0,
    #.width = c(.95, .80, .50),
    # alpha = 0.3,
    #show.legend = FALSE
  ) +
  #scale_fill_brewer(palette="Greys")+
  scale_color_brewer(palette = "Paired")+
  scale_linetype_manual(values = c(
    "solid", "dashed", "dotted", "dotdash",
    "longdash", "twodash", "solid", "dashed", "dotted"
  ))+
  guides(fill = "none")+
  ylim(0,NA)+
  labs(
    title = "Model-predicted taxon richness over time by river basin district",
    y = "Taxon richness",
    caption = paste0("Lines indicate model predictions\n",
                     "Model formula: ",fit3.5$formula,"\n",
                     "Family: ",fit3.5$family,"\n",
                     "Chains: ",summ_fit3.5$chains,"; Iterations: ",summ_fit3.5$iter,"\n",
                     "R-sq: ",round(bayes_R2(fit3.5)[,"Estimate"]*100,2),"%"),
    linetype = "River basin district",
    colour = "River basin district"
  )+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2),
    strip.text = element_text(face=2),
    axis.text = element_text(face=2),
    legend.title = element_text(face=2)
  )
toc(log=TRUE)

#### fit3.5 Plot second derivatives ####
tic("Extract posterior smooths")

##### Extract posterior smooths ####
year_seq <- seq(min(df_S_Annual$year), max(df_S_Annual$year), length.out = 100)

newdata <- expand_grid(
  year = year_seq,
  rbd = unique(df_S_Annual$rbd)
)

fitted_draws <- fit3.5 |>
  add_fitted_draws(newdata = newdata, re_formula = NA)
toc(log = TRUE)

tic("Compute derivatives efficiently")
##### Compute derivatives efficiently ####

################
# convert fitted_draws to data table
tic("Convert to DT")
setDT(fitted_draws)
setorder(fitted_draws, rbd, .draw, year)
fitted_draws[
  , first_deriv := c(NA_real_, diff(.value) / diff(year)), 
  by = .(rbd, .draw)
]
fitted_draws[
  , second_deriv := c(NA_real_, diff(first_deriv) / diff(year[-1])), 
  by = .(rbd, .draw)
]
fitted_draws <- fitted_draws[.draw %% 10 == 0]  # keep 1 in 10 draws
toc(log=TRUE)

###############
get_derivatives <- function(df) {
  df |>
    arrange(year) |>
    group_by(.draw) |>
    mutate(
      first_deriv  = c(NA, diff(.value) / diff(year)),
      second_deriv = c(NA, diff(first_deriv) / diff(year[-1]))
    ) |>
    ungroup()
}
tic("calc derivatives")
deriv_df <- fitted_draws |>
  group_by(rbd) |>
  group_modify(~ get_derivatives(.x))
toc(log=TRUE)
toc(log=TRUE)

#### Summarise posterior for curvature ####
tic("Summarise posterior for curvature")
curvature_summary <- deriv_df |>
  group_by(rbd, year) |>
  summarise(
    mean_curv = mean(second_deriv, na.rm = TRUE),
    lower     = quantile(second_deriv, 0.025, na.rm = TRUE),
    upper     = quantile(second_deriv, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

toc(log = TRUE)

tic("Identify changepoints")
#### Identify changepoints ####
changepoints <- curvature_summary |>
  group_by(rbd) |>
  filter(sign(mean_curv) != lag(sign(mean_curv), default = first(sign(mean_curv)))) |>
  filter(!between(0, lower, upper))
toc(log=TRUE)

tic("Plot for diagnosis")
#### Plot for diagnosis ####
ggplot(curvature_summary, aes(year, mean_curv)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ rbd)
toc(log=TRUE)


### fit 3.6 ####
tic("fit3.6") # add factor-smoothed term to allow different trends by rbd
# fit3.6 <- brms::brm(bf(S|trunc(lb = 0) ~ s(I(year - 2008),
#                                            rbd,
#                                            bs = "fs"),
#                        sigma ~ rbd),
#                     data = df_S_Annual,
#                     family = gaussian(),
#                     chains = 4,
#                     iter = 10000,
#                     warmup = 3000)
# saveRDS(fit3.6, file = "outputs/models/03.6brmsSimple.Rdat")
fit3.6 <- readRDS("outputs/models/03.6brmsSimple.Rdat")
toc(log=TRUE)

tic("fit3.6 summaries")
(summ_fit3.6 <- summary(fit3.6))
pp_check(fit3.6,
         # ndraws = 800,
         ndraws = 300,
         alpha=0.25
         )

toc(log=TRUE)

tic("fit3.6 plots")
df_S_Annual %>%
  group_by(rbd) %>%
  add_predicted_draws(fit3.6) %>%
  ggplot(aes(x = year, y = S,
             # color = ordered(rbd),
             # fill = ordered(rbd),
             # shape = rbd
  )) +
  stat_lineribbon(
    aes(y = .prediction),
    .width = c(.95, .80, .50),
    alpha = 0.3,
    show.legend = FALSE
  ) +
  scale_fill_brewer(palette = "Greys") +
  geom_point(data = df_S_Annual) +
  labs(
    title = "Mean taxon richness over time by river basin district",
    y = "Taxon richness",
    caption = paste0("Bold lines indicate model prediction, shaded bands represent 50, 80, and 95% credible intervals\n",
                     "Model formula: (",fit3.6$formula,", ",as.character(fit3.6$formula)[2],")\n",
                     "Family: ",fit3.6$family,"\n",
                     "Chains: ",summ_fit3.6$chains,"; Iterations: ",summ_fit3.6$iter,"\n",
                     "R-sq: ",round(bayes_R2(fit3.6)[,"Estimate"]*100,2),"%")
  )+
  facet_wrap(. ~ rbd)+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2),
    strip.text = element_text(face=2),
    axis.text = element_text(face=2)
  )

df_S_Annual %>%
  group_by(rbd) %>%
  add_predicted_draws(fit3.6) %>%
  ggplot(aes(x = year, y = S,
             color = ordered(rbd),
             linetype = ordered(rbd)
  )) +
  stat_lineribbon(
    aes(y = .prediction),
    .width = 0,
    #.width = c(.95, .80, .50),
    # alpha = 0.3,
    #show.legend = FALSE
  ) +
  #scale_fill_brewer(palette="Greys")+
  scale_color_brewer(palette = "Paired")+
  scale_linetype_manual(values = c(
    "solid", "dashed", "dotdash", "dotted",
    "longdash", "twodash", "solid", "dashed", "dotted"
  ))+
  guides(fill = "none")+
  ylim(0,NA)+
  labs(
    title = "Model-predicted taxon richness over time by river basin district",
    y = "Taxon richness",
    caption = paste0("Lines indicate model predictions\n",
                     "Model formula: (",fit3.6$formula,", ",as.character(fit3.6$formula)[2],")\n",
                     "Family: ",fit3.6$family,"\n",
                     "Chains: ",summ_fit3.6$chains,"; Iterations: ",summ_fit3.6$iter,"\n",
                     "R-sq: ",round(bayes_R2(fit3.6)[,"Estimate"]*100,2),"%"),
    linetype = "River basin district",
    colour = "River basin district"
  )+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2),
    strip.text = element_text(face=2),
    axis.text = element_text(face=2),
    legend.title = element_text(face=2)
  )
toc(log=TRUE)

#### fit3.6 Plot second derivatives ####
tic("Extract posterior smooths")

##### Extract posterior smooths ####
year_seq <- seq(min(df_S_Annual$year), max(df_S_Annual$year), length.out = 100)

newdata <- expand_grid(
  year = year_seq,
  rbd = unique(df_S_Annual$rbd)
)

fitted_draws <- fit3.6 |>
  add_fitted_draws(newdata = newdata, re_formula = NA)
toc(log = TRUE)

tic("Compute derivatives efficiently")
##### Compute derivatives efficiently ####

# convert fitted_draws to data table
tic("Convert to DT")
setDT(fitted_draws)
setorder(fitted_draws, rbd, .draw, year)
fitted_draws[
  , first_deriv := c(NA_real_, diff(.value) / diff(year)), 
  by = .(rbd, .draw)
]
fitted_draws[
  , second_deriv := c(NA_real_, diff(first_deriv) / diff(year[-1])), 
  by = .(rbd, .draw)
]
fitted_draws <- fitted_draws[.draw %% 10 == 0]  # keep 1 in 10 draws
toc(log=TRUE)

get_derivatives <- function(df) {
  df |>
    arrange(year) |>
    group_by(.draw) |>
    mutate(
      first_deriv  = c(NA, diff(.value) / diff(year)),
      second_deriv = c(NA, diff(first_deriv) / diff(year[-1]))
    ) |>
    ungroup()
}
tic("calc derivatives")
deriv_df <- fitted_draws |>
  group_by(rbd) |>
  group_modify(~ get_derivatives(.x))
toc(log=TRUE)
toc(log=TRUE)

#### Summarise posterior for curvature ####
tic("Summarise posterior for curvature")
curvature_summary <- deriv_df |>
  group_by(rbd, year) |>
  summarise(
    # mean_curv = mean(first_deriv, na.rm = TRUE),
    # lower     = quantile(first_deriv, 0.025, na.rm = TRUE),
    # upper     = quantile(first_deriv, 0.975, na.rm = TRUE),
    mean_curv = mean(second_deriv, na.rm = TRUE),
    lower     = quantile(second_deriv, 0.025, na.rm = TRUE),
    upper     = quantile(second_deriv, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

toc(log = TRUE)

tic("Identify changepoints")
#### Identify changepoints ####
changepoints <- curvature_summary |>
  group_by(rbd) |>
  filter(sign(mean_curv) != lag(sign(mean_curv), default = first(sign(mean_curv)))) |>
  filter(!between(0, lower, upper))
toc(log=TRUE)

tic("Plot for diagnosis")

#### Plot for diagnosis ####
ggplot(curvature_summary, aes(year, mean_curv)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    ### first derivatives
    title = "First derivatives of taxon richness over time",
    #first derivatives
    subtitle = "Visualisation of rates of change over time.\nPositive first derivatives indicate increasing taxon richness, negative values indicate decreasing taxon richness"
    ## second derivatives
    # title = "Second derivatives of taxon richness over time",
    # ## second derivatives
    # subtitle = "Visualisation of rates of change over time.\nPositive second derivatives indicate accelerating change, negative values indicate decelerating change, and zero values indicate a constant rate of change"
    )+
  facet_wrap(~ rbd)
toc(log=TRUE)


# compare models ####
# Model with elpd of 0 shows the best fit to our data
loo(fit3,fit3.1,fit3.2,fit3.3,fit3.4, fit3.5,fit3.6,
    compare = TRUE)
# fit3.5 gives best fit

x0 <- as.data.frame(bayes_R2(fit3))
x1 <- as.data.frame(bayes_R2(fit3.1))
x2 <- as.data.frame(bayes_R2(fit3.2))
x3 <- as.data.frame(bayes_R2(fit3.3))
x4 <- as.data.frame(bayes_R2(fit3.4))
x5 <- as.data.frame(bayes_R2(fit3.5))
x6 <- as.data.frame(bayes_R2(fit3.6))

row.names(x0) <- "fit3"
row.names(x1) <- "fit3.1"
row.names(x2) <- "fit3.2"
row.names(x3) <- "fit3.3"
row.names(x4) <- "fit3.4"
row.names(x5) <- "fit3.5"
row.names(x6) <- "fit3.6"
(fit3x_R2 <- rbind(x0,x1,x2,x3,x4,x5,x6))
rm(x0,x1,x2,x3,x4,x5,x6)

# Model 3.5 seems a good enough fit.
# identify changepoints for fit 3.5 ####

############
# DENSITY DATA ####
# throw away WBs with fewer than this number of samples:
number <- 500
dfw %>% 
  dplyr::filter(N > 0) %>% # remove zero S values
  dplyr::filter(sample_date > "2008-01-01") %>% 
  dplyr::group_by(wb_name) %>%
  dplyr::filter(n() >= number) %>% ungroup() %>% 
  ggplot(aes(
    x = sample_date,
    y = log10(N)
  )
  ) +
  geom_point()+
  facet_wrap(.~wb_name)+
  geom_smooth(
    method = "gam"
  )

### fit 3.5N ####
tic("fit3.5N") #truncated gaussian
# fit3.5N <- brms::brm(bf(log(N)|trunc(lb = 0) ~ s(I(year - 2008))+rbd,
#                        sigma ~ rbd),
#                     data = df_S_Annual,
#                     family = gaussian(),
#                     chains = 4,
#                     iter = 10000,
#                     warmup = 3000)
# saveRDS(fit3.5N, file = "outputs/models/03.5brmsSimple_N.Rdat")
fit3.5N <- readRDS("outputs/models/03.5brmsSimple_N.Rdat")
toc(log=TRUE)

tic("fit3.5 summaries")
(summ_fit3.5N <- summary(fit3.5N))
pp_check(fit3.5N,
         # ndraws = 800,
         ndraws = 300,
         alpha=0.15
)
toc(log=TRUE)

tic("fit3.5N plots")
df_S_Annual %>%
  group_by(rbd) %>%
  add_predicted_draws(fit3.5N) %>%
  ggplot(aes(x = year, y = log(N),
             # color = ordered(rbd),
             # fill = ordered(rbd),
             # shape = rbd
  )) +
  stat_lineribbon(
    aes(y = .prediction),
    .width = c(.95, .80, .50),
    alpha = 0.3,
    show.legend = FALSE
  ) +
  scale_fill_brewer(palette = "Greys") +
  geom_point(data = df_S_Annual) +
  labs(
    title = "Mean log(phytoplankton density) over time by river basin district",
    y = "Log10(Phytoplankton density)",
    caption = "Bold lines indicate model prediction, shaded bands represent 50, 80, and 95% credible intervals"
  )+
  facet_wrap(. ~ rbd)+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2),
    strip.text = element_text(face=2),
    axis.text = element_text(face=2)
  )

df_S_Annual %>%
  group_by(rbd) %>%
  add_predicted_draws(fit3.5N) %>%
  ggplot(aes(x = year, y = log(N),
             color = ordered(rbd),
             linetype = ordered(rbd)
  )) +
  stat_lineribbon(
    aes(y = .prediction),
    .width = 0,
    #.width = c(.95, .80, .50),
    # alpha = 0.3,
    #show.legend = FALSE
  ) +
  #scale_fill_brewer(palette="Greys")+
  scale_color_brewer(palette = "Paired")+
  scale_linetype_manual(values = c(
    "solid", "dashed", "dotted", "dotdash",
    "longdash", "twodash", "solid", "dashed", "dotted"
  ))+
  guides(fill = "none")+
  # ylim(0,NA)+
  labs(
    title = "Model-predicted log phytoplankton density over time by river basin district",
    y = "Log10(Phytoplankton density)",
    caption = paste0("Lines indicate model predictions\n",
                     "Model formula: ",fit3.5$formula,"\n",
                     "Family: ",fit3.5$family,"\n",
                     "Chains: ",summ_fit3.5N$chains,"; Iterations: ",summ_fit3.5N$iter,"\n",
                     "R-sq: ",round(bayes_R2(fit3.5)[,"Estimate"]*100,2),"%"),
    linetype = "River basin district",
    colour = "River basin district"
  )+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2),
    strip.text = element_text(face=2),
    axis.text = element_text(face=2),
    legend.title = element_text(face=2)
  )
toc(log=TRUE)

#### fit3.5N Plot second derivatives ####
tic("Extract posterior smooths")

##### Extract posterior smooths ####
year_seq <- seq(min(df_S_Annual$year), max(df_S_Annual$year), length.out = 100)

newdata <- expand_grid(
  year = year_seq,
  rbd = unique(df_S_Annual$rbd)
)

fitted_draws <- fit3.5N |>
  add_fitted_draws(newdata = newdata, re_formula = NA)
toc(log = TRUE)

tic("Compute derivatives efficiently")
##### Compute derivatives efficiently ####

################
# convert fitted_draws to data table
tic("Convert to DT")
setDT(fitted_draws)
setorder(fitted_draws, rbd, .draw, year)
fitted_draws[
  , first_deriv := c(NA_real_, diff(.value) / diff(year)), 
  by = .(rbd, .draw)
]
fitted_draws[
  , second_deriv := c(NA_real_, diff(first_deriv) / diff(year[-1])), 
  by = .(rbd, .draw)
]
fitted_draws <- fitted_draws[.draw %% 10 == 0]  # keep 1 in 10 draws
toc(log=TRUE)

###############
get_derivatives <- function(df) {
  df |>
    arrange(year) |>
    group_by(.draw) |>
    mutate(
      first_deriv  = c(NA, diff(.value) / diff(year)),
      second_deriv = c(NA, diff(first_deriv) / diff(year[-1]))
    ) |>
    ungroup()
}
tic("calc derivatives")
deriv_df <- fitted_draws |>
  group_by(rbd) |>
  group_modify(~ get_derivatives(.x))
toc(log=TRUE)
toc(log=TRUE)

#### Summarise posterior for curvature ####
tic("Summarise posterior for curvature")
curvature_summary <- deriv_df |>
  group_by(rbd, year) |>
  summarise(
    mean_curv = mean(second_deriv, na.rm = TRUE),
    lower     = quantile(second_deriv, 0.025, na.rm = TRUE),
    upper     = quantile(second_deriv, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

toc(log = TRUE)

tic("Identify changepoints")
#### Identify changepoints ####
changepoints <- curvature_summary |>
  group_by(rbd) |>
  filter(sign(mean_curv) != lag(sign(mean_curv), default = first(sign(mean_curv)))) |>
  filter(!between(0, lower, upper))
toc(log=TRUE)

tic("Plot for diagnosis")
#### Plot for diagnosis ####
ggplot(curvature_summary, aes(year, mean_curv)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ rbd)
toc(log=TRUE)

##################
### fit 3.6N ####
tic("fit3.6N") # add factor-smoothed term to allow different trends by rbd
# fit3.6N <- brms::brm(bf(log(N)|trunc(lb = 0) ~ s(I(year - 2008),
#                                            rbd,
#                                            bs = "fs"),
#                        sigma ~ rbd),
#                     data = df_S_Annual,
#                     family = gaussian(),
#                     chains = 4,
#                     iter = 10000,
#                     warmup = 3000)
# saveRDS(fit3.6N, file = "outputs/models/03.6brmsSimple_N.Rdat")
fit3.6N <- readRDS("outputs/models/03.6brmsSimple_N.Rdat")
toc(log=TRUE)

tic("fit3.6N summaries")
(summ_fit3.6N <- summary(fit3.6N))
pp_check(fit3.6N,
         # ndraws = 800,
         ndraws = 300,
         alpha=0.25
)

toc(log=TRUE)

tic("fit3.6N plots")
df_S_Annual %>%
  group_by(rbd) %>%
  add_predicted_draws(fit3.6N) %>%
  ggplot(aes(x = year, y = log(N),
             # color = ordered(rbd),
             # fill = ordered(rbd),
             # shape = rbd
  )) +
  stat_lineribbon(
    aes(y = .prediction),
    .width = c(.95, .80, .50),
    alpha = 0.3,
    show.legend = FALSE
  ) +
  scale_fill_brewer(palette = "Greys") +
  geom_point(data = df_S_Annual) +
  labs(
    title = "Mean log phytoplankton density over time by river basin district",
    y = "Log(Phytoplankton density)",
    caption = paste0("Bold lines indicate model prediction, shaded bands represent 50, 80, and 95% credible intervals\n",
                     "Model formula: (",fit3.6N$formula,", ",as.character(fit3.6N$formula)[2],")\n",
                     "Family: ",fit3.6N$family,"\n",
                     "Chains: ",summ_fit3.6N$chains,"; Iterations: ",summ_fit3.6N$iter,"\n",
                     "R-sq: ",round(bayes_R2(fit3.6N)[,"Estimate"]*100,2),"%")
  )+
  facet_wrap(. ~ rbd)+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2),
    strip.text = element_text(face=2),
    axis.text = element_text(face=2)
  )

df_S_Annual %>%
  group_by(rbd) %>%
  add_predicted_draws(fit3.6N) %>%
  ggplot(aes(x = year, y = log(N),
             color = ordered(rbd),
             linetype = ordered(rbd)
  )) +
  stat_lineribbon(
    aes(y = .prediction),
    .width = 0,
    #.width = c(.95, .80, .50),
    # alpha = 0.3,
    #show.legend = FALSE
  ) +
  #scale_fill_brewer(palette="Greys")+
  scale_color_brewer(palette = "Paired")+
  scale_linetype_manual(values = c(
    "solid", "dashed", "dotdash", "dotted",
    "longdash", "twodash", "solid", "dashed", "dotted"
  ))+
  guides(fill = "none")+
  # ylim(0,NA)+
  labs(
    title = "Model-predicted phytoplankton density over time by river basin district",
    y = "Log(phytoplankton density)",
    caption = paste0("Lines indicate model predictions\n",
                     "Model formula: (",fit3.6N$formula,", ",as.character(fit3.6N$formula)[2],")\n",
                     "Family: ",fit3.6N$family,"\n",
                     "Chains: ",summ_fit3.6N$chains,"; Iterations: ",summ_fit3.6N$iter,"\n",
                     "R-sq: ",round(bayes_R2(fit3.6N)[,"Estimate"]*100,2),"%"),
    linetype = "River basin district",
    colour = "River basin district"
  )+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2),
    strip.text = element_text(face=2),
    axis.text = element_text(face=2),
    legend.title = element_text(face=2)
  )
toc(log=TRUE)

#### fit3.6 Plot second derivatives ####
tic("Extract posterior smooths")

##### Extract posterior smooths ####
year_seq <- seq(min(df_S_Annual$year), max(df_S_Annual$year), length.out = 100)

newdata <- expand_grid(
  year = year_seq,
  rbd = unique(df_S_Annual$rbd)
)

fitted_draws <- fit3.6N |>
  add_fitted_draws(newdata = newdata, re_formula = NA)
toc(log = TRUE)

tic("Compute derivatives efficiently")
##### Compute derivatives efficiently ####

# convert fitted_draws to data table
tic("Convert to DT")
setDT(fitted_draws)
setorder(fitted_draws, rbd, .draw, year)
fitted_draws[
  , first_deriv := c(NA_real_, diff(.value) / diff(year)), 
  by = .(rbd, .draw)
]
fitted_draws[
  , second_deriv := c(NA_real_, diff(first_deriv) / diff(year[-1])), 
  by = .(rbd, .draw)
]
fitted_draws <- fitted_draws[.draw %% 10 == 0]  # keep 1 in 10 draws
toc(log=TRUE)

get_derivatives <- function(df) {
  df |>
    arrange(year) |>
    group_by(.draw) |>
    mutate(
      first_deriv  = c(NA, diff(.value) / diff(year)),
      second_deriv = c(NA, diff(first_deriv) / diff(year[-1]))
    ) |>
    ungroup()
}
tic("calc derivatives")
deriv_df <- fitted_draws |>
  group_by(rbd) |>
  group_modify(~ get_derivatives(.x))
toc(log=TRUE)
toc(log=TRUE)

#### Summarise posterior for curvature ####
tic("Summarise posterior for curvature")
curvature_summary <- deriv_df |>
  group_by(rbd, year) |>
  summarise(
    mean_curv = mean(first_deriv, na.rm = TRUE),
    lower     = quantile(first_deriv, 0.025, na.rm = TRUE),
    upper     = quantile(first_deriv, 0.975, na.rm = TRUE),
    # mean_curv = mean(second_deriv, na.rm = TRUE),
    # lower     = quantile(second_deriv, 0.025, na.rm = TRUE),
    # upper     = quantile(second_deriv, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

toc(log = TRUE)

tic("Identify changepoints")
#### Identify changepoints ####
changepoints <- curvature_summary |>
  group_by(rbd) |>
  filter(sign(mean_curv) != lag(sign(mean_curv), default = first(sign(mean_curv)))) |>
  filter(!between(0, lower, upper))
toc(log=TRUE)

tic("Plot for diagnosis")

#### Plot for diagnosis ####
ggplot(curvature_summary, aes(year, mean_curv)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    ### first derivatives
    title = "First derivatives of taxon richness over time",
    #first derivatives
    subtitle = "Visualisation of rates of change over time.\nPositive first derivatives indicate increasing taxon richness, negative values indicate decreasing taxon richness"
    ## second derivatives
    # title = "Second derivatives of taxon richness over time",
    # ## second derivatives
    # subtitle = "Visualisation of rates of change over time.\nPositive second derivatives indicate accelerating change, negative values indicate decelerating change, and zero values indicate a constant rate of change"
  )+
  facet_wrap(~ rbd)
toc(log=TRUE)














unlist(tictoc::tic.log())


