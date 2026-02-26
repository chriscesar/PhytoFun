# 202602_phyto_figs_for_rep.R ####
## generate figures for 2026 Inshore Plankton Report

# load packages ####
ld_pkgs <- c("tidyverse","ggthemes")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log
source("R/000setup.R")

outfol <- "//prodds.ntnl/Shared/AN/KFH/Groups/N_Marine/02 Projects_Tasks/04 E&B/Ecology and Ecosystems/Planktonic_Food_Webs/202602_FigsForReport/"
tictoc::tic("Load and prepare data")
df0 <- readRDS("outputs/Phyto_2000_2025_with_Lifeforms_USE_updateCalcs_updateWB.Rdat")

# generate factor levels
fac_tmp <- sprintf("%04d%02d", rep(2000:2025, each = 12), 1:12)

## keep only 'interesting' variables
df <- df0 %>% 
  select(sample_id,rbd,wb_name,sample_date,phyto_lf,name_use,
         cells_per_litre:md_carbon_tot_mg_c_per_m3
         ) %>% 
  # create yyymm variable & move it
  mutate(yyyymm = sprintf("%d%02d",
                          year(sample_date),
                          month(sample_date))) %>% 
  relocate(yyyymm, .after = sample_date) %>% 
  mutate(date_plot = make_date(year = year(sample_date),
                               month = month(sample_date),
                               day = 1)) %>% 
  relocate(date_plot, .after=yyyymm)

# assign dates to factors
df$yyyymm_fac <- factor(df$yyyymm, levels = fac_tmp)
rm(fac_tmp)
# 
# ### make a 'complete' version so that 'missing' surveys will be included in the charts
# df_complete <- %>%
#   group_by(BIOSYS.Code) %>%  # Group by site
#   # tidyr::complete(yyyy_mm, fill = list(mn_carbTot_m3 = 0)) %>% # Fill missing months with 0 or NA
#   # tidyr::complete(DJF, fill = list(mn_carbTot_m3 = 0)) #Fill missing seasons
#   # tidyr::complete(yyyy_mm, DJF, fill = list(mn_carbTot_m3 = 0))
#   tidyr::complete(yyyy_mm_fac, DJF, fill = list(mn_carbTot_m3 = 0))


tictoc::toc()

# summarise counts by sampling event & plot ####
# png(file = paste0(outfol,"phyto_abund_log10.png"),
png(file = paste0(outfol,"phyto_abund_raw.png"),
    width=18*ppi, height=12*ppi, res=ppi)
df %>%
  ## total abundance: cells per litre
  select(sample_id,date_plot,cells_per_m3) %>% #View()
  group_by(across(-cells_per_m3)) %>% 
  # get sum per sample
  summarise(cells_per_m3 = sum(cells_per_m3,na.rm = TRUE),
            .groups = "drop") %>% 
  # get mean sum by date
  ## drop sample ID
  select(-sample_id) %>% group_by(date_plot) %>% 
  summarise(mn_cells_per_m3 = mean(cells_per_m3,na.rm = TRUE),
            sd_cells_per_m3 = sd(cells_per_m3,na.rm=TRUE),
            n = n(),
            .groups = "drop"
            ) %>% 
  ggplot(.,
         aes(
           x = date_plot,
           # y = cells_per_litre
           # y = log10(mn_cells_per_m3+1)
           y=mn_cells_per_m3
         )
         )+
  geom_rug(sides = "b")+
  geom_vline(xintercept = as.Date(c("2000-01-01","2005-01-01",
                                    "2010-01-01","2015-01-01",
                                    "2020-01-01","2025-01-01"
                                    )
                                  ),
  lty = 2, col="grey")+
  geom_vline(xintercept = as.Date(c("2001-01-01","2002-01-01","2003-01-01","2004-01-01",
                                    "2006-01-01","2007-01-01","2008-01-01","2009-01-01",
                                    "2011-01-01","2012-01-01","2013-01-01","2014-01-01",
                                    "2016-01-01","2017-01-01","2018-01-01","2019-01-01",
                                    "2021-01-01","2022-01-01","2023-01-01","2024-01-01"
                                    )
                                  ),
             lty = 3, col="grey")+
  geom_point(aes(
    fill = n
    ),
    size=4,pch=21)+
  geom_smooth(method = "gam")+
  labs(
    # title = "Log10 mean phytoplanton abundances by year_month",
    title = "Mean phytoplanton abundances by year_month",
    # y = "Log10(n+1) mean cell abundances (cells per m3) per sample",
    y = "Mean cell abundances (cells per m3) per sample",
    caption=paste0("Values represent mean total cell counts across all samples gathered in a particular month","<br>",
                   "Point shading reflects number of samples contributing to that mean",
                   "<br>","Blue line represents generalised additive model trend"
                   ),
    fill = "Num.<br>samples"
    )+
  # scale_colour_binned() +
  scale_fill_stepsn(
    # colours = c("red", "yellow", "green", "yellow", "red"),
    # breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)
    colours = c("#DEEBF7", "#9ECAE1", "#3182BD"),
    breaks = c(100, 200, 300)
  )+
  theme(
    palette.color.continuous = c("#DEEBF7", "#9ECAE1", "#3182BD"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2,size=12),
    axis.text = element_text(face=2),
    axis.text.x = element_text(size = 12),
    plot.caption = ggtext::element_markdown(face=2,size=12),
    plot.title = element_text(face=2,size = 14),
    legend.title = ggtext::element_markdown(face=2),
    )
rm(pad,ylims,df_sum_tmp)
dev.off()

# summarise carbon by sampling event & plot ####
# png(file = paste0(outfol,"phyto_carbon_raw.png"),
png(file = paste0(outfol,"phyto_carbon_log10.png"),
    width=18*ppi, height=12*ppi, res=ppi)
df %>% #names()
  ## total carbon: cells per litre
  select(sample_id,date_plot,md_carbon_tot_ug_c_per_m3) %>% 
  ##remove zero values
  dplyr::filter(!is.na(md_carbon_tot_ug_c_per_m3)) %>% 
  group_by(across(-md_carbon_tot_ug_c_per_m3)) %>% 
  # get sum per sample
  summarise(md_carbon_tot_ug_c_per_m3 = sum(md_carbon_tot_ug_c_per_m3,
                                               na.rm = TRUE),
            .groups = "drop") %>% ungroup() %>% 
  dplyr::filter(md_carbon_tot_ug_c_per_m3 !=0) %>% 
  # get mean sum by date
  ## drop sample ID
  select(-sample_id) %>% group_by(date_plot) %>%
  # filter(!is.na(md_carbon_tot_mg_C_per_litre))
  summarise(md_carbon_tot_ug_c_per_m3 = median(md_carbon_tot_ug_c_per_m3),
            n=n(),
            .groups = "drop") %>% ungroup() %>% 
  ggplot(.,
         aes(
           x = date_plot,
           # y = cells_per_litre
           # y = md_carbon_tot_ug_c_per_m3,
           y = log10(md_carbon_tot_ug_c_per_m3+1)
         )
  )+
  geom_rug(sides = "b")+
  geom_vline(xintercept = as.Date(c("2000-01-01","2005-01-01",
                                    "2010-01-01","2015-01-01",
                                    "2020-01-01","2025-01-01"
                                    )
                                  ),
             lty = 2, col="grey")+
  geom_vline(xintercept = as.Date(c("2001-01-01","2002-01-01","2003-01-01","2004-01-01",
                                    "2006-01-01","2007-01-01","2008-01-01","2009-01-01",
                                    "2011-01-01","2012-01-01","2013-01-01","2014-01-01",
                                    "2016-01-01","2017-01-01","2018-01-01","2019-01-01",
                                    "2021-01-01","2022-01-01","2023-01-01","2024-01-01"
                                    )
                                  ),
             lty = 3, col="grey")+
  geom_point(aes(fill = n),size=4,pch=21)+
  geom_smooth(method = "gam")+
  labs(
    title = "Log10 median phytoplanton carbon content by year_month",
    # title = "Median phytoplanton carbon content by year_month",
    y = "Log10(n+1) median carbon content (um/m3) per sample",
    # y = "Median carbon content (um/m3) per sample",
    caption=paste0("Values represent median total carbon contents across all samples gathered in a particular month","<br>",
                   "Point shading reflects number of samples contributing to that mean",
                   "<br>","Blue line represents generalised additive model trend"
                   ),
    fill = "Num.<br>samples"
    )+
  # scale_colour_binned() +
  scale_fill_stepsn(
    # colours = c("red", "yellow", "green", "yellow", "red"),
    # breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)
    colours = c("#FEE0D2", "#FC9272", "#DE2D26"),
    breaks = c(100, 200, 300)
  )+
  theme(
    # palette.color.continuous = c("#FEE0D2", "#FC9272", "#DE2D26"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2,size=12),
    axis.text = element_text(face=2),
    axis.text.x = element_text(size = 12),
    plot.caption = ggtext::element_markdown(face=2,size=12),
    plot.title = element_text(face=2,size = 14),
    legend.title = ggtext::element_markdown(face=2),
  )
dev.off()

# chlorophyll content ####
## load data
dfwims0 <- readRDS("//prodds.ntnl/Shared/AN/KFH/Groups/N_Marine/Temporary only - delete in 24 days/For Chris/DATA/Phyto_longTerm/phyto_wims.Rdat")

## turn into one dataframe

dfwims <- dfwims0$meta
dfwims$chla <- dfwims0$wims$`Chlorophyll_ug/l`

# generate plot
png(file = paste0(outfol,"phyto_chla_log10.png"),
# png(file = paste0(outfol,"phyto_chla_raw.png"),
    width=18*ppi, height=12*ppi, res=ppi)
dfwims %>% #names()
  select(date,month, year,chla) %>% 
  mutate(date_plot = make_date(year = year,
                               month = month,
                               day =1)) %>% 
  select(date_plot,chla) %>% 
  filter(!is.na(chla)) %>% 
  group_by(date_plot) %>% 
  summarise(chla = mean(chla, na.rm = TRUE),
            n=n(),
            .groups = "drop",
            ) %>% 
  ungroup() %>% 
  # filter(date_plot >"2006-12-31") %>% 
  ggplot(., aes(
    x = date_plot,
    # y = chla,
    y = log10(chla),
    )
    )+
  geom_rug(sides = "b")+
  geom_vline(xintercept = as.Date(c("2000-01-01","2005-01-01",
                                    "2010-01-01","2015-01-01",
                                    "2020-01-01","2025-01-01"
  )
  ),
  lty = 2, col="grey")+
  geom_vline(xintercept = as.Date(c("2001-01-01","2002-01-01","2003-01-01","2004-01-01",
                                    "2006-01-01","2007-01-01","2008-01-01","2009-01-01",
                                    "2011-01-01","2012-01-01","2013-01-01","2014-01-01",
                                    "2016-01-01","2017-01-01","2018-01-01","2019-01-01",
                                    "2021-01-01","2022-01-01","2023-01-01","2024-01-01"
  )
  ),
  lty = 3, col="grey")+
  geom_point(aes(fill = n),size=4,pch=21)+
  geom_smooth(method = "gam")+
  labs(
    # title = "Chlorophyll content by year_month",
    title = "Log10 chlorophyll content by year_month",
    # y = "Mean chlorophyll content (ug/l) per sample",
    y = "Log10 mean chlorophyll content (ug/l) per sample",
    caption=paste0("Values represent mean chlorophyll contents across all samples gathered in a particular month",
                   "<br>",
                   "Point shading reflects number of samples contributing to that mean",
                   "<br>","Blue line represents generalised additive model trend"
                   ),
    fill = "Num.<br>samples"
    )+
  # scale_colour_binned() +
  scale_fill_stepsn(
    # colours = c("red", "yellow", "green", "yellow", "red"),
    # breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)
    colours = c("#D7FAD8", "#7BBD7D", "#134F15"),
    breaks = c(100, 200, 300)
  )+
  theme(
    # palette.color.continuous = c("#D7FAD8", "#7BBD7D", "#134F15"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2),
    axis.text = element_text(face=2),
    axis.text.x = element_text(size = 12),
    plot.caption = ggtext::element_markdown(face=2,size=12),
    plot.title = element_text(face=2,size = 14),
    legend.title = ggtext::element_markdown(face=2),
  )
dev.off()
