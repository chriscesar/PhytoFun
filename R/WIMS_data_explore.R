# WIMS_data_explore.R ----
## AIMS ----
# import wims marine data
# for each unique site id, determine the first and last sample dates we have

ld_pkgs <- c("tidyverse","tictoc","stringr", "arrow")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

source("R/000setup.R")

# load data ####
tic("load data")
# dfsummary <- read.csv(file = "outputs/archive/Phyto_by_sample.csv")

dfwimsSites <- readxl::read_xlsx("data/all.saline.sites_WBs_EastNorth.xlsx",
                                 sheet = "SalineSitesJoined",
                                 guess_max = 10000
                                 )
toc(log=TRUE)

# Pull data from PARQUET files ####
tic("Pull data from PARQUET files")
# estuarine
wims.est <- open_dataset("C:/WIMS_cache/WIMS_est")
# 'O:/National_Hydroecology/WIMS_cache/WIMS_freshwater')

# sea
wims.sea <- open_dataset("C:/WIMS_cache/WIMS_sea")
toc(log = TRUE)

tic("Pull estuary data")
# est
est_phyto_site <- wims.est |>
  collect() |>
  dplyr::select(SAMP_SMPT_USER_REFERENCE, res) |>
  dplyr::distinct()
toc(log=TRUE)

tic("Pull marine data")
# sea
sea_phyto_site <- wims.sea |>
  collect() |>
  dplyr::select(SAMP_SMPT_USER_REFERENCE, res) |>
  dplyr::distinct()
toc(log=TRUE)
toc()

sites <- rbind(est_phyto_site, sea_phyto_site)

# read data for identified sites

site.list <- data.frame(WIMSCode = sites$SAMP_SMPT_USER_REFERENCE,
                        Region = sites$res)

site.list <- 
  site.list %>% 
  mutate(notation = paste(Region, WIMSCode, sep = '-'))

# extract data for selected sites ----
tic()
hms.data.est <- 
  wims.est |>
  filter(notation %in% site.list$notation) %>%
  # filter(time_period != 1) %>% 
  # filter(DATE_TIME > "2000-01-01 00:00:00") %>% 
  # & MEAS_DETERMINAND_CODE %in%  c('0076', '0180')) %>%
  collect()

hms.data.sea <- 
  wims.sea |>
  filter(notation %in% site.list$notation) %>%
  # filter(time_period != 1) %>% 
  # filter(DATE_TIME > "2000-01-01 00:00:00") %>% 
  # filter(time_period != 1 & MEAS_DETERMINAND_CODE %in%  c('0076', '0180')) %>%
  collect()

wimsdat <- bind_rows(hms.data.est, hms.data.sea) |>
  mutate(
    sample_date = as_date(
      parse_date_time(
        as.character(DATE_TIME),
        orders = c("ymd HMS", "ymd")
      )
    )
  )
rm(hms.data.est,hms.data.sea)
toc()

# calculate number of samples, first date & last date for monitoring
tic("summarise site data")
wimsdat %>% 
  dplyr::group_by(notation) %>% 
  dplyr::summarise(
    min_date = min(DATE_TIME),
    max_date = max(DATE_TIME),
    n = n()
  ) %>% 
  left_join(
    ., dfwimsSites %>% dplyr::select(
      notation,SMPT_SHORT_NAME,Easting,
      Northing,RIVER_BASI,WBID,WBName, Type,SMPT_TYPE,
      MODIFIED_W,COUNTRY,SMPT_STATUS
      ),
    by = "notation"
  ) -> site_ranges0

# remove mis-assigned (typically river) sites
site_ranges0 %>% 
  dplyr::filter(
    !is.na(SMPT_TYPE)
  ) -> site_ranges

saveRDS(site_ranges0,file="outputs/site_dates_unfilt.Rdat")
write.csv(site_ranges0,file="outputs/site_dates_unfilt.csv",row.names = FALSE)
saveRDS(site_ranges,file="outputs/site_dates.Rdat")
write.csv(site_ranges,file="outputs/site_dates.csv",row.names = FALSE)
toc()

site_ranges <- readRDS("outputs/site_dates.Rdat")

# Decade boundaries
decade_lines <- seq(
  as.Date("1980-01-01"),
  as.Date("2020-01-01"),
  by = "10 years"
  )

# 5-year boundaries (excluding decades to avoid double plotting)

fiveyear_lines <- seq(
  as.Date("1985-01-01"),
  as.Date("2025-01-01"),
  by = "10 years")


site_ranges %>% 
  dplyr::filter(!is.na(RIVER_BASI)) %>% 
  dplyr::mutate(
    range_days = as.numeric(
      difftime(max_date, min_date, units = "days")
      )) %>% 
  dplyr::relocate(range_days, .after = notation) %>% 
  arrange(desc(range_days)) %>% # longest date range first
  mutate(
    notation = factor(notation, levels = notation)
    ) %>%
  mutate(status = case_when(
    SMPT_STATUS == "C" ~ "Closed",
    SMPT_STATUS == "O" ~ "Open"
  )) %>% 
  mutate(status = factor(status, levels = c("Open","Closed"))) %>% 
  mutate(RBD = factor(RIVER_BASI, levels = c(
    "Northumbria",
    "Humber",
    "Anglian",
    "Thames",
    "South East",
    "South West",
    "Severn",
    "North West",
    "Solway Tweed"
  ))) %>% 
  ggplot(.,
       aes(y = notation)) +
  geom_segment(
    aes(x = min_date, xend = max_date, yend = notation,
        col = status,
        ),
    alpha = 0.2,
    linewidth = 1
  ) +
  geom_vline(xintercept = decade_lines,
             colour = "grey40",
             linewidth = 0.6,
             linetype = 2
  ) +
  geom_vline(xintercept = fiveyear_lines,
             colour = "grey70",
             linewidth = 0.4,
             linetype = 3
  ) +
  # geom_point(aes(x = min_date,
  #                fill = status
  #                ), size = 3,pch=21,
  #            alpha = 0.4
  #            ) +
  geom_point(aes(x = max_date,
                 fill = status
                 ), size = 3,pch=21,
             alpha = 0.4
             ) +
  scale_fill_manual(values = c("Open" = "azure2","Closed" = "coral1"))+
  scale_colour_manual(values = c("Open" = "deepskyblue2","Closed" = "coral1"))+
  labs(
    x = "Date",
    y = "Site"
  ) +
  ggthemes::theme_few()+
  facet_wrap(.~RBD)+
  labs(
    title = "Timescale over which individual EA WIMS sampling sites in estuarine and coastal waters have been monitored",
    subtitle = "Points indicate most recent sampling event, with tails indicating timescale over which monitoring has occurred.
Shading indicates whether a sample is currently flagged as 'Open' or 'Closed'",
    caption = "Each line/point item represents a separate sampling station as defined by WIMS SAMP_SMPT_USER_REFERENCE values"
  ) +
  theme(
    strip.text = element_text(face = 2, size=14),
    plot.title = element_text(face = 2),
    plot.subtitle = element_text(face = 2),
    plot.caption = element_text(face = 2, size = 11, lineheight = 0.5),
    axis.title.x = element_text(face = 2, size=12),
    axis.title.y = element_blank(),
    axis.text.x = element_text(face = 2),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(face = 2, size = 11),
  )

# create BIG data for each sampling event
wimsdat %>% 
  left_join(
    ., dfwimsSites %>% dplyr::select(
      notation,SMPT_SHORT_NAME,Easting,
      Northing,RIVER_BASI,WBID,WBName, Type,SMPT_TYPE,
      MODIFIED_W,COUNTRY
    ),
    by = "notation"
  ) %>% #names()
  dplyr::select("notation","SAMP_SMPT_USER_REFERENCE","SAMP_ID",
                "DATE_TIME","SAMP_MATERIAL","SAMP_PURPOSE_CODE",
                "res", "wims_region","time_period","material_type",
                "sample_date","SMPT_SHORT_NAME","Easting", "Northing",
                "RIVER_BASI","WBID", "WBName","Type","SMPT_TYPE",
                "MODIFIED_W","COUNTRY") %>% 
  dplyr::distinct() %>% 
  dplyr::filter(
    !is.na(SMPT_TYPE)
    ) -> smps_all

smps_all %>% 
  ggplot(.,
         aes(
           x = DATE_TIME,
           y = notation
         ))+
  geom_point(alpha = 0.2)+
  # facet_grid(.~RIVER_BASI, rows = 1)
  facet_wrap(.~ RIVER_BASI, nrow = 1, scales = "free_y")+
  theme(
    axis.text.y = element_blank(),
  )

library(dplyr)
library(lubridate)
library(ggplot2)

# Count unique samples per month and region
monthly_samples <- smps_all %>%
  dplyr::mutate(
    year_month = lubridate::floor_date(sample_date, unit = "month")
  ) %>%
  dplyr::group_by(wims_region, year_month) %>%
  dplyr::summarise(
    n_samples = n_distinct(notation),
    .groups = "drop"
  ) %>% 
  dplyr::mutate(region_lb = case_when(
    wims_region == "wimsnepr" ~ "North East",
    wims_region == "wimsanpr" ~ "Anglian",
    wims_region == "wimsthpr" ~ "Thames",
    wims_region == "wimssopr" ~ "Southern",
    wims_region == "wimsswpr" ~ "South West",
    wims_region == "wimsmipr" ~ "Midlands",
    wims_region == "wimsnwpr" ~ "North West",
  )) %>% 
  dplyr::mutate(region_lb = factor(
    region_lb,
    levels = c(
      "North East",
      "Anglian",
      "Thames",
      "Southern",
      "South West",
      "Midlands",
      "North West"
    )
  )) %>% 
  dplyr::mutate(region = factor(wims_region,
                                levels = c(
                                  "wimsnepr",
                                  "wimsanpr",
                                  "wimsthpr",
                                  "wimssopr",
                                  "wimsswpr",
                                  "wimsmipr",
                                  "wimsnwpr"
                                  )
                                ))
  
# Visualise temporal trends ----
# Decade boundaries
decade_lines <- seq(
  as.Date("1980-01-01"),
  as.Date("2020-01-01"),
  by = "10 years"
  )

# 5-year boundaries (excluding decades to avoid double plotting)

fiveyear_lines <- seq(
  as.Date("1985-01-01"),
  as.Date("2015-01-01"),
  by = "10 years")

ggplot(monthly_samples,
       aes(x = year_month, y = n_samples)) +
  geom_vline(xintercept = decade_lines,
             colour = "grey40",
             linewidth = 0.6,
             linetype = 2
             ) +
  geom_vline(xintercept = fiveyear_lines,
             colour = "grey70",
             linewidth = 0.4,
             linetype = 3
             ) +
  geom_point(colour = "steelblue") +
  coord_cartesian(ylim = c(0, NA)) +
  facet_wrap(~ region_lb,
             scales = "free_y"
             ) +
  labs(
    x = "Year-Month",
    y = "Number of unique samples",
    title = "Monthly sampling effort by WIMS region",
    caption = "Points indicate the number of unique sampling station ID values in WIMS data within a given year-month combination.\n
    Note that ordinates (y-axes) are not universally scaled.\n
    Apparent dichotomy in values (e.g., South West) due to reduced successful sampling events in winter months compared to other seasons."
  ) +
  ggthemes::theme_few() +
  geom_smooth(
    method = "gam",
    colour = 2
  ) +
  theme(
    strip.text = element_text(face = 2, size=14),
    plot.title = element_text(face = 2),
    plot.caption = element_text(face = 2, size = 10, lineheight = 0.5),
    axis.title = element_text(face = 2, size=12),
    axis.text = element_text(face = 2),
  )
