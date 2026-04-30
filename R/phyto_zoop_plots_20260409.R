# phyto_zoop_plots_20260409.R ####

# load packages & data ####
ld_pkgs <- c("tidyverse","ggthemes","ggh4x")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log
source("R/000setup.R")

outfol <- "//prodds.ntnl/Shared/AN/KFH/Groups/N_Marine/02 Projects_Tasks/04 E&B/Ecology and Ecosystems/Planktonic_Food_Webs/202602_FigsForReport/"

df0 <- readRDS(file = "outputs/ZoopPhytoMatchingBIOSYS_long.Rdat")

# generate factor levels
fac_tmp <- sprintf("%04d%02d", rep(2000:2025, each = 12), 1:12)

## keep only 'interesting' variables
df <- df0 %>% 
  select(sample_id,rbd,wb_name,sample_date,data_set,taxon,lifeform,
         abundance:md_carb_ug_c_per_m3
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

# FOR COMPARISON WITH OTHER METHOD ####
df %>%
  ## total abundance: cells per litre
  select(sample_id,rbd,wb_name,sample_date,
         date_plot,md_carb_ug_c_per_m3) %>% #View()
  group_by(across(-md_carb_ug_c_per_m3)) %>% 
  # get sum per sample
  summarise(md_carb_ug_c_per_m3 = sum(md_carb_ug_c_per_m3,
                                            na.rm = TRUE),
            .groups = "drop") %>% 
  # get mean sum by date
  ## drop sample ID
  select(-sample_id) %>% group_by(date_plot) %>% 
  summarise(mn_md_carb_ug_c_per_m3 = mean(md_carb_ug_c_per_m3,na.rm = TRUE),
            sd_md_carb_ug_c_per_m3 = sd(md_carb_ug_c_per_m3,na.rm=TRUE),
            n = n(),
            .groups = "drop"
  ) %>%
  write.csv(.,
            file="data/TEST_phyto_from_phyto_zoop_plots_20260409.R.csv",
            row.names = FALSE)

# initial plot: all data ####
df %>% 
  ggplot(.,
         aes(
           x = date_plot,
           y = log10(md_carb_ug_c_per_m3+1),
           colour = data_set,
           # fill = NA,
           group = data_set,
         )
  )+
  geom_point(alpha=0.01)+
  geom_smooth(se=FALSE)+
  # scale_fill_manual(values = c("Phytoplankton" = "green",
  #                              "Zooplankton" = "blue"))+
  scale_colour_manual(values = c("Phytoplankton" = "darkgreen",
                                 "Zooplankton" = "darkblue"))+
  facet_wrap(.~rbd)

# Sum by sample: all data ####
df %>% #names() %>% 
  dplyr::select(sample_id, rbd,wb_name, sample_date, date_plot,data_set,
                md_carb_ug_c_per_m3) %>% 
  group_by(across(-md_carb_ug_c_per_m3)) %>% 
  summarise(md_carb_ug_c_per_m3 = sum(md_carb_ug_c_per_m3, na.rm = TRUE),
            .groups = "drop") %>% ungroup() %>% 
  filter(md_carb_ug_c_per_m3>0) %>% 
  ggplot(.,
         aes(
           x = date_plot,
           y = log10(md_carb_ug_c_per_m3+1),
           colour = data_set,
           fill = data_set,
           group = data_set,
         )
  )+
  geom_point(alpha=0.1)+
  geom_smooth(se=FALSE)+
  scale_fill_manual(values = c("Phytoplankton" = "green",
                               "Zooplankton" = "blue"))+
  scale_colour_manual(values = c("Phytoplankton" = "darkgreen",
                                 "Zooplankton" = "darkblue"))+
  facet_wrap(.~rbd)

# Sum by sample: median by year_month all years ####
df %>% #names() %>% 
  dplyr::select(sample_id, rbd,wb_name, date_plot,data_set,
                md_carb_ug_c_per_m3) %>% 
  group_by(across(-md_carb_ug_c_per_m3)) %>% 
  summarise(md_carb_ug_c_per_m3 = sum(md_carb_ug_c_per_m3, na.rm = TRUE),
            .groups = "drop") %>% ungroup() %>% 
  filter(md_carb_ug_c_per_m3>0) %>%
  group_by(across(-md_carb_ug_c_per_m3)) %>% 
  summarise(md_carb_ug_c_per_m3 = median(md_carb_ug_c_per_m3,na.rm=TRUE),
            .groups = "drop") %>% ungroup() %>% 
  ggplot(.,
         aes(
           x = date_plot,
           y = log10(md_carb_ug_c_per_m3+1),
           fill = data_set,
           group = data_set,
           shape = data_set,
         )
  )+
  geom_point(alpha=0.1,
             aes(
               colour = data_set,
             ))+
  geom_smooth(method = "gam",
              formula = y ~ s(x, bs = "ds",k=10),# for date
              # formula = y ~ s(x, bs = "cc",k=10),# for day of year
              se=FALSE,
              show.legend = FALSE,
              # aes(colour = data_set)
              colour = 2,
              )+
  scale_shape_manual(values = c(21,24))+
  scale_fill_manual(values = c("Phytoplankton" = "darkgreen",
                               "Zooplankton" = "blue"))+
  scale_colour_manual(values = c("Phytoplankton" = "green",
                                 "Zooplankton" = "blue"))+
  facet_wrap(.~rbd)+
  labs(
    y = "Log10(n+1) Total carbon content (as micrograms carbon per m3)"
      )+
  theme(
    strip.text = element_text(face=2),
    legend.text = element_text(face=2),
    legend.title = element_blank(),
    axis.title = element_text(face=2),
    axis.title.x = element_blank(),
    axis.text = element_text(face=2),
  )

# Sum by sample: median by year_month zoop WBs only ####
zp_wbs <- df %>% filter(data_set=="Zooplankton") %>%
  pull(wb_name) %>%
  unique()

mn_zoop <- df %>% filter(data_set == "Zooplankton") %>% 
  pull(sample_date) %>% min()

df %>% #names() %>% 
  filter(wb_name %in% zp_wbs) %>% 
  filter(sample_date >mn_zoop-1) %>% 
  dplyr::select(sample_id, rbd,wb_name, date_plot,data_set,
                md_carb_ug_c_per_m3) %>% 
  group_by(across(-md_carb_ug_c_per_m3)) %>% 
  summarise(md_carb_ug_c_per_m3 = sum(md_carb_ug_c_per_m3, na.rm = TRUE),
            .groups = "drop") %>% ungroup() %>% 
  filter(md_carb_ug_c_per_m3>0) %>%
  group_by(across(-md_carb_ug_c_per_m3)) %>% 
  summarise(md_carb_ug_c_per_m3 = median(md_carb_ug_c_per_m3,na.rm=TRUE),
            .groups = "drop") %>% ungroup() %>% 
  # organise order of WBs clockwise from north east
  mutate(wb_name = factor(wb_name, levels = c(
    "Northumberland North","Farne Islands to Newton Haven", "TEES",
    "Yorkshire South",
    "Lincolnshire","Wash Outer","Blackwater Outer",
    "THAMES LOWER",
    "Kent South","Isle of Wight East","Solent","SOUTHAMPTON WATER",
    "Cornwall North", "Barnstaple Bay","Bristol Channel Inner South",
    "Mersey Mouth","Solway Outer South"))) %>% 
  ggplot(.,
         aes(
           x = date_plot,
           y = log10(md_carb_ug_c_per_m3+1),
           # colour = data_set,
           fill = data_set,
           group = data_set,
           shape = data_set,
         )
  )+
  geom_point(
    # alpha=0.1
    )+
  ylim(0,NA)+
  geom_smooth(method = "gam",
              formula = y ~ s(x, bs = "ds",k=10),# for date
              # formula = y ~ s(x, bs = "cc",k=10),# for day of year
              se=FALSE, aes(colour = data_set))+
  scale_shape_manual(values = c(21,24))+
  scale_fill_manual(values = c("Phytoplankton" = "green",
                               "Zooplankton" = "blue"))+
  scale_colour_manual(values = c("Phytoplankton" = "green",
                                 "Zooplankton" = "blue"))+
  # facet_wrap(.~wb_name)+
  facet_wrap2(
    vars(rbd,wb_name),
    strip = strip_nested(bleed = TRUE)
  )+
  # facet_nested(~ ,cols = 5)+
  labs(
    y = "Log10(n+1) Total carbon content (as micrograms carbon per m3)"
  )+
  theme(
    strip.text = element_text(face=2),
    legend.text = element_text(face=2),
    legend.title = element_blank(),
    axis.title = element_text(face=2),
    axis.title.x = element_blank(),
    axis.text = element_text(face=2),
    strip.background = element_rect(color = 1)
  )

# summarise carbon by sampling event & plot ####
## which data? ####
### Choose Phyto ####
nom <- "Phytoplankton"
line_dates_bold <- as.Date(c("2000-01-01","2005-01-01",
                             "2010-01-01","2015-01-01",
                             "2020-01-01","2025-01-01"
                             ))
line_dates_ann <- as.Date(c("2001-01-01","2002-01-01","2003-01-01",
                            "2004-01-01","2006-01-01","2007-01-01",
                            "2008-01-01","2009-01-01","2011-01-01",
                            "2012-01-01","2013-01-01","2014-01-01",
                            "2016-01-01","2017-01-01","2018-01-01",
                            "2019-01-01","2021-01-01","2022-01-01",
                            "2023-01-01","2024-01-01")
                          )

### Choose Zoops ####
nom <- "Zooplankton"
line_dates_bold <- as.Date(c("2023-01-01","2024-01-01","2025-01-01")
                           )
line_dates_ann <- as.Date(c(
  "2022-06-01"
)
                          )
### Run ####
df %>% 
  as_tibble() %>% 
  filter(data_set==nom) %>% #names() %>%
  ## total carbon: cells per litre
  dplyr::select(sample_id, rbd,wb_name, date_plot,data_set,
                md_carb_ug_c_per_m3) %>% 
  # calculate total C per sampling event
  group_by(across(-md_carb_ug_c_per_m3)) %>% 
  summarise(md_carb_ug_c_per_m3 = sum(md_carb_ug_c_per_m3,na.rm = TRUE),
            .groups = "drop") %>% ungroup() %>% 
  dplyr::select(date_plot,data_set,md_carb_ug_c_per_m3) %>% 
  group_by(across(-md_carb_ug_c_per_m3)) %>% 
  summarise(md_carb_ug_c_per_m3 = median(md_carb_ug_c_per_m3),
            n=n(),
            .groups = "drop") %>% ungroup() %>%
  ggplot(.,
         aes(
           x = date_plot,
           y = log10(md_carb_ug_c_per_m3+1),
           # group = data_set,
           # colour = data_set,
         )
  )+
  geom_rug(sides = "b")+
  geom_vline(xintercept =line_dates_bold,lty = 2, col="grey")+
  geom_vline(xintercept = line_dates_ann,lty = 3, col="grey")+
  geom_point(aes(fill = n),size=4,pch=21)+
  # geom_hline(yintercept = 4,col=2,lwd=1.5)+
  geom_smooth(method = "gam")+
  # facet_wrap(.~rbd)+
  # ylim(0,NA)+
  labs(
    # title = "Log10 median phytoplanton carbon content by year_month",
    title = paste0("Median ",tolower(nom)," carbon content by year_month"),
    # y = "Log10(n+1) median carbon content (ug/m3) per sample",
    y = "Median carbon content (ug/m3) per sample",
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
    # breaks = c(100, 200, 300)
  )+
  theme(
    palette.color.continuous = c("#FEE0D2", "#FC9272", "#DE2D26"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2,size=12),
    axis.text = element_text(face=2),
    strip.text = element_text(face=2,size=12),
    axis.text.x = element_text(size = 12),
    plot.caption = ggtext::element_markdown(face=2,size=12),
    plot.title = element_text(face=2,size = 14),
    legend.title = ggtext::element_markdown(face=2),
  )
ggsave(plot = get_last_plot(),
       filename = paste0(outfol,nom,"_carbon_log10.png"),
       device = "png",
       width=14, height=8,
       # width=18, height=12,
       # res=ppi
       )
