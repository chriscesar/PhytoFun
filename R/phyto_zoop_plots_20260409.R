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

# initial plot: all data ####
df %>% 
  ggplot(.,
         aes(
           x = date_plot,
           y = log10(md_carb_ug_c_per_m3+1),
           colour = data_set,
           fill = NA,
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
           colour = data_set,
           fill = data_set,
           group = data_set,
           shape = data_set,
         )
  )+
  geom_point(alpha=0.1)+
  geom_smooth(method = "gam",
              formula = y ~ s(x, bs = "ds",k=10),# for date
              # formula = y ~ s(x, bs = "cc",k=10),# for day of year
              se=FALSE, aes(colour = data_set))+
  scale_shape_manual(values = c(21,24))+
  scale_fill_manual(values = c("Phytoplankton" = "green",
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
