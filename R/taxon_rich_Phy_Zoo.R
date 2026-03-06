## taxon_rich_Phy_Zoo.R ####

# load packages ####
ld_pkgs <- c("tidyverse","tictoc","ggthemes","ggh4x","vegan")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

source("R/000setup.R")
outfol <- "//prodds.ntnl/Shared/AN/KFH/Groups/N_Marine/02 Projects_Tasks/04 E&B/Ecology and Ecosystems/Planktonic_Food_Webs/202602_FigsForReport/"

# load data ####
df0 <- readRDS("outputs/ZoopPhytoMatchingBIOSYS.Rdat")
df0$region <- factor(df0$region, levels =c(
  "Northumbria",
  "Humber",
  "Anglian",
  "Thames",
  "South East",
  "South West",
  "North West"
  ))

df0$wb <- factor(df0$wb, levels = c(
  "Northumberland North","Farne Islands to Newton Haven","TEES",
  "Yorkshire South",
  "Lincolnshire","Lincs Offshore","Wash Outer","Blackwater Outer",
  "THAMES LOWER",
  "Kent South","Isle of Wight East","Solent","SOUTHAMPTON WATER",
  "Barnstaple Bay","Bristol Channel Inner South","Cornwall North",
  "Mersey Mouth","Solway Outer South"
  ))

df <- list()
df$zoop <- df0 %>% filter(data_set=="Zooplankton")
df$phyto <- df0 %>% filter(data_set=="Phytoplankton")

df$zoop %>% as_tibble() %>% 
  select(region,wb_id,wb,wb_type,
         biosys_short,sample_date,taxon,abundance) %>% 
  group_by(across(-abundance)) %>% 
  summarise(abundance = sum(abundance),.groups = "drop") %>%
  ungroup() %>% 
  filter(!is.na(abundance)) %>% 
  filter(abundance != 0) %>% 
  pivot_wider(names_from = taxon,
              values_from = abundance,
              values_fill = 0)-> df$zoop_w

df$phyto %>% as_tibble() %>%
  select(region,wb_id,wb,wb_type,
         biosys_short,sample_date,taxon,abundance) %>% 
  group_by(across(-abundance)) %>% 
  summarise(abundance = sum(abundance),.groups = "drop") %>%
  ungroup() %>% 
  filter(!is.na(abundance)) %>% 
  filter(abundance != 0) %>% 
  filter(sample_date>min(df$zoop_w$sample_date)-1) %>% 
  pivot_wider(names_from = taxon,
              values_from = abundance,
              values_fill = 0)-> df$phyto_w


df$zoop_w %>%
  select(-c(region,wb_id,wb,
            wb_type,biosys_short,
            sample_date)) %>% 
  mutate(S = vegan::specnumber(.)) %>% 
  pull(S)->zp_S

df$phyto_w %>%
  select(-c(region,wb_id,wb,
            wb_type,biosys_short,
            sample_date)) %>% 
  mutate(S = vegan::specnumber(.)) %>% 
  pull(S)->ph_S

df$zoop_w$S <- zp_S
df$phyto_w$S <- ph_S

df$zoop_w %>% 
  ggplot(.,aes(
    x = sample_date,
    y = S,
    ),
  )+
  geom_point(aes(colour = region),show.legend = FALSE)+
  facet_wrap(.~wb)+
  theme(
    strip.text = element_text(face=2),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2),
    axis.text = element_text(face=2),
    axis.text.x = element_text(face=2, size = 14),
    )


df$phyto_w %>% 
  select(region,wb,sample_date,S) %>% 
  mutate(type = "Phytoplankton")->df_p

df$zoop_w %>% 
  select(region,wb,sample_date,S) %>% 
  mutate(type = "Zooplankton")->df_z

rbind(df_p,df_z) -> df_pl

# Build a strip background list that matches your facet order
zones <- unique(df_pl$region)
bg_cols <- cbPaletteFill[c(1,2,3,4,6,7,9)]
bg_cols <- setNames(bg_cols[seq_along(zones)], zones)

strip_elems <- lapply(zones, function(z)
  element_rect(fill = bg_cols[[z]], color = "black", linewidth = 1)
)

df_pl %>% 
  ggplot(.,aes(
    x = sample_date,
    y = S,
  ),
  )+
  geom_point(
    aes(
      colour = type,
      shape = type
      ),
    size=3,
    show.legend = TRUE)+
  scale_shape_manual(values = c(21,24))+
  # facet_wrap(.~wb,
  #            ) +
  facet_wrap2(. ~ wb, 
              strip = strip_themed(background_x = strip_elems)) +
  geom_smooth(method = "gam",
              aes(group=type,colour=type),
              se=FALSE
              )+
  labs(
    title = "Taxon richness in phytoplankton and zooplankton assemblages"
  )+
  theme(
    plot.title = element_text(size=16,face=2),
    strip.text = element_text(face=2),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face=2),
    axis.text = element_text(face=2),
    axis.text.x = element_text(face=2, size = 14),
    legend.title = element_blank(),
  )
########

library(ggh4x)
library(dplyr)

# named palette for region colours
bg_cols <- cbPaletteFill[c(1,2,3,4,6,7,9)]
bg_cols <- setNames(bg_cols, levels(df_pl$region))

# (1) facet order
facet_order <- levels(df_pl$wb)

# (2) map each facet to its region
facet_regions <- sapply(facet_order, function(w) {
  reg <- df_pl$region[df_pl$wb == w][1]
  as.character(reg)
})

# (3) create strip elements in facet order
strip_elems <- lapply(facet_regions, function(reg) {
  element_rect(
    fill = alpha(bg_cols[[reg]],0.2),
    colour = "black",
    linewidth = 1
  )
})

png(file = paste0(outfol,"phyto_zoop_richness.png"),
    width=18*ppi, height=12*ppi, res=ppi)
# (4) plot
df_pl %>%
  ggplot(aes(
    x = sample_date,
    y = S
  )) +
  geom_point(
    aes(
      colour = type,
      shape = type
    ),
    size = 3
  ) +
  scale_shape_manual(values = c(21, 24)) +
  facet_wrap2(
    . ~ wb,
    strip = strip_themed(background_x = strip_elems)
  ) +
  geom_smooth(
    method = "gam",formula = y ~ s(x, k = 7),
    aes(group = type, colour = type),
    se = FALSE
  ) +
  labs(
    title = "Taxon richness in phytoplankton and zooplankton assemblages",
    caption = "Lines represent GAM predictions. Strip colours reflect EA regions") +
  theme(
    legend.position ="inside",
    legend.position.inside = c(0.8,0.1),
    plot.title = element_text(size = 16, face = 2),
    plot.caption = element_text(size = 12, face = 2),
    strip.text = element_text(face = 2,size=14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = 2),
    axis.text = element_text(face = 2),
    axis.text.x = element_text(face = 2, size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size=12,face=2),
    )
dev.off()
