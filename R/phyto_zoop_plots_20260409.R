# phyto_zoop_plots_20260409.R ####

# load packages & data ####
ld_pkgs <- c("tidyverse","ggthemes")
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

df %>% 
  ggplot(.,
         aes(
           x = date_plot,
           y = log10(md_carb_ug_c_per_m3+1),
           colour = data_set,
           group = data_set,
         )
         )+
  geom_point(alpha=0.1)+
  geom_smooth()+
  facet_wrap(.~rbd)
