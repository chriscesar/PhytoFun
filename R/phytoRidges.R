# phytoRidges.R ####
## Generate ridgeplot

# load packages ####
ld_pkgs <- c("tidyverse","ggplot2","sf","maps",
             "ggpubr", "ggspatial","ggrepel","tictoc",
             "ggridges", "vegan", "lubridate")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)

rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log
theme_set(ggthemes::theme_few())

df0 <- readRDS("data/Phyto_WB_wide.rdat")
df0$Region <- ifelse(df0$AGENCY_AREA == "ANGLIAN - CENTRAL","Anglian",
                     ifelse(df0$AGENCY_AREA == "ANGLIAN - EASTERN","Anglian",
                            ifelse(df0$AGENCY_AREA == "ANGLIAN - NORTHERN","Anglian",
                                   ifelse(df0$AGENCY_AREA == "NORTH EAST - NORTH EAST","NORTH EAST",
                                          ifelse(df0$AGENCY_AREA == "NORTH EAST - YORKSHIRE","NORTH EAST",
                                                 ifelse(df0$AGENCY_AREA == "NORTH WEST - NORTH","NORTH WEST",
                                                        ifelse(df0$AGENCY_AREA == "NORTH WEST - SOUTH","NORTH WEST",
                                                               ifelse(df0$AGENCY_AREA == "SOUTH WEST - CORNWALL","SOUTH WEST",
                                                                      ifelse(df0$AGENCY_AREA == "SOUTH WEST - DEVON","SOUTH WEST",
                                                                             ifelse(df0$AGENCY_AREA == "SOUTH WEST - NORTH WESSEX","SOUTH WEST",
                                                                                    ifelse(df0$AGENCY_AREA == "SOUTH WEST - SOUTH WESSEX","SOUTH WEST",
                                                                                           ifelse(df0$AGENCY_AREA == "SOUTHERN - KENT & E. SUSSEX","SOUTHERN",
                                                                                                  ifelse(df0$AGENCY_AREA == "SOUTHERN - SOLENT & S. DOWNS","SOUTHERN","ERROR")
                                                                                           ))))))))))))                                          
                                                                      
                            )

names(df0)

df0$S <- df0 %>% 
  select(starts_with("NBNSYS")) %>%
  vegan::specnumber(.)

ggplot(df0, aes(x=S,y=as.factor(Year),height = S))+
  geom_point()+
  facet_wrap(.~AGENCY_AREA)

df0 %>% 
  filter(.,AGENCY_AREA !="NATIONAL MARINE") %>% 
  ggplot(., aes(x=S,y=as.factor(Year)))+
  geom_density_ridges(#stat="identity",
    scale=1.4,alpha=0.7,
    show.legend = FALSE,
    aes(fill=Region),
    colour="darkgrey"
      )+
  facet_wrap(.~AGENCY_AREA)+
  scale_y_discrete(limits=rev)



