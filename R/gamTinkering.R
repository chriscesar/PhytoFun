#gamTinkering.R
ppi <- 300
library(mgcv);library(gratia)
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



df0$date <- as.Date(df0$SAMPLE_DATE)

df0$month <- month(df0$date)

df0 %>% 
  dplyr::select(.,c(1:12),NHMSYS0020749164_,Year:month) -> df1


table(df1$WATER_BODY)

df1 %>% filter(.,WATER_BODY=="TEES")


df1$date_num <- scale(as.numeric(df1$date))

fit1 <- mgcv::gam(
  NHMSYS0020749164_ ~ 
    s(date_num, bs = "gp") + #longer term process
    s(month, bs = "cc"), # account for seasonality
  data = df1,
  family = gaussian(),# poisson() fails as values expressed as numbers per ml(?)
  method = "REML"
)
draw(fit1) -> pl1
saveRDS(fit1,file="outputs/TEES_Time_NHMSYS0020749164.Rdat")
png(filename = "outputs/TEES_Time_NHMSYS0020749164.png",width=12*ppi, height=6*ppi, res=ppi)
print(pl1)
dev.off()
