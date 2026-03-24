# expl_NP.R ####

# exploration of N to P ratios in WIMS data ####
# load packages ####
ld_pkgs <- c("tidyverse","tictoc","ggthemes","ggh4x")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

tictoc::tic("load data")
# load data ####

wimsdat_wb0 <- readRDS(file="outputs/wims_all.RDat")
tictoc::toc(log=TRUE)

# remove reactive, waste, external & unspecified monitoring
rmv <- c("UF","UI","WA","WF","WI","WO","XC","XO","ZZ")

wimsdat_wb <- wimsdat_wb0 %>% filter(!(SAMP_PURPOSE_CODE %in% rmv))
rm(rmv)

dets <- c("Phosphate","OrthophsFilt","Orthophospht",
          "Nitrate Filt","Nitrate-N",
          "Nitrite-N","Nitrite Filt","Ammonia(N)")

dets <- c(
  # Ammonia dets
  "9081","7374","8958","3603","3406","6599","5282","3163","3656","9803",
  "9150","6496","119","6382","1077","1167","7677","7004","7404","9259",
  "9260","3970","4027","7532","3561","6003","8508","110","9720","7403",
  "9397","6021","5092","111","9875","5506","9057","3403","9993","7186",
  "8078","8729","6872","8153",
  # Nitrate dets
  "5507","7194","8601","6004","8501","9232","6490","7633","6022","9395",
  "117","6380","9880","9853","9055","5567","9607","8070","5280","9881",
  # Nitrite dets
  "7634", "6005","8502", "9886","7635", "118", "9396","6023", "6381","3402",
  "9879","5508", "9056","3404", "6485","7891", "8079","5279",
  # Orthophosphate dets
  "188", "6008", "9876", "9398", "191", "6393", "9856", "9058", "189", 
  "8755", "8068", "8536", "6392", "180",
  # Phosphate dets
  "7315", "1174", "8504", "4127", "192", "6025", "9449", "9752", "190",
  "3507", "6394", "4190", "3845" )

wimsdat_wb %>% 
  filter(MEAS_DETERMINAND_CODE %in% dets) -> df_PN
  # filter(DETE_SHORT_DESC %in% dets) -> df_PN
rm(dets)



df_PN %>% 
  select(DETE_SHORT_DESC,UNIT_SHORT_DESC) %>% 
  distinct() 

# most compounds as mg/l, some as ug/l. Convert all to mg/l
df_PN %>% 
  mutate(val_mgl = case_when(
    UNIT_SHORT_DESC == "mg/l" ~ MEAS_RESULT,
    UNIT_SHORT_DESC == "ug/l" ~ MEAS_RESULT/1000,
    TRUE ~ NA
  )) %>% 
  relocate(val_mgl, .after = MEAS_RESULT) %>% 
  select(-MEAS_RESULT) -> df_PN

# widen data to calculate N:P ratio
df_PN %>% 
  select(
    RIVER_BASI,WBID,WBName,Type,EA_AREA_NA,
    material_type,sample_date,SAMP_ID,
    DETE_SHORT_DESC,MEAS_SIGN,val_mgl) %>% 
  mutate(val_mgl_sign = case_when(
    is.na(MEAS_SIGN) ~ as.character(val_mgl),
    MEAS_SIGN == "<" ~ paste0("<",as.character(val_mgl)),
    MEAS_SIGN == ">" ~ paste0(">",as.character(val_mgl)),
    MEAS_SIGN == " " ~ as.character(val_mgl),
    TRUE ~ NA)) %>% 
  select(-c(MEAS_SIGN,val_mgl))
    
table(df_PN$DETE_SHORT_DESC)
df_PN %>% 
  select(
    RIVER_BASI,WBID,WBName,Type,EA_AREA_NA,
    material_type,sample_date,SAMP_ID,
    DETE_SHORT_DESC,MEAS_SIGN,val_mgl) %>% 
  # ignore signs
  select(-MEAS_SIGN) %>% 
  group_by(across(-val_mgl)) %>% 
  summarise(val_mgl = mean(val_mgl),.groups = "drop") %>% ungroup() %>% 
  
  # 
  # # remove annoying dets
  # filter(DETE_SHORT_DESC != "Nitrate-N") %>% 
  # filter(DETE_SHORT_DESC != "Nitrite-N") %>% 
  

  pivot_wider(names_from = DETE_SHORT_DESC,
              values_from = val_mgl) -> df_PN_w

View(df_PN_w)  

df_PN_w %>% 
  ggplot(.,
         aes(
           x = sample_date,
           # y = OrthophsFilt,
           y=df_PN_w$`Nitrate Filt`,
         ),
         )+
  geom_point()+
  ggthemes::theme_few()+
  geom_smooth(method = "gam")+
  facet_wrap(.~RIVER_BASI,scales = "free_y")+
  theme(
    strip.text = element_text(face=2),
  )

df_PN %>% 
  ggplot(.,
         aes(
           x = sample_date,
           y = log10(val_mgl),
           colour = DETE_SHORT_DESC
         ))+
  ggthemes::theme_few()+
  geom_point(
    aes(
      group=DETE_SHORT_DESC
    ))+
  facet_wrap(.~RIVER_BASI,scales = "free_y")+
  theme(
    strip.text = element_text(face=2),
  )


## extract samples which contain values for BOTH OrthophsFilt & Nitrate Filt
df_PN_w %>% 
  filter(!is.na(OrthophsFilt), !is.na(`Nitrate Filt`)) %>% 
  mutate(N_to_P = `Nitrate Filt`/OrthophsFilt) %>% 
  ggplot(.,
         aes(
           x=sample_date,
           y = N_to_P,
         )
         )+
  ggthemes::theme_few()+
  geom_hline(yintercept = 16, lty=2,col=2)+
  geom_point()+
  geom_smooth(method = "gam")+
  facet_wrap(.~RIVER_BASI,scales = "free_y")+
  theme(
    strip.text = element_text(face=2),
  )


# mol values ####
df_PN_w %>% 
  mutate(
    N = rowSums(across(c(`NH3 filt N`, `Nitrate Filt`, `Nitrite Filt`)), na.rm = TRUE),
    N_mol_l = N/14010,
    P_mol_l = OrthophsFilt/30970,
    N_to_P = N_mol_l/P_mol_l,
    yday = lubridate::yday(sample_date),
  ) %>% 
  filter(sample_date>"2006-12-31") %>%
  ## remove likely typo Nitrate values
  filter(`Nitrate Filt` < 10) %>% 
  ggplot(.,
         aes(
           x=sample_date,
           # x=yday,
           y = N_to_P,
         )
  )+
  ggthemes::theme_few()+
  geom_hline(yintercept = 16, lty=2,col=2)+
  geom_hline(yintercept = 32, lty=2,col=2)+
  geom_point(alpha = 0.05,
             show.legend = FALSE,
             aes(
               colour = Type,
             ))+
  geom_smooth(method = "gam")+
  # facet_wrap(RIVER_BASI~Type,
  #            # scales = "free_y"
  #            )+
  facet_wrap2(
    vars(RIVER_BASI,Type), ncol = 4,
    strip = strip_nested(bleed = FALSE),
    scales = "free_y"
  )+
  labs(
    title = "Nitrogen to Phosphorus ratios in coastal and estuarine waters",
    y = "N:P ratio",
    caption = "Values represent proportion of N (as mol/l) to P (as mol/l).
    N values are the summed mg/l concentrations of `NH3 filt N`, `Nitrate Filt`, `Nitrite Filt` divided by 14,101.
    P values are the concentrations of `Orthophosphate, Filtered as P` (mg/l) divided by 30,970.
    Blue lines indicate generalised additive model estimates.
    Dashed red lines show indicative N:P ratios of 16 and 32."
    ) +
  ylim(0,200)+
  theme(
    strip.text = element_text(face=2,size = 12),
    strip.background = element_rect(fill = NA, colour=1),
    axis.title.x = element_blank(),
    axis.text.x= element_text(face =2, size =14),
    axis.text.y= element_text(face =2, size =12),
    axis.title.y = element_text(face=2),
    plot.caption = element_text(face=1,size=12),
    plot.title = element_text(face=2,size=16),
  )

# ID wacky values ####
df_PN_w %>% 
  mutate(
    N = rowSums(across(c(`NH3 filt N`, `Nitrate Filt`, `Nitrite Filt`)), na.rm = TRUE),
    N_to_P = N/OrthophsFilt,
    yday = lubridate::yday(sample_date),
  ) %>% 
  filter(sample_date>"2006-12-31") %>% 
  ggplot(.,aes(
    x = RIVER_BASI,
    # y = `NH3 filt N`
    y = `Nitrate Filt`
  ))+
  geom_boxplot()
