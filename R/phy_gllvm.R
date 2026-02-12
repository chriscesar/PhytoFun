# phy_gllvm.R ####
## Initial assessments of phyto and nutrient data

# load packages & metadata ####
ld_pkgs <- c("tidyverse",
             "tictoc",
             "gllvm",
             "patchwork"
             )
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
source("R/000setup.R")
source("R/helperFunctions.R")

tictoc::tic.clearlog()

rm(ld_pkgs)

# load data####
df0 <- readRDS(
        file = "//prodds.ntnl/Shared/AN/KFH/Groups/N_Marine/Temporary only - delete in 24 days/For Chris/DATA/Phyto_longTerm/phyto_wims.Rdat")
tictoc::toc(log=TRUE)

names(df0$meta)

## identify dates after 2019
kp <- df0$meta$date > "2019-01-01"

## keep only those
df <- lapply(df0, function(x) x[kp, ])
rm(kp)

# select water bodies of interest ####
wbs <- c(
  "Mersey Mouth", "SOUTHAMPTON WATER", "TEES", "Solent",
  "CARRICK ROADS INNER", "Wash Outer", "POOLE HARBOUR", "MERSEY",
  "TAW / TORRIDGE", "CAMEL", "GREAT OUSE", "Kent South", "EXE",
  "Northumberland North", "RIBBLE", "DART",
  "BURE & WAVENEY & YARE & LOTHING", "LUNE",
  "Farne Islands to Newton Haven","THAMES LOWER"
  )

kp_wb <- df$meta$WB_NAME %in% wbs
## keep only those of interest
df <- lapply(df, function(x) x[kp_wb, ])
rm(kp_wb)

# Try with 1 WB ####
# kp_wb <- "Mersey Mouth"
# kp_wb <- "SOUTHAMPTON WATER"
# kp_wb <- "TEES"
# kp_wb <- "Solent" ### causes issues ###
# kp_wb <- "CARRICK ROADS INNER" ### causes issues ###
# kp_wb <- "Wash Outer"
# kp_wb <- "POOLE HARBOUR"
# kp_wb <- "MERSEY"
# kp_wb <- "TAW / TORRIDGE"
# kp_wb <- "CAMEL"
# kp_wb <- "GREAT OUSE"
# kp_wb <- "Kent South"
# kp_wb <- "EXE"
# kp_wb <- "Northumberland North"
# kp_wb <- "RIBBLE"
# kp_wb <- "DART"
# kp_wb <- "BURE & WAVENEY & YARE & LOTHING"
# kp_wb <- "LUNE"
# kp_wb <- "Farne Islands to Newton Haven"
kp_wb <- "THAMES LOWER"

kp_wb <- df$meta$WB_NAME %in% kp_wb
## keep only those of interest
df_wb <- lapply(df, function(x) x[kp_wb, ])
rm(kp_wb)

## drop SPM & DIN data ####
df_wb$wims <- df_wb$wims %>% 
  ## remove SPM (lots of NA values)
  select(-c(`SPM_mg/l`)) %>% 
  select(-c("DIN"))

# identify WIMS rows with 3 or more NA values
kpwims <- !(rowSums(is.na(df_wb$wims)) >= 3)

df_wb <- lapply(df_wb, function(x) x[kpwims, ])
rm(kpwims)

## rename variables ####
df_wb$wims %>% 
  rename(
    nh4 = "Ammonia_umol/l",
    chla = "Chlorophyll_ug/l",
    sal = "Salinity",
    temp = "Temperature",
    no3 = "Nitrate Nitrogen_umol/l",
    no2 = "Nitrite Nitrogen_umol/l",
    o2_dis_mll = "DO_ml/l",
    # din = "DIN"
    ) -> df_wb$wims

## replace NA values with mean values for respective column ####
## Leave non-numeric cols unchanged
df_wb$wims %>% 
  mutate(across(where(is.numeric),
                ~replace(.,
                         is.na(.),
                         mean(.,
                              na.rm = TRUE)
                         )
                )
         ) -> df_wb$wims_use

## tidy up phytotax data ####
# non-zero column sums
phy_col <- colSums(df_wb$abund)>0
df_wb$abund[,phy_col] -> df_wb$abund_use

## identify taxa with more than n occurrences
n <- 20
vars_kp <- as.vector(colSums(df_wb$abund_use !=0)>n-1)
df_wb$abund_use[,vars_kp] -> df_wb$abund_use

# FIT GLLVMs ####
# By WB ####
## Unconstrained w/ Random ####
### Neg Bin ####
tic("Unconstrained Negative Binomial model")
sDsn <- data.frame(season = df_wb$meta$season)
m_lvm_0 <- gllvm(df_wb$abund_use, # unconstrained model
                 studyDesign = sDsn,
                 # row.eff = ~(1|season),
                 family = "negative.binomial",
                 method = "LA"
                 )
saveRDS(m_lvm_0,
        file = paste0("outputs/models/gllvm_traits_uncon_negbin_",
                      vegan::make.cepnames(unique(df_wb$meta$WB_NAME)),
                      ".Rdat"))
toc(log=TRUE)
# m_lvm_0 <- readRDS("outputs/models/gllvm_traits_uncon_negbin.Rdat")
# plot(m_lvm_0)

## Constrained with random ####
tic("Constrained Negative Binomial model")
## scale variables
df_wb$wims_use_scaled <- as.data.frame(scale(df_wb$wims_use))
### Neg Bin ####
sDsn <- data.frame(season = df_wb$meta$season)
m_lvm_4 <- gllvm(y=df_wb$abund_use, # model with environmental parameters
                 X=df_wb$wims_use_scaled, #scaled
                 formula = ~ nh4+
                   chla+
                   sal+
                   temp+
                   no3+
                   no2+
                   o2_dis_mll,#+
                   # din,
                 studyDesign = sDsn,
                 # row.eff = ~(1|season),
                 method = "LA",
                 family="negative.binomial"
                 )
saveRDS(m_lvm_4,
        file = paste0("outputs/models/gllvm_traits_con_negbin_",
                      vegan::make.cepnames(unique(df_wb$meta$WB_NAME)),
                      ".Rdat"))
toc(log=TRUE)
# m_lvm_4 <- readRDS("outputs/models/gllvm_traits_con_negbin.Rdat")

# tic("Model plots")
# cr <- getResidualCor(m_lvm_4)
# pdf(file = "outputs/figs/m_lvm_4_trt_corrplot.pdf",width=14,height=14)
# corrplot::corrplot(cr, diag = FALSE, type = "lower", method = "square",
#                    tl.srt = 25)
# dev.off()

AIC(m_lvm_0,m_lvm_4)
# anova(m_lvm_0,m_lvm_4)

# pdf(file = "outputs/figs/coef_trt_all_unordered.pdf",width=16,height=8)
# coefplot(m_lvm_4,cex.ylab = 0.7,
#          order=FALSE)
# dev.off()

## Model explore ####
## extract 'significant' model/species terms from model
ci_mod_all <- as.data.frame(confint(m_lvm_4))
ci_mod_var <- ci_mod_all[grep("^X", rownames(ci_mod_all)), ]
rownames(ci_mod_var) <- substring(rownames(ci_mod_var), 7)
ci_mod_var$varTrt <- rownames(ci_mod_var)

sigterms_all <- summary(m_lvm_4)
sigterms_all <- as.data.frame(sigterms_all$Coef.tableX)
sigterms_all$variable <- sub(":.*","",row.names(sigterms_all))
sigterms_all$trt <- sub(".*:","",row.names(sigterms_all))
sigterms_all$varTrt <- rownames(sigterms_all)
sigterms_all <- left_join(sigterms_all, ci_mod_var, by = "varTrt")
sigterms_all$sig <- sigterms_all$`2.5 %`*sigterms_all$`97.5 %`>0

sigterms_sig <- sigterms_all[sigterms_all$`Pr(>|z|)`>0.05,]

## Generate plot ####
plot_list <- list()
sigterms_all$variable <- as.factor(sigterms_all$variable)
ntrt <- length(unique(sigterms_all$trt))-.5
sigterms_all %>% 
  mutate(flag = case_when(
    !sig ~ "NONE",
    sig & Estimate >0 ~ "pos",
    sig & Estimate <0 ~ "neg"
  )) -> sigterms_all

# Iterate over each level of the factor 'trt'
for (level in levels(sigterms_all$variable)) {
  # Subset the data for the current level
  subset_data <- sigterms_all[sigterms_all$variable == level, ]
  
  # Calculate the mean value of 'Estimate' for the current level
  mean_estimate <- mean(subset_data$Estimate)
  
  # Determine the color of the vertical line based on the mean estimate
  line_color <- ifelse(mean_estimate > 0, "blue", ifelse(mean_estimate < 0, "red", "grey"))
  
  # Create a plot for the current level
  current_plot <- ggplot(subset_data,
                         aes(x=Estimate,
                             y=trt,
                             xmin=`2.5 %`,
                             xmax=`97.5 %`,
                             colour=flag,
                             fill=flag)) +
    geom_hline(yintercept = seq(1.5,ntrt,by=1),col="lightgrey",lty=3)+
    geom_vline(xintercept = 0)+
    # Add vertical line for mean estimate
    # geom_vline(xintercept = mean_estimate, color = line_color,linetype="dashed") +
    geom_linerange()+
    labs(title = paste0(level))+
    geom_point(shape=21) +
    scale_y_discrete(limits = rev(levels(as.factor(sigterms_all$trt))))+
    scale_colour_manual(values = c("neg" = "red",#negative
                                   "NONE" = "grey",#null
                                   "pos" = "blue"#positive
    ))+
    scale_fill_manual(values = c("neg" = "red",#negative
                                 "NONE" = "white",#null
                                 "pos" = "blue"#positive
    ))+
    guides(colour="none",
           fill="none")+
    theme(axis.title = element_blank(),
          plot.title = element_text(hjust=0.5),
          )
  
  # Add the current plot to the list
  plot_list[[as.character(level)]] <- current_plot
  
  # Iterate over each plot in the list
  for (i in seq_along(plot_list)) {
    # If it's not the first plot, hide y-axis labels
    if (i > 1) {
      plot_list[[i]] <- plot_list[[i]] + theme(axis.text.y = element_blank())
    }
  }
}

rcov0 <- getResidualCov(m_lvm_0, adjust = 1) # 'null' model
rcov1 <- getResidualCov(m_lvm_4, adjust = 1) # model with env variables #REGION
btwn <- 100 - (rcov1$trace / rcov0$trace*100)
print(paste0("Including environmental parameters in the model explains ",round(btwn,2),"% of the (co)variation in zooplankton abundances"))

### Combine all the individual plots into a single plot ####
(final_plot <- wrap_plots(plotlist = plot_list,
                          ncol = nlevels(sigterms_all$variable))+  # Adjust the number of columns as needed
    plot_annotation(title=paste0("Caterpillar plot of generalised linear latent variable model outputs: ",unique(df_wb$meta$WB_NAME)),
                    # subtitle = bquote("Point estimates & 95% confidence intervals of lifeform-specific model coefficients "~italic(hat(beta)[j])~". Based on phytoplankton taxon densities and scaled water quality parameters"), #scaled
                    subtitle = paste0("Including environmental variables explains ", "<b>",round(btwn,2),"%</b>"," of the (co)variation in taxon abundances compared to the null (taxa-only) model"),
                    caption = paste0(
                      "Based on <b>",nrow(df_wb$meta),"</b> samples gathered between ",format(min(df_wb$meta$date),"%d/%m/%Y")," & ",format(max(df_wb$meta$date),"%d/%m/%Y"),"<br>",
                      "Points indicate model estimate. Bars indicate lifeform 95% confidence intervals which do (grey) or do not (red/blue) include zero","\n","<br>",
                      "Dashed vertical lines indicate mean point estimate values","<br>",
                      "Taxa recorded in fewer than ",n," samples removed from data prior to model estimations","<br>",
                      "Model call: ~",as.character(m_lvm_4$formula)[2],"<br>",
                      "\nDistribution family: ",as.character(m_lvm_4$family),". "
                      # "Random row effects: ",as.character(m_lvm_4$call)[7]
                      ),
                    theme = theme(
                      plot.title = element_text(size = 16, face="bold"),
                      plot.subtitle = ggtext::element_markdown(size=11),
                      plot.caption = ggtext::element_markdown(),
                      )
                    )
 )

# pdf(file = "figs/coef_trt_all_unordered_v2_unscaled.pdf",width=16,height=8)#unscaled
pdf(file = paste0("outputs/figs/coef_trt_all_unordered_v2_scaled_",vegan::make.cepnames(unique(df_wb$meta$WB_NAME)),".pdf"),
    width=16+8,height=8+4) #scaled
print(final_plot)
dev.off()

AIC(m_lvm_0,m_lvm_4)
rm(df_wb, m_lvm_0,m_lvm_4)

# ALL WBs ####
## prep data ####
df_all <- df

### drop SPM data ####
df_all$wims <- df_all$wims %>% 
  select(-c(`SPM_mg/l`))

# identify WIMS rows with 3 or more NA values
kpwims <- !(rowSums(is.na(df_all$wims)) >= 3)

df_all <- lapply(df_all, function(x) x[kpwims, ])
rm(kpwims)

### rename variables ####
df_all$wims %>% 
  rename(
    nh4 = "Ammonia_umol/l",
    chla = "Chlorophyll_ug/l",
    sal = "Salinity",
    temp = "Temperature",
    no3 = "Nitrate Nitrogen_umol/l",
    no2 = "Nitrite Nitrogen_umol/l",
    o2_dis_mll = "DO_ml/l",
    din = "DIN"
  ) -> df_all$wims

### replace NA values with mean values for respective column ####
## Leave non-numeric cols unchanged
df_all$wims %>% 
  mutate(across(where(is.numeric),
                ~replace(.,
                         is.na(.),
                         mean(.,
                              na.rm = TRUE)
                )
  )
  ) -> df_all$wims_use

### tidy up phytotax data ####
# non-zero column sums
phy_col <- colSums(df_all$abund)>0
df_all$abund[,phy_col] -> df_all$abund_use

## identify taxa with more than n occurrences
n <- 20
vars_kp <- as.vector(colSums(df_all$abund_use !=0)>n)
df_all$abund_use[,vars_kp] -> df_all$abund_use

## Unconstrained w/ Random ####
### Neg Bin ####
tic("Unconstrained Negative Binomial model")
sDsn <- data.frame(WB_NAME = df_all$meta$WB_NAME)
m_lvm_0 <- gllvm(df_all$abund_use, # unconstrained model
                 studyDesign = sDsn,
                 row.eff = ~(1|WB_NAME),
                 family = "negative.binomial",
                 method = "LA"
                 )
saveRDS(m_lvm_0,
        file = "outputs/models/gllvm_traits_uncon_negbin_all_WBs.Rdat"
        )
toc(log=TRUE)
# m_lvm_0 <- readRDS("outputs/models/gllvm_traits_uncon_negbin_all_WBs.Rdat")
# plot(m_lvm_0)

## Constrained with random ####
tic("Constrained Negative Binomial model")
## scale variables
df_all$wims_use_scaled <- as.data.frame(scale(df_all$wims_use))
### Neg Bin ####
sDsn <- data.frame(WB_NAME = df_all$meta$WB_NAME)
m_lvm_4 <- gllvm(y=df_all$abund_use, # model with environmental parameters
                 X=df_all$wims_use_scaled, #scaled
                 formula = ~ nh4+chla+sal+temp+no3+no2+o2_dis_mll+din,
                 studyDesign = sDsn,
                 row.eff = ~(1|WB_NAME),
                 method = "LA",
                 family="negative.binomial"
                 )
saveRDS(m_lvm_4,
        file = "outputs/models/gllvm_traits_con_negbin_all_WBs.Rdat"
        )
toc(log=TRUE)
# m_lvm_4 <- readRDS("outputs/models/gllvm_traits_con_negbin_all_WBs.Rdat")
