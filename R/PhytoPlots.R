# PhytoPlots.R ####

# Investigate Solent data

# load packages ####
ld_pkgs <- c("tidyverse","tictoc","stringr","gllvm")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

source("R/000setup.R")
outfol <- "//prodds.ntnl/Shared/AN/KFH/Groups/N_Marine/02 Projects_Tasks/04 E&B/Ecology and Ecosystems/Planktonic_Food_Webs/202602_FigsForReport/"
## define random vector function ####
random_logical_vector <- function(n, k) {
  if (k > n) stop("k cannot be greater than n")
  
  # Create vector with k TRUEs and (n-k) FALSEs
  vec <- c(rep(TRUE, k), rep(FALSE, n - k))
  
  # Randomly shuffle the vector
  sample(vec)
}

# load data ####
tic("Load data")
df0 <- readxl::read_xlsx("outputs/OSPAR/phyto_wims_wb.xlsx",
                         sheet="phyto_wims_wb",guess_max = 10000)
df0_tx_lk <- readxl::read_xlsx("outputs/OSPAR/phyto_wims_wb.xlsx",sheet="tx")
toc(log=TRUE)

tic("Create list object for data")
# split data into meta, wims & biol ####
df0 %>% 
  dplyr::select(STNNO,LATIT.x,LONGI.x,STATN,SDATE,SMPNO,
                LATIT.y,LONGI.y,WB_ID,WB_NAME,WB_CAT,WB_HMWB
                ,WB_TYPOL,OP_CATCH,OP_CATCH1,MGT_CATCH,MGT_CATCH1,EA_AREA,
                RBD,COUNTRY,NITR_DIR,SHELLF_DIR,UWW_DIR,CONSERVATI,
                HABITATS_A,OVERALL__1,WATER_BO_5,WATER_BO_6,WATER_BO_7,
                SHAPE_Leng,SHAPE_Area,EXPORT_DAT,Distance) %>%
  mutate(sample.date = as.Date(as.character(SDATE), format = "%Y%m%d")) %>% 
  relocate(sample.date, .after = SDATE) -> df0_meta

df0 %>% 
  dplyr::select(STNNO,SMPNO,"Ammonia_umol/l","Chlorophyll_ug/l",
                "Salinity","Temperature","Nitrate Nitrogen_umol/l",
                "Nitrite Nitrogen_umol/l","SPM_mg/l","DO_ml/l",
                "DIN") -> df0_wims

df0 %>% 
  dplyr::select(STNNO,SMPNO,
                "Cyclotella atomus":Micractinium) -> df0_tx

## swap species for LFs
df0_tx %>% tidyr::pivot_longer(cols = -c(STNNO,SMPNO)) %>% 
  left_join(.,df0_tx_lk,
            by = c("name" = "TAX Report")) -> df0_tx_l

df0_tx_l %>% 
  dplyr::select(STNNO,SMPNO, phyto_lf, value) %>% 
  group_by(across(-value)) %>% 
  summarise(sum=sum(value, na.rm=TRUE),.groups = "drop") %>% 
  ungroup() %>% 
  pivot_wider(names_from = phyto_lf,values_from = sum) -> df0_tx_w

# create list object ###
df <- list("meta0"=df0_meta,
           "abund0"=df0_tx_w,
           "wims0"=df0_wims)

rm(df0_meta,
   df0_tx_w,
   df0_wims, df0_tx,df0_tx_l,df0_tx_lk)

toc(log=TRUE)

tic("prep for model")
# prep for model ####
df$wims0 %>% 
  rename(
    nh4 = "Ammonia_umol/l",
    chla = "Chlorophyll_ug/l",
    no3 = "Nitrate Nitrogen_umol/l",
    no2 = "Nitrite Nitrogen_umol/l",
    din = "DIN",
    o2_dis_mll = "DO_ml/l",
    sal_ppt = "Salinity",
    tempC = "Temperature",
    spm = "SPM_mg/l"
  ) %>% #names() %>% 
  dplyr::select(-STNNO,-SMPNO) -> df$wims

row_flag <- df$wims %>%
  transmute(flag = !if_all(everything(), is.na)) %>%
  pull(flag)

# keep only samples which are not all NA values
df$meta0[row_flag,] -> df$meta_use
df$wims[row_flag,] -> df$wims_use
df$abund0[row_flag,] -> df$abund_use

## replace NA values with mean values for respective column.
## Leave non-numeric cols unchanged
df$wims_use %>% 
  mutate(across(where(is.numeric),
                ~replace(.,
                         is.na(.),
                         mean(.,
                              na.rm = TRUE)
                )
  )
  ) -> X

### create scaled version for comparison of effects on model
X_scaled <- as.data.frame(scale(X))

# define species data
Y <- df$abund_use %>% dplyr::select(-STNNO,-SMPNO)
toc(log=TRUE)

# Generate subsample of data for analysis ####
# set seed for replicability
set.seed(78678)
n <- nrow(Y) # total number of rows in data
k <- 1500 # how many samples should I extract?
rowkeep <- random_logical_vector(n = n, k=k)

## subsample data
X_scaled_sub <- X_scaled[rowkeep,]

# define minimum number of samples for species to be retained
n_samp <- 100
Y_sub <- Y[rowkeep,]
# Y_sub <- Y_sub %>% select(where(~ sum(.x) != 0))#remove empty columns
keep <- colSums(Y_sub != 0, na.rm = TRUE) >= n_samp
Y_sub <- Y_sub[, keep, drop = FALSE]

meta_sub <- df$meta_use[rowkeep,]

# remove empty taxa from subsample

tic("Fit Unconstrained model")
## Fit unconstrained GLLVM ####
## set model options
runs <- 5

### Negative Binomial ####
sDsn <- data.frame(WB = meta_sub$WB_NAME)
m_lvm_0 <- gllvm(as.matrix(Y_sub), # unconstrained model
                 studyDesign = sDsn,
                 row.eff = ~(1|WB),
                 family = "negative.binomial",
                 starting.val="res",
                 n.init = runs, #re-run model to get best fit
                 trace=TRUE,
                 # seed = 123,
                 num.lv = 2
)

saveRDS(
  m_lvm_0,
  file="outputs/models/gllvm_longterm/gllvm_phyto_abund_negbin_uncond.Rdat"
  )#404.73sec
toc(log=TRUE)

m_lvm_0 <- readRDS(
  file="outputs/models/gllvm_longterm/gllvm_phyto_abund_negbin_uncond.Rdat"
  )

tic("Fit Constrained model")
## Fit Constrained GLLVM ####
### Negative Binomial ####
sDsn <- data.frame(WB = meta_sub$WB_NAME)
m_lvm_4 <- gllvm(as.matrix(Y_sub), # unconstrained model
                 X = X_scaled_sub,
                 formula = ~chla + nh4 + no3 + no2 + spm + sal_ppt +
                   o2_dis_mll + tempC,
                 studyDesign = sDsn,
                 row.eff = ~(1|WB),
                 family = "negative.binomial",
                 starting.val="res",
                 n.init = runs, #re-run model to get best fit
                 trace=TRUE,
                 # seed = 123,
                 num.lv = 2
)
saveRDS(
  m_lvm_4,
  file="outputs/models/gllvm_longterm/gllvm_phyto_abund_negbin_cond.Rdat"
)

m_lvm_4 <- readRDS(
  "outputs/models/gllvm_longterm/gllvm_phyto_abund_negbin_cond.Rdat"
  )
toc(log=TRUE)

cr <- getResidualCor(m_lvm_4)
  
corrplot::corrplot(cr, diag = FALSE, type = "lower", method = "square",
                     tl.srt = 25)

## compare models
AIC(m_lvm_0,m_lvm_4)

coefplot(m_lvm_4)
saveRDS(unlist(tictoc::tic.log()),
        file="outputs/models/gllvm_longterm/tic_log.Rdat"
        )

#########

## GLLVM plots ####

for(i in 1:ncol(m_lvm_4$X.design)){
  pdf(file = paste0("outputs/figs/gllvm_all_data/coef_nb_trt_",i,".pdf"),width = 7, height = 14)
  coefplot(m_lvm_4,mfrow = c(1,1),which.Xcoef = i, cex.ylab = 0.6,
           main=colnames(m_lvm_4$X.design)[i])
  dev.off()
}

### extract model terms for plotting ####
mod_coefs <- as.data.frame(m_lvm_4$params$Xcoef)
mod_coefs$LF <- row.names(mod_coefs)
mod_coefs <- mod_coefs %>% 
  pivot_longer(.,cols=!LF,names_to = "coefficient", values_to = "Estimate")
sdXcoef <- as.data.frame(m_lvm_4$sd$Xcoef[,  drop = FALSE])
sdXcoef$LF <- row.names(sdXcoef)
sdXcoef <- sdXcoef %>% 
  pivot_longer(.,cols=!LF,names_to = "coefficient", values_to = "sd")

mod_coefs <- dplyr::left_join(mod_coefs,sdXcoef, by=c("LF","coefficient"))
mod_coefs <- mod_coefs %>% 
  mutate(.,lower = Estimate-1.96*sd,
         upper = Estimate+1.96*sd) %>% 
  mutate(.,varTrt=paste0(coefficient,"_",LF))

sigterms_all <- summary(m_lvm_4)
sigterms_all <- as.data.frame(sigterms_all$Coef.tableX)
sigterms_all$variable <- sub(":.*","",row.names(sigterms_all))
sigterms_all$trt <- sub(".*:","",row.names(sigterms_all))
sigterms_all$varTrt <- rownames(sigterms_all)
sigterms_all$varTrt <- gsub(":","_",sigterms_all$varTrt)

mod_coefs <- dplyr::left_join(mod_coefs,
                              sigterms_all[,c("varTrt","Std. Error","z value","Pr(>|z|)")],
                              by="varTrt")
rm(sigterms_all,sdXcoef)
## create flag if CI doesn't cross 0
mod_coefs$sig <- mod_coefs$lower*mod_coefs$upper>0

### plot! ####
plot_list <- list()
mod_coefs$coefficient <- as.factor(mod_coefs$coefficient)
# sigterms_all$variable <- as.factor(sigterms_all$variable)
ntrt <- length(unique(mod_coefs$LF))-.5
mod_coefs %>% 
  mutate(flag = case_when(
    !sig ~ "NONE",
    sig & Estimate >0 ~ "pos",
    sig & Estimate <0 ~ "neg"
  )) %>% 
  mutate(clr = case_when(
    !sig ~ "grey",
    sig & Estimate >0 ~ "blue",
    sig & Estimate <0 ~ "red"
  )) -> mod_coefs

# Iterate over each level of the factor 'LF'
for (level in levels(mod_coefs$coefficient)) {
  # Subset the data for the current level
  subset_data <- mod_coefs[mod_coefs$coefficient == level, ]
  #subset_data <- mod_coefs[mod_coefs$coefficient == "sal_ppt", ]
  
  # Calculate the mean value of 'Estimate' for the current level
  mean_estimate <- mean(subset_data$Estimate)
  
  # Determine the color of the vertical line based on the mean estimate
  line_color <- ifelse(mean_estimate > 0, "blue",
                       ifelse(mean_estimate < 0, "red", "grey"))
  
  # Create a plot for the current level
  current_plot <- ggplot(subset_data,
                         aes(x=Estimate, y=LF,
                             xmin=lower,
                             xmax=upper,
                             colour=clr,
                             fill=clr)) +
    geom_hline(yintercept = seq(1.5,ntrt,by=1),col="lightgrey",lty=3)+
    geom_vline(xintercept = 0, lty=2, col="darkgrey")+
    # Add vertical line for mean estimate
    # geom_vline(xintercept = mean_estimate, color = line_color,linetype="dashed") +
    geom_linerange()+
    labs(title = paste0(level))+
    geom_point() +
    scale_y_discrete(limits = rev(levels(as.factor(mod_coefs$LF))))+
    scale_colour_identity()+
    guides(colour="none",
           fill="none")+
    theme(axis.title = element_blank(),
          plot.title = element_text(hjust=0.5,size=14,face=2),
          axis.text = element_text(face=2),
          axis.text.y = element_text(size=12),
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
print(paste0("Including environmental parameters in the model explains ",round(btwn,1),"% of the (co)variation in phytoplankton abundances"))
AIC(m_lvm_0,m_lvm_4)

# Combine all the individual plots into a single plot
(final_plot <-
    patchwork::wrap_plots(
      plotlist = plot_list,
      ncol = nlevels(mod_coefs$coefficient)
    ) +
    patchwork::plot_annotation(
      title="Caterpillar plot of generalised linear latent variable model outputs",
      subtitle = paste0("Including environmental variables explains <b>",round(btwn,1),"% </b>of the (co)variation in phytoplankton lifeform carbon content compared to the null (lifeforms-only) model"),
      caption = paste0(
        "Colours indicate lifeform 95% confidence intervals which do (grey) or do not (red/blue) include zero<br>",
        "Taxa recorded in fewer than <b>",n_samp,"</b> samples removed from data prior to model estimations<br>",
        "<i>Model call:</i> ~",as.character(m_lvm_4$formula)[2],"",
        "<br><i>Distribution family:</i> ",as.character(unique(m_lvm_4$family)),".<br>",
        "<i>Random row effects: </i>",as.character(m_lvm_4$call)[8],". <i>Number of model iterations: </i> ",m_lvm_4$n.init,"<br>",
        "Based on a random subsample of <b>",length(unique(meta_sub$SMPNO)),"</b> samples gathered across <b>",length(unique(meta_sub$WB_NAME)),"</b> water bodies between <b>",format(min(meta_sub$sample.date), "%d/%m/%Y")," & ",
        format(max(meta_sub$sample.date), "%d/%m/%Y"),"</b>",
        "<br><i>AIC(null_mod):</i>",round(AIC(m_lvm_0),2),
        "<br><i>AIC(env_mod):</i>",round(AIC(m_lvm_4),2)
      ),
      # >>> Theme ONLY for the annotation <<<
      theme = theme(
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = ggtext::element_markdown(size = 14),
        plot.caption = ggtext::element_markdown(size = 11)
      )
    ))

png(file = paste0(outfol,"phyto_gllvm_coef_all.png"),
    width=18*ppi, height=12*ppi, res=ppi)
print(final_plot)
dev.off()


par(mfrow=c(1,1))
VP <- varPartitioning(m_lvm_4)
plotVarPartitioning(VP,col=palette(hcl.colors(10,"Roma")))

# plot with ggplot2
as.data.frame(VP$PropExplainedVarSp) %>%
  mutate(name = row.names(.)) %>% 
  pivot_longer(.,cols = nh4:"Row random effect: WB",
               names_to = "var",values_to = "val") %>% 
  mutate(var = as.factor(var)) %>% 
  mutate(var = fct_relevel(var,
                           "Row random effect: WB",
                           "LV1","LV2")) -> VP_plot

(VP_plot %>% 
    ggplot(., aes(fill = var, y=val, x= name)) + 
    geom_bar(col=1,position = "fill", stat= "identity")+
    scale_x_discrete(limits=rev)+
    coord_flip()+
    labs(title = "Variance partitioning of phtoplankton life form abundances",
         y="Variance explained",
         caption = paste0("Samples gathered between ",
                          min(meta_sub$sample.date)," & ",
                          max(meta_sub$sample.date)))+
    scale_fill_manual(values = c("Row random effect: WB" = "#000000",
                                 "LV1" = "grey50",
                                 "LV2" = "grey70",
                                 chla = cbPalette2[3],
                                 nh4 = cbPalette2[4],
                                 no3 = cbPalette2[5],
                                 no2 = cbPalette2[6],
                                 spm = cbPalette2[7],
                                 sal_ppt = cbPalette2[8],
                                 o2_dis_mll = cbPalette2[9],
                                 tempC = cbPalette2[10]
    ))+
    # scale_colour_manual(values=c("red","red","red",rep("black",7)))+
    theme(
      axis.title.y = element_blank(),
      axis.text = element_text(face = 2),
      axis.title.x = element_text(face = 2),
      legend.text = element_text(face = 2),
      plot.title = element_text(face = 2),
      plot.caption = element_text(face = 2),
      legend.title = element_blank(),
    ) -> pl_VP)
##########################################
# TO DO ####

# re-run gllvms on lifeforms after calculating month-year means across all data