# PhytoPlots.R ####

# Investigate Solent data

# load packages ####
ld_pkgs <- c("tidyverse","tictoc","stringr","gllvm")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

source("R/000setup.R")

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

Y <- df$abund_use %>% dplyr::select(-STNNO,-SMPNO)

toc(log=TRUE)

# Generate subsample of data for analysis ####
# set seed for replicability
set.seed(pi)
n <- nrow(Y) # total number of rows in data
k <- 700 # how many samples should I extract?
rowkeep <- random_logical_vector(n = n, k=k)

## subsample data
X_scaled_sub <- X_scaled[rowkeep,]
Y_sub <- Y[rowkeep,]
Y_sub <- Y_sub %>% select(where(~ sum(.x) != 0))#remove empty columns
meta_sub <- df$meta_use[rowkeep,]

# remove empty taxa from subsample

tic("Fit Unconstrained model")
## Fit unconstrained GLLVM ####
## set model options
runs <- 1

### Negative Binomial ####
sDsn <- data.frame(WB = meta_sub$WB_NAME)
m_lvm_0 <- gllvm(as.matrix(Y_sub), # unconstrained model
                 studyDesign = sDsn,
                 row.eff = ~(1|WB),
                 family = "negative.binomial",
                 starting.val="res",
                 n.init = runs, #re-run model to get best fit
                 # trace=TRUE,
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
                 # trace=TRUE,
                 # seed = 123,
                 num.lv = 2
)
saveRDS(
  m_lvm_4,
  file="outputs/models/gllvm_longterm/gllvm_phyto_abund_negbin_cond.Rdat"
)
toc(log=TRUE)

saveRDS(unlist(tictoc::tic.log()),
        file="outputs/models/gllvm_longterm/tic_log.Rdat"
        )
