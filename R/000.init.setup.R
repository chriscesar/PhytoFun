# 000.init.setup.R ####
## set up project, join site data to shapefiles

# load packages ####
ld_pkgs <- c("tidyverse","ggplot2","sf","maps",
             "ggpubr", "ggspatial","ggrepel","tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log

tic("DATA IMPORT: Downloading, importing & loading data");print("Downloading, importing & loading data")
# load data source locations ####
source("R/dataFolders.R")

# Check if the "data" folder exists in the project directory
if (!file.exists("data")) {
  # If it doesn't exist, create the "data" folder
  dir.create("data")
  print("The 'data' folder has been created.")
} else {
  print("The 'data' folder already exists.")
}

## Check if data file already exists
data_csv <- "data/PHYT_OPEN_DATA_TAXA.csv"
if (!file.exists(data_csv)) {
  # Download data from Open Gov if CSV is missing
  url <- "https://environment.data.gov.uk/ecology/explorer/downloads/PHYT_OPEN_DATA_TAXA.zip"
  data_file <- "data/phyto.zip"
  download.file(url,dest=data_file)
  print("Downloaded data file: phyto.zip")
} else {
  print("Data file (PHYT_OPEN_DATA_TAXA.csv) already exists in the 'data' folder.")
}

## Unzip data only if downloaded and CSV doesn't exist (prevents overwriting)
if (file.exists(data_file) & !file.exists(data_csv)) {
  unzip(data_file, exdir = "data")
  print("Unzipped data file.")
} else {
  print("Skipping unzip as data file wasn't downloaded or CSV already exists.")
}

# ### load WER Cycle3 Trac Shapefile
# base_WBs <- st_read(paste0(GISfol,"C3_WFD_Waterbody/","EnglandTRAC_C3_ReducedFields.shp"))
# base_WBs_df <- fortify(base_WBs)
# plot(base_WBs)

# load data ####
df_phyto0 <- read.csv("data/PHYT_OPEN_DATA_TAXA.csv")
# convert to spatial df
df_phyto0_sf <- sf::st_as_sf(df_phyto0,
                            coords = c("SITE_FULL_EASTING","SITE_FULL_NORTHING"),
                            crs = 27700)

base_WBs <- sf::read_sf(paste0(GISfol,"C3_WFD_Waterbody/","EnglandTRAC_C3_ReducedFields.shp"))
toc(log=TRUE)

# check CRS
# sf::st_crs(base_WBs)

tic("Trim & widen data, convert to SF")
### Widen data to 1 row per sample for joining
df_phyto0 %>%
  as_tibble(.) %>% 
  filter(.,SAMPLE_METHOD_DESCRIPTION=="MARINE: Surface water sample") %>% 
  mutate(.,SAMPLE_REASON = str_replace(SAMPLE_REASON, "National Monitoring ", "National Monitoring")) %>% 
  filter(.,SAMPLE_REASON=="National Monitoring") %>% 
  mutate(taxonKey_taxQual=paste0(TAXON_LIST_ITEM_KEY,"_",TAXON_QUALIFIER_DESC)) %>% 
  mutate(CELLS_LITRE = ifelse(is.na(CELLS_LITRE), 0, CELLS_LITRE),
         COLS_LITRE = ifelse(is.na(COLS_LITRE), 0, COLS_LITRE)) %>% 
  mutate(abund_per_l=(CELLS_LITRE+COLS_LITRE)) %>%
  dplyr::select(
    AGENCY_AREA,REPORTING_AREA,SEA_AREA,WATERBODY_TYPE,
    WFD_WATERBODY_ID,WATER_BODY,SITE_ID,
    SITE_FULL_EASTING,SITE_FULL_NORTHING,
    WIMS_SITE_ID, SAMPLE_ID,
    SAMPLE_DATE,
    taxonKey_taxQual,abund_per_l
    ) %>% 
  group_by(across(-abund_per_l)) %>%
  summarise(abund_per_l=sum(abund_per_l),.groups = "drop") %>% 
  ungroup() %>% 
  pivot_wider(.,names_from = taxonKey_taxQual, values_from = abund_per_l,values_fill = 0) -> df_phytow

### format dates
df_phytow$SAMPLE_DATE <- as.Date(df_phytow$SAMPLE_DATE, format = "%d/%m/%Y")
df_phytow$Year <- year(df_phytow$SAMPLE_DATE)

df_phytow_sf <- sf::st_as_sf(df_phytow,
                             coords = c("SITE_FULL_EASTING","SITE_FULL_NORTHING"),
                             crs = 27700)
saveRDS(df_phytow, file="data/Phyto_WB_wide.rdat")
saveRDS(df_phytow_sf, file="data/Phyto_WB_wide_SF.rdat")

toc(log=TRUE)

# join taxon and spatial data ####
### not necessary: use WBID info in gov.uk download

# tic("join taxon and spatial data"); print("Joining taxon and spatial data")
# joined <- sf::st_join(df_phytow_sf, base_WBs, join = st_within)
# saveRDS(joined, file="data/Phyto_WB_Join.rdat")
# toc(log=TRUE)

(x <- unlist(tic.log()))
