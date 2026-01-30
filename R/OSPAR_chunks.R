# OSPAR_chunks.R ####
# split Phytoplankton data for OSPAR into defined size chunks
## aim: keep individual sample & year data together

# load packages ####
ld_pkgs <- c("tidyverse","tictoc","data.table")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog() ##clear log
tictoc::tic("TOTAL")
# ####################################
# # look into splitting into specific-sized chunks ####
# 
# # Sample size check
# sample_file <- tempfile()
# write.csv(head(dfout_exp, 1000), sample_file, row.names = FALSE)
# sample_size <- file.info(sample_file)$size
# row_size <- sample_size / 1000
# rows_per_file <- floor((100 * 1024^2) / row_size)  # 100 MB limit

tictoc::tic("Load data")
# load data ####
## generated in GenerateOSPAR.v2.R
df <- readRDS(
  #file = "outputs/20251201_OSPAR_Phyto_export_MASTER.Rdat"
  # file = "outputs/20251210_OSPAR_Phyto_export_MASTER.Rdat"
  # file = "outputs/OSPAR/20251219_OSPAR_Phyto_export_MASTER.Rdat"
  # file = "outputs/OSPAR/20260107_OSPAR_Phyto_export_MASTER.Rdat"
  # file = "outputs/OSPAR/20260108_OSPAR_Phyto_export_MASTER.Rdat"
  # file = "outputs/OSPAR/20260109_OSPAR_Phyto_export_MASTER.Rdat"
  # file = "outputs/OSPAR/20260112_OSPAR_Phyto_export_MASTER.Rdat"
  file = "outputs/OSPAR/20260116_OSPAR_Phyto_export_MASTER.Rdat"
  )
tictoc::toc(log = TRUE)

tictoc::tic("Group by year")
# Group by MYEAR ####
groups <- df %>%
  group_by(
    # SMPNO,
    MYEAR
    ) %>%
  group_split()
tictoc::toc(log = TRUE)

tictoc::tic("Estimate file size by year")
# Estimate size of each group in bytes ####
group_sizes <- sapply(groups, function(g) {
  # Rough estimate: number of characters when written as CSV
  sum(nchar(capture.output(write.csv(g, row.names = FALSE))))
})

## Convert to MB ####
group_sizes_mb <- group_sizes / (1024 * 1024)
tictoc::toc(log = TRUE)

tictoc::tic("Split into chunks of ~50 MB")
# Split into chunks of ~50 MB ####
chunks <- list()
current_chunk <- list()
current_size <- 0
chunk_index <- 1

for (i in seq_along(groups)) {
  if (current_size + group_sizes_mb[i] > 48 && length(current_chunk) > 0) {
    chunks[[chunk_index]] <- bind_rows(current_chunk)
    chunk_index <- chunk_index + 1
    current_chunk <- list()
    current_size <- 0
  }
  current_chunk[[length(current_chunk) + 1]] <- groups[[i]]
  current_size <- current_size + group_sizes_mb[i]
}

# Add last chunk
if (length(current_chunk) > 0) {
  chunks[[chunk_index]] <- bind_rows(current_chunk)
}
tictoc::toc(log=TRUE)

tictoc::tic("Write data for export")
# Write chunks to UTF-8 CSV files
for (j in seq_along(chunks)) {
  file_name <- paste0("outputs/OSPAR/",format(Sys.Date(), format="%Y%m%d"),"EA_OSPAR_phytoplankton_", j, ".csv")
  write.csv(chunks[[j]], file_name, row.names = FALSE, fileEncoding = "UTF-8")
}
tictoc::toc(log=TRUE)
tictoc::toc(log=TRUE)

unlist(tictoc::tic.log())
