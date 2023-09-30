# 000.init.setup.R
## set up project

# create data folder:
# Check if the "data" folder exists in the project directory
if (!file.exists("data")) {
  # If it doesn't exist, create the "data" folder
  dir.create("data")
  print("The 'data' folder has been created.")
} else {
  print("The 'data' folder already exists.")
}


## DL data
url <- "https://environment.data.gov.uk/ecology/explorer/downloads/PHYT_OPEN_DATA_TAXA.zip"
dest <- "data/phyto.zip"
download.file(url,dest)

## unzip file
unzip("data/phyto.zip", exdir = "data")
