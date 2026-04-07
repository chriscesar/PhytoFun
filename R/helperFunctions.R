# helperFunctions.R ####
# Full helper: accepts either numeric easting/northing or an OS grid ref string.
ngr_to_wgs84 <- function(ngr) {
  # If 'ngr' is a list/data.frame with easting & northing, convert directly
  if (is.data.frame(ngr) && all(c("easting", "northing") %in% names(ngr))) {
    return(bng_to_wgs84(ngr$easting, ngr$northing))
  }
  
  # If package rnrfa is available, use its osg_parse which handles many formats
  if (requireNamespace("rnrfa", quietly = TRUE)) {
    parsed <- rnrfa::osg_parse(ngr, coord_system = "WGS84")
    # rnrfa::osg_parse(..., coord_system="WGS84") returns lon/lat in WGS84
    return(parsed)
  }
  
  # Otherwise fallback to a compact parser -> easting/northing -> sf transform
  # Minimal parser: handles letter+digits like "SU387148", "SU 38714 14807", "TQ123456"
  ngr_clean <- toupper(gsub("[^A-Z0-9]", "", ngr))
  letter_pair <- substr(ngr_clean, 1, 2)
  digits       <- substring(ngr_clean, 3)
  
  if (nchar(letter_pair) != 2 || nchar(digits) %% 2 != 0) {
    stop("Unrecognised NGR format. Try installing package 'rnrfa' or pass easting/northing.")
  }
  
  # parse numeric part into easting/northing offsets within the 100km square
  half <- nchar(digits) / 2
  e_part <- substr(digits, 1, half)
  n_part <- substr(digits, half + 1, nchar(digits))
  
  # scale to metres depending on precision (e.g., '387' -> 38700 m at 3-digit)
  factor <- 10^(5 - half)
  e_offset <- as.numeric(e_part) * factor
  n_offset <- as.numeric(n_part) * factor
  
  # mapping of 100km grid letters -> 100km easting/northing offsets
  # (based on OS grid: 1st letter selects 500km square, 2nd selects 100km square)
  # We implement the standard mapping as in many OS implementations.
  letter_idx <- function(ch) {
    letters_no_i <- unlist(strsplit("ABCDEFGHJKLMNOPQRSTUVWXYZ", ""))
    which(letters_no_i == ch)
  }
  
  # The 100km grid (for the second letter) uses a 5x5 layout; we follow
  # the established OS scheme:
  letters_no_i <- unlist(strsplit("ABCDEFGHJKLMNOPQRSTUVWXYZ", ""))
  idx1 <- letter_idx(substr(letter_pair, 1, 1))
  idx2 <- letter_idx(substr(letter_pair, 2, 2))
  if (length(idx1) == 0 || length(idx2) == 0) {
    stop("Invalid grid letters in NGR.")
  }
  
  # compute 100km-square origins (see OS algorithm)
  # First letter gives 500km square origin:
  e100k1 <- ((idx1 - 1) %% 5) * 500000
  n100k1 <- (4 - ((idx1 - 1) %/% 5)) * 500000
  
  # Second letter gives 100km square inside that:
  e100k2 <- ((idx2 - 1) %% 5) * 100000
  n100k2 <- (4 - ((idx2 - 1) %/% 5)) * 100000
  
  easting <- e100k1 + e100k2 + e_offset
  northing <- n100k1 + n100k2 + n_offset
  
  # Now transform the BNG easting/northing to WGS84 lat/lon
  bng_to_wgs84(easting, northing)
}

#################
# OSGB Grid Reference to WGS84 Lat/Long Converter
# Requires sf and units packages

library(sf)
library(units)

# Function to convert OSGB grid reference to WGS84 coordinates
osgb_to_wgs84 <- function(grid_ref) {
  
  # Grid square lookup table (bottom-left corners in meters)
  grid_squares <- list(
    "HP" = c(400000, 1200000), "HT" = c(300000, 1100000), "HU" = c(400000, 1100000),
    "HW" = c(100000, 1000000), "HX" = c(200000, 1000000), "HY" = c(300000, 1000000), "HZ" = c(400000, 1000000),
    "NA" = c(0, 900000), "NB" = c(100000, 900000), "NC" = c(200000, 900000), "ND" = c(300000, 900000), "NE" = c(400000, 900000),
    "NF" = c(0, 800000), "NG" = c(100000, 800000), "NH" = c(200000, 800000), "NJ" = c(300000, 800000), "NK" = c(400000, 800000),
    "NL" = c(0, 700000), "NM" = c(100000, 700000), "NN" = c(200000, 700000), "NO" = c(300000, 700000), "NP" = c(400000, 700000),
    "NQ" = c(0, 600000), "NR" = c(100000, 600000), "NS" = c(200000, 600000), "NT" = c(300000, 600000), "NU" = c(400000, 600000),
    "NW" = c(100000, 500000), "NX" = c(200000, 500000), "NY" = c(300000, 500000), "NZ" = c(400000, 500000), "OV" = c(500000, 500000),
    "SB" = c(100000, 400000),
    "SC" = c(200000, 400000), "SD" = c(300000, 400000), "SE" = c(400000, 400000), "TA" = c(500000, 400000),
    "SH" = c(200000, 300000), "SJ" = c(300000, 300000), "SK" = c(400000, 300000), "TF" = c(500000, 300000), "TG" = c(600000, 300000),
    "SM" = c(100000, 200000), "SN" = c(200000, 200000), "SO" = c(300000, 200000), "SP" = c(400000, 200000), "TL" = c(500000, 200000), "TM" = c(600000, 200000),
    "SR" = c(100000, 100000), "SS" = c(200000, 100000), "ST" = c(300000, 100000), "SU" = c(400000, 100000), "TQ" = c(500000, 100000), "TR" = c(600000, 100000),
    "SW" = c(100000, 0), "SX" = c(200000, 0), "SY" = c(300000, 0), "SZ" = c(400000, 0), "TV" = c(500000, 0)
  )
  
  # Extract grid square (first 2 characters)
  grid_square <- substr(grid_ref, 1, 2)
  
  # Extract remaining digits
  digits <- substr(grid_ref, 3, nchar(grid_ref))
  
  # Split digits in half for easting and northing
  n_digits <- nchar(digits)
  if (n_digits %% 2 != 0) {
    stop("Grid reference digits must be even in length")
  }
  
  half_length <- n_digits / 2
  easting_digits <- substr(digits, 1, half_length)
  northing_digits <- substr(digits, half_length + 1, n_digits)
  
  # Convert to numbers and scale based on precision
  precision <- 10^(5 - half_length)  # 5 digits = 1m precision
  local_easting <- as.numeric(easting_digits) * precision
  local_northing <- as.numeric(northing_digits) * precision
  
  # Get grid square origin
  if (!grid_square %in% names(grid_squares)) {
    stop(paste("Unknown grid square:", grid_square))
  }
  
  origin <- grid_squares[[grid_square]]
  
  # Calculate full OSGB coordinates
  full_easting <- origin[1] + local_easting
  full_northing <- origin[2] + local_northing
  
  # Create sf point in British National Grid (EPSG:27700)
  point_bng <- st_sfc(st_point(c(full_easting, full_northing)), crs = 27700)
  
  # Transform to WGS84 (EPSG:4326)
  point_wgs84 <- st_transform(point_bng, crs = 4326)
  
  # Extract coordinates
  coords <- st_coordinates(point_wgs84)
  
  return(data.frame(
    longitude = coords[1],
    latitude = coords[2],
    easting = full_easting,
    northing = full_northing
  ))
}

# Function to convert multiple grid references
convert_multiple_osgb <- function(grid_refs) {
  result <- data.frame()
  
  for (i in 1:length(grid_refs)) {
    tryCatch({
      converted <- osgb_to_wgs84(grid_refs[i])
      converted$grid_ref <- grid_refs[i]
      result <- rbind(result, converted)
    }, error = function(e) {
      warning(paste("Error converting", grid_refs[i], ":", e$message))
    })
  }
  
  return(result)
}

