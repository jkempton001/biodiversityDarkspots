# GBIF Occurrence Data Download for Vascular Plants in New Guinea

# Load required packages ----
library(tidyverse)
library(rgbif)
library(magrittr)
library(bit64)
library(keyring)

get_gbif_credentials <- function(service_user = "jkempton001") {
  
  # Helper: Ensure a credential exists, else prompt to set it
  ensure_credential <- function(service_name) {
    existing <- tryCatch(
      key_get(service_name, username = service_user),
      error = function(e) NULL
    )
    
    if (is.null(existing)) {
      message(sprintf("No credential found for '%s'. Please enter it now.", service_name))
      key_set(service_name, username = service_user)
      existing <- key_get(service_name, username = service_user)
    }
    
    return(existing)
  }
  
  list(
    user     = ensure_credential("gbif_user"),
    password = ensure_credential("gbif_password"),
    email    = ensure_credential("gbif_email")
  )
}


# Load credentials securely
creds <- get_gbif_credentials()
GBIF_USER <- creds$user
GBIF_PWD <- creds$password
GBIF_EMAIL <- creds$email

cat("GBIF credentials loaded securely from keychain\n")

# Define administrative regions ----
# Using GADM codes for precise geographic filtering
# Note: Island-level filtering on GBIF returns inaccurate results, 
# so we filter by administrative districts instead

regions <- list(
  # New Guinea - Indonesian provinces
  papua_barat = "IDN.22_1",
  papua = "IDN.23_1",
  
  # New Guinea - PNG mainland provinces (17 provinces)
  png_mainland = c("PNG.18_1", "PNG.5_1", "PNG.11_1", "PNG.7_1", "PNG.21_1", 
                   "PNG.9_1", "PNG.10_1", "PNG.3_1", "PNG.19_1", "PNG.6_1", 
                   "PNG.14_1", "PNG.8_1", "PNG.22_1", "PNG.17_1", "PNG.2_1", 
                   "PNG.15_1", "PNG.13_1")
)
 

# Create standardized download function ----
download_gbif_occurrences <- function(region_codes, region_name) {
  # Build predicates for the download
  predicates <- list(
    pred("taxonKey", 7707728),  # Tracheophyta (vascular plants)
    pred("basisOfRecord", "PRESERVED_SPECIMEN"),  # herbarium specimens only
    pred("hasGeospatialIssue", FALSE),  # remove records with location issues
    pred("hasCoordinate", TRUE),  # coordinates required
    pred("occurrenceStatus", "PRESENT")  # presence records only
  )
  
  # Add GADM predicates - handle both single and multiple regions
  if (length(region_codes) == 1) {
    predicates <- append(predicates, list(pred("gadm", region_codes)), after = 1)
  } else {
    # For multiple regions (like PNG mainland), add each as separate predicate
    gadm_preds <- map(region_codes, ~pred("gadm", .x))
    predicates <- append(predicates, list(pred_in("gadm", region_codes)), after = 1)
  }
  
  # Execute download
  download_key <- do.call(occ_download, c(
    predicates,
    format = "DWCA",  # Darwin Core Archive for reproducibility and citation
    user = GBIF_USER,
    pwd = GBIF_PWD,
    email = GBIF_EMAIL
  ))
  
  cat("Download initiated for", region_name, "- Key:", download_key, "\n")
  return(download_key)
}

# Execute downloads for each region ----
# Note: Each region must be downloaded separately due to rgbif limitations
# This approach ensures precise geographic filtering for each administrative unit

download_keys <- list()

cat("Initiating GBIF downloads for all regions...\n")

# New Guinea regions
download_keys$papua_barat <- download_gbif_occurrences(regions$papua_barat, "Papua Barat")
download_keys$papua <- download_gbif_occurrences(regions$papua, "Papua")
download_keys$png_mainland <- download_gbif_occurrences(regions$png_mainland, "PNG Mainland")

# Check download status ----
check_download_status <- function(download_keys) {
  cat("\nChecking download status...\n")
  for (region in names(download_keys)) {
    status <- occ_download_meta(download_keys[[region]])$status
    cat(region, ":", status, "\n")
  }
}

# Uncomment to check status
check_download_status(download_keys)

# Import downloaded data ----
# Replace these with your actual download keys once downloads are complete
actual_download_keys <- list(
  papua_barat = "",
  papua = "", 
  png_mainland = ""
)

# Function to import and process regional data
import_and_process_region <- function(download_key, island, country) {
  # Import data
  data <- occ_download_get(download_key, overwrite = TRUE) %>%
    occ_download_import(na.strings = c("", NA))
  
  # Select relevant columns and add geographic identifiers
  processed_data <- data %>%
    select(gbifID, license, institutionCode, collectionCode, occurrenceID, 
           catalogNumber, recordNumber, recordedBy, eventDate, habitat, 
           eventRemarks, locality, decimalLatitude, decimalLongitude, 
           coordinateUncertaintyInMeters, coordinatePrecision, identifiedBy, 
           scientificName, kingdom, phylum, class, order, family, genus, 
           taxonRank, taxonomicStatus, elevation, elevationAccuracy, issue, 
           taxonKey, acceptedTaxonKey, speciesKey, species, acceptedScientificName, 
           verbatimScientificName, iucnRedListCategory, countryCode, 
           individualCount, year, basisOfRecord, datasetName, establishmentMeans, 
           references) %>%
    # COMPREHENSIVE data type standardization - convert ALL columns to expected types
    mutate(
      # Convert ALL potentially problematic columns to character first, then to proper type
      across(everything(), as.character)
    ) %>%
    # Now convert to proper types where needed
    mutate(
      # Numeric columns 
      decimalLatitude = as.numeric(decimalLatitude),
      decimalLongitude = as.numeric(decimalLongitude),
      coordinateUncertaintyInMeters = as.numeric(coordinateUncertaintyInMeters),
      coordinatePrecision = as.numeric(coordinatePrecision),
      elevation = as.numeric(elevation),
      elevationAccuracy = as.numeric(elevationAccuracy),
      individualCount = as.numeric(individualCount),
      year = as.numeric(year),
      
      # Keep ID columns as character (safer than bit64)
      gbifID = as.character(gbifID),
      taxonKey = as.character(taxonKey),
      acceptedTaxonKey = as.character(acceptedTaxonKey),
      speciesKey = as.character(speciesKey),
      
      # All other columns remain as character (including issue, license, eventDate, etc.)
    ) %>%
    add_column(ISLAND = island, COUNTRY = country)
  
  cat("Processed", island, "-", country, ":", nrow(processed_data), "records\n")
  return(processed_data)
}

# Define island and country assignments
region_metadata <- list(
  papua_barat = list(island = "NewGuinea", country = "Indonesia"),
  papua = list(island = "NewGuinea", country = "Indonesia"),
  png_mainland = list(island = "NewGuinea", country = "PNG")
)

# Import and process all regional datasets
cat("Importing and processing regional datasets...\n")

# Process each dataset individually with error handling
regional_datasets <- list()
for (region_name in names(download_keys)) {
  tryCatch({
    regional_datasets[[region_name]] <- import_and_process_region(
      download_keys[[region_name]], 
      region_metadata[[region_name]]$island,
      region_metadata[[region_name]]$country
    )
  }, error = function(e) {
    cat("Error processing", region_name, ":", e$message, "\n")
    regional_datasets[[region_name]] <<- NULL
  })
}

# Remove any NULL entries (failed downloads)
regional_datasets <- regional_datasets[!sapply(regional_datasets, is.null)]

cat("Successfully processed", length(regional_datasets), "regional datasets\n")

# Combine all regional datasets
cat("Combining datasets from all regions...\n")

# Additional safety check - examine data types before combining
cat("Checking data types for consistency...\n")
for (i in seq_along(regional_datasets)) {
  region_name <- names(regional_datasets)[i]
  cat("Dataset", i, "(", region_name, "):\n")
  
  # Check for problematic columns that have caused issues
  problematic_cols <- c("issue", "license", "eventDate")
  for (col in problematic_cols) {
    if (col %in% names(regional_datasets[[i]])) {
      col_type <- class(regional_datasets[[i]][[col]])[1]
      cat("  ", col, ":", col_type, "\n")
    }
  }
}

# Combine datasets (should work smoothly now with standardized types)
gbif_combined <- bind_rows(regional_datasets)

cat("Successfully combined", nrow(gbif_combined), "records from", length(regional_datasets), "regions\n")

# Optional: Check for any remaining data type issues
cat("Data summary by column type:\n")
gbif_combined %>% 
  summarise(across(everything(), ~class(.)[1])) %>% 
  pivot_longer(everything(), names_to = "column", values_to = "type") %>%
  count(type) %>%
  print()

# Clean and standardize GBIF data ----
gbif_cleaned <- gbif_combined %>%
  # Parse accepted scientific name into components
  # Note: Many records will only have genus, or genus + species, without authority
  # This generates expected warnings for records missing some components
  separate(acceptedScientificName, 
           into = c("gen", "sp", "authority"), 
           sep = " ", 
           extra = "merge",
           fill = "right") %>%  # Fill missing pieces with NA (suppresses warnings)
  # Standardize column names to match internal dataset
  rename(collector = recordedBy,
         number = recordNumber,
         lat = decimalLatitude,
         long = decimalLongitude) %>%
  # Select core columns for modeling PLUS geographic metadata for analyses
  select(family, gen, sp, authority, collector, number, lat, long, 
         countryCode, year, coordinateUncertaintyInMeters, coordinatePrecision,
         taxonRank, ISLAND, COUNTRY)

# Quick data quality check
cat("GBIF data summary after cleaning:\n")
cat("Total records:", nrow(gbif_cleaned), "\n")
cat("Records with genus:", sum(!is.na(gbif_cleaned$gen)), "\n")
cat("Records with species:", sum(!is.na(gbif_cleaned$sp)), "\n") 
cat("Records with authority:", sum(!is.na(gbif_cleaned$authority)), "\n")
cat("Records with coordinates:", sum(!is.na(gbif_cleaned$lat) & !is.na(gbif_cleaned$long)), "\n")

# Import and integrate Kew Asia Team data ----
kew_data <- read_csv("AsiaTeam_georefs_forSHAWN.csv") %>%
  rename(family = FAMILY,
         gen = GENUS,
         sp = SP1,
         collector = COLLECTOR,
         number = NUMBER,
         lat = LAT,
         long = LONG) %>%
  # Add columns to match GBIF data structure
  mutate(
    authority = NA_character_,           # Kew data likely lacks authority info
    countryCode = NA_character_,         # Will be determined spatially
    year = NA_real_,                     # May not have collection year
    coordinateUncertaintyInMeters = NA_real_,  # Precision info not available
    coordinatePrecision = NA_real_,      # Precision info not available  
    taxonRank = "SPECIES",               # Assume species-level (can verify later)
    ISLAND = NA_character_,              # Will be determined spatially
    COUNTRY = NA_character_              # Will be determined spatially
  ) %>%
  # Reorder columns to match GBIF data
  select(family, gen, sp, authority, collector, number, lat, long, 
         countryCode, year, coordinateUncertaintyInMeters, coordinatePrecision,
         taxonRank, ISLAND, COUNTRY)

# Combine GBIF and Kew datasets
combined_occurrences <- bind_rows(
  gbif_cleaned %>% mutate(source = "GBIF"),
  kew_data %>% mutate(source = "Kew")
)

# Summary statistics ----
cat("Dataset summary:\n")
cat("Total GBIF records:", nrow(gbif_cleaned), "\n")
cat("Total Kew records:", nrow(kew_data), "\n") 
cat("Combined records:", nrow(combined_occurrences), "\n")
cat("Unique species:", combined_occurrences %>% 
      filter(!is.na(gen), !is.na(sp)) %>% 
      distinct(gen, sp) %>% 
      nrow(), "\n")

# Geographic distribution summary
cat("\nGeographic distribution:\n")
combined_occurrences %>% 
  filter(!is.na(ISLAND)) %>% 
  count(ISLAND, COUNTRY, sort = TRUE) %>%
  print()

# Preview the data structure for rWCVP workflow
cat("\nData structure preview for rWCVP name resolution:\n")
combined_occurrences %>% 
  filter(!is.na(gen), !is.na(sp)) %>%
  select(family, gen, sp, authority, ISLAND, COUNTRY) %>%
  slice_head(n = 5) %>%
  print()

# Generate date stamp and save cleaned dataset ----
date_stamp <- format(Sys.Date(), "%Y%m%d")
output_filename <- paste0("vascular_plants_occurrences_new_guinea_borneo_", date_stamp, ".csv")

write_csv(combined_occurrences, output_filename)
cat("Combined dataset saved as:", output_filename, "\n")
