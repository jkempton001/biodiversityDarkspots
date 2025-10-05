# One-time requirements

# ── Libraries (your list as-is) ─────────────────────────────────────────────────
library(terra); library(maps); library(mapdata); library(letsR); library(sf)
library(dplyr); library(ggplot2); library(scico); library(rnaturalearth)
library(purrr); library(smoothr); library(readr); library(tidyverse)
library(CoordinateCleaner); library(countrycode); library(ggrastr); library(rgbif)
library(magrittr); library(bit64); library(keyring); library(geodata)
library(stringr); library(tidyr); library(units)

# Basis-of-record

VALID_BOR <- c(
  "OBSERVATION","MACHINE_OBSERVATION","HUMAN_OBSERVATION",
  "MATERIAL_SAMPLE","MATERIAL_CITATION","PRESERVED_SPECIMEN","OCCURRENCE"
)

bor_tag <- function(bor_vec) {
  if (is.null(bor_vec)) return("BOR-all")                # not used in this proposal, but safe
  bor_vec <- toupper(bor_vec)
  if (length(bor_vec) == 1) paste0("BOR-", bor_vec) else "BOR-multi"
}


# ── Credentials ────────────────────────────────────────────────────────────────
get_gbif_credentials <- function(service_user = "jkempton001") {
  ensure <- function(service_name){
    tryCatch(key_get(service_name, username = service_user),
             error = function(e) { message("Enter ", service_name); key_set(service_name, username = service_user); 
               key_get(service_name, username = service_user)})
  }
  list(user = ensure("gbif_user"),
       password = ensure("gbif_password"),
       email = ensure("gbif_email"))
}
creds <- get_gbif_credentials()
GBIF_USER  <- creds$user; GBIF_PWD <- creds$password; GBIF_EMAIL <- creds$email
cat("GBIF credentials loaded from keychain\n")

# ── Region definitions (GADM filters for GBIF) ────────────────────────────────
REGIONS <- list(
  papua_barat = "IDN.22_1",
  papua       = "IDN.23_1",
  kalimantan  = c("IDN.12_1","IDN.13_1","IDN.14_1","IDN.34_1","IDN.35_1"),
  sumatra     = c("IDN.1_1","IDN.30_1","IDN.31_1","IDN.32_1","IDN.16_1","IDN.24_1","IDN.8_1","IDN.5_1","IDN.17_1","IDN.3_1"),
  java        = c("IDN.9_1","IDN.10_1","IDN.11_1","IDN.7_1","IDN.33_1","IDN.4_1"),
  bali        = "IDN.2_1",
  sulawesi    = c("IDN.25_1","IDN.26_1","IDN.27_1","IDN.28_1","IDN.29_1","IDN.6_1"),
  maluku      = c("IDN.18_1","IDN.19_1"),
  lesser_sunda_islands = c("IDN.20_1","IDN.21_1","TLS"),
  png         = "PNG",
  phillipines = "PHL",
  brunei      = "BRN",
  sabah       = "MYS.13_1",
  sarawak     = "MYS.14_1"
)

# For tagging records with island & country labels in your outputs
REGION_META <- list(
  papua       = list(island="New Guinea", country="Indonesian Papua"),
  papua_barat = list(island="New Guinea", country="Indonesian Papua"),
  kalimantan  = list(island="Borneo", country="Indonesian Borneo"),
  brunei      = list(island="Borneo", country="Brunei"),
  phillipines = list(island="Phillipines", country="Phillipines"),
  sabah       = list(island="Borneo", country="Malaysian Borneo"),
  sarawak     = list(island="Borneo", country="Malaysian Borneo"),
  java        = list(island="Java", country="Java"),
  sumatra     = list(island="Sumatra", country="Sumatra"),
  bali        = list(island="Bali", country="Bali"),
  maluku      = list(island="Maluku", country="Maluku"),
  sulawesi    = list(island="Sulawesi", country="Sulawesi"),
  lesser_sunda_islands = list(island="Lesser Sunda Islands", country="Lesser Sunda Islands"),
  png         = list(island="New Guinea", country="Papua New Guinea")
)

# ── Plot region (same window you used) ────────────────────────────────────────

make_region_poly <- function(
    region_geom = NULL,         # <- union of polygons (sf) to define the map window
    region_bbox = NULL,         # <- alternatively, pass an explicit bbox (named vector)
    map_margin_deg = c(0.5, 0.5, 0.5, 0.5)  # c(left, right, bottom, top) padding (deg)
) {
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  malaysia <- rnaturalearth::ne_countries(country = "Malaysia", returnclass = "sf")
  east_my <- sf::st_crop(malaysia, c(xmin = 108, xmax = 131, ymin = -7, ymax = 8))
  
  # Fallback: original “Indo-Pacific” selection if no geom/bbox was provided
  default_keep <- world |>
    dplyr::filter(name %in% c("Indonesia","Brunei","Philippines","Papua New Guinea","Timor-Leste")) |>
    dplyr::bind_rows(east_my) |>
    sf::st_make_valid() |>
    sf::st_union()
  
  # Determine the region (land) polygon we’ll display
  base_region <- sf::st_make_valid(region_geom %||% default_keep)
  
  # Determine bbox: either from user, or from the region_geom/default_keep
  if (is.null(region_bbox)) {
    bb <- sf::st_bbox(base_region)
  } else {
    stopifnot(all(c("xmin","xmax","ymin","ymax") %in% names(region_bbox)))
    bb <- region_bbox
  }
  
  # Apply margin (degrees)
  if (length(map_margin_deg) == 1) map_margin_deg <- rep(map_margin_deg, 4)
  bb["xmin"] <- bb["xmin"] - map_margin_deg[1]
  bb["xmax"] <- bb["xmax"] + map_margin_deg[2]
  bb["ymin"] <- bb["ymin"] - map_margin_deg[3]
  bb["ymax"] <- bb["ymax"] + map_margin_deg[4]
  
  # Clip region to that bbox
  bbox_sfc <- sf::st_as_sfc(bb)
  regionPoly <- suppressWarnings(sf::st_intersection(base_region, bbox_sfc)) |>
    sf::st_make_valid() |>
    sf::st_as_sf()
  
  # Boundary lines, clipped and filtered to the plotting area
  borders <- rnaturalearth::ne_download(
    scale = "large", type = "admin_0_boundary_lines_land",
    category = "cultural", returnclass = "sf"
  ) |> sf::st_make_valid()
  
  borders_clip <- suppressWarnings(sf::st_intersection(borders, bbox_sfc))
  borders_clip <- borders_clip[sf::st_intersects(borders_clip, regionPoly, sparse = FALSE), ]
  
  list(regionPoly = regionPoly, borders_clip = borders_clip, bbox = bb)
}
# ---------- NEW: utils to construct a region geometry from REGIONS codes -----

# Safe null-coalescing helper
`%||%` <- function(x, y) if (is.null(x)) y else x

# Pull the GADM/ISO3 codes for a set of scope names (REGIONS entries)
codes_for_scope <- function(scope_names) {
  stopifnot(all(scope_names %in% names(REGIONS)))
  unlist(REGIONS[scope_names], use.names = FALSE)
}

# Given a single code, return an sf polygon
# - "IDN.22_1" style => GADM level 1 feature
# - "PNG" style (3-letter ISO) => admin0 country polygon
geom_from_code <- function(code, gadm_cache = new.env(parent = emptyenv())) {
  # Returns an sf polygon feature in EPSG:4326 for either:
  # - ISO3 country code (e.g., "PNG"), or
  # - GADM L1 code (e.g., "IDN.22_1")
  
  out <- if (grepl("^[A-Z]{3}$", code)) {
    # ---- ISO3 COUNTRY BRANCH ----
    w <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
    code <- toupper(code)
    
    # pick the first ISO3-like column that exists
    iso_cols <- intersect(c("iso_a3", "adm0_a3", "gu_a3"), names(w))
    if (length(iso_cols) == 0) stop("Could not find any ISO3 column in Natural Earth data.")
    iso_col <- iso_cols[1]
    
    cand <- w[toupper(w[[iso_col]]) == code, , drop = FALSE]
    if (nrow(cand) == 0) stop("Could not find country with ISO3=", code)
    cand
    
  } else if (grepl("^[A-Z]{3}\\.\\d+_\\d+$", code)) {
    # ---- GADM CODE BRANCH (e.g., "IDN.22_1") ----
    iso3  <- sub("\\..*$", "", code)
    level <- as.integer(sub(".*_(\\d+)$", "\\1", code))
    key   <- paste0(iso3, "_L", level)
    
    gadm_obj <- if (!exists(key, envir = gadm_cache)) {
      x <- geodata::gadm(country = iso3, level = level, path = tempdir())
      assign(key, x, envir = gadm_cache); x
    } else {
      get(key, envir = gadm_cache)
    }
    
    # coerce SpatVector -> sf if needed
    if (inherits(gadm_obj, "SpatVector")) gadm_obj <- sf::st_as_sf(gadm_obj)
    
    gid_col <- paste0("GID_", level)
    if (!gid_col %in% names(gadm_obj))
      stop("GADM table for ", iso3, " level ", level, " lacks column ", gid_col)
    
    sel <- gadm_obj[gadm_obj[[gid_col]] == code, ]
    if (nrow(sel) == 0) stop("Could not find GADM feature for code: ", code)
    sel
    
  } else {
    stop("Unrecognised region code format: ", code)
  }
  
  # ---- Normalise to valid EPSG:4326 sf ----
  out <- sf::st_make_valid(out)
  if (is.na(sf::st_crs(out))) sf::st_crs(out) <- sf::st_crs(4326)
  if (!identical(sf::st_crs(out)$epsg, 4326L)) out <- sf::st_transform(out, 4326)
  
  out
}



# Union a set of codes (IDN.22_1, PNG, etc.) into one geometry
# ---- robust: union a set of codes into one MULTIPOLYGON sf in EPSG:4326 ----
region_geom_from_codes <- function(codes) {
  stopifnot(length(codes) > 0)
  
  # standardise one sf to 2D MULTIPOLYGON in EPSG:4326 and return geometry (sfc)
  standardise <- function(sfobj) {
    obj <- sf::st_make_valid(sfobj)
    # enforce CRS
    if (is.na(sf::st_crs(obj))) sf::st_crs(obj) <- sf::st_crs(4326)
    if (!identical(sf::st_crs(obj)$epsg, 4326L)) obj <- sf::st_transform(obj, 4326)
    
    # drop non-area bits & force 2D MULTIPOLYGON
    geom <- sf::st_geometry(obj)
    geom <- sf::st_zm(geom, drop = TRUE, what = "ZM")                       # drop Z/M
    geom <- sf::st_collection_extract(geom, "POLYGON", warn = FALSE)        # only polygons
    geom <- sf::st_cast(geom, "MULTIPOLYGON", warn = FALSE)                 # unify types
    geom
  }
  
  # get each code’s geometry (sfc)
  parts <- lapply(codes, function(cd) standardise(geom_from_code(cd)))
  
  # guard against empties
  lens <- vapply(parts, length, integer(1))
  if (any(lens == 0))
    stop("At least one code resolved to an empty geometry: ",
         paste(codes[which(lens==0)], collapse = ", "))
  
  # concatenate sfc safely and union
  combined <- do.call(c, parts)                      # c.sfc
  sf::st_crs(combined) <- sf::st_crs(4326)           # ensure CRS sticks
  unioned  <- sf::st_union(sf::st_combine(combined)) # dissolve internal borders
  
  # wrap back into sf
  sf::st_as_sf(data.frame(id = 1), geometry = sf::st_sfc(unioned, crs = 4326))
}



#=============================== CORE BUILDING BLOCKS =========================#

# 1) Download: by regions for a taxonKey
occ_download_by_regions <- function(
    taxon_key,
    scope_names,
    format = "DWCA",
    basis_of_record = VALID_BOR,     # if you added BOR param earlier
    max_concurrent = 3,
    poll_every_sec = 120,
    retry_max = 6
) {
  stopifnot(all(scope_names %in% names(REGIONS)))
  if (!all(toupper(basis_of_record) %in% VALID_BOR)) {
    bad <- setdiff(toupper(basis_of_record), VALID_BOR)
    stop("Unknown basis_of_record value(s): ", paste(bad, collapse = ", "))
  }
  
  # helper to launch one download, with retry if GBIF says "too many simultaneous"
  launch_one <- function(region_name) {
    region_codes <- REGIONS[[region_name]]
    
    preds <- list(
      pred("taxonKey", taxon_key),
      if (length(basis_of_record) == 1) pred("basisOfRecord", toupper(basis_of_record))
      else pred_in("basisOfRecord", toupper(basis_of_record)),
      pred("hasGeospatialIssue", FALSE),
      pred("hasCoordinate", TRUE),
      pred("occurrenceStatus", "PRESENT")
    )
    
    if (length(region_codes) == 1) {
      preds <- append(preds, list(pred("gadm", region_codes)), after = 1)
    } else {
      preds <- append(preds, list(pred_in("gadm", region_codes)), after = 1)
    }
    
    attempt <- 1
    repeat {
      key <- tryCatch({
        do.call(
          occ_download,
          c(preds, format = format, user = GBIF_USER, pwd = GBIF_PWD, email = GBIF_EMAIL)
        )
      }, error = function(e) {
        msg <- tolower(conditionMessage(e))
        if (grepl("download limitation is exceeded|too many simultaneous downloads", msg) &&
            attempt < retry_max) {
          wait <- poll_every_sec * attempt
          message("GBIF limit hit starting '", region_name, "'. Retrying in ", wait, "s (attempt ",
                  attempt + 1, "/", retry_max, ") …")
          Sys.sleep(wait)
          return(NULL)  # signal retry
        }
        stop(e)
      })
      if (!is.null(key)) break
      attempt <- attempt + 1
    }
    message("Started ", region_name, " => ", key)
    key
  }
  
  # status helper
  get_status <- function(key) {
    meta <- rgbif::occ_download_meta(key)
    if (is.null(meta$status)) return(NA_character_)
    meta$status
  }
  
  # queue/active/results
  queue <- scope_names
  active <- list()   # region_name -> key
  results <- list()  # region_name -> key (finished)
  
  cat("Started downloads for:", paste(scope_names, collapse = ", "),
      "\nBasis of record filter:", paste(basis_of_record, collapse = ", "),
      "\nMax concurrent:", max_concurrent, "\n")
  
  while (length(queue) > 0 || length(active) > 0) {
    
    # fill up to concurrency limit
    while (length(queue) > 0 && length(active) < max_concurrent) {
      nm <- queue[[1]]; queue <- queue[-1]
      key <- launch_one(nm)
      active[[nm]] <- key
    }
    
    # poll current actives
    if (length(active) > 0) {
      statuses <- purrr::imap_chr(active, function(key, nm) {
        st <- get_status(key)
        message("  [", nm, "] status: ", st)
        st
      })
      
      # move finished jobs to results (keep FAILED/KILLED/CANCELLED but warn)
      finished <- names(statuses)[statuses %in% c("SUCCEEDED","FAILED","KILLED","CANCELLED")]
      if (length(finished) > 0) {
        for (nm in finished) {
          if (statuses[[nm]] != "SUCCEEDED")
            warning("Download for '", nm, "' ended with status ", statuses[[nm]],
                    ". Key: ", active[[nm]])
          results[[nm]] <- active[[nm]]
          active[[nm]] <- NULL
        }
      }
      
      # if still have running jobs, wait before next poll
      if (length(active) > 0) Sys.sleep(poll_every_sec)
    }
  }
  
  # return keys in the same naming style as before
  results[scope_names]
}

# ---- NEW: robust DWCA reader (handles quoted newlines etc.) -------------------
read_dwca_occurrence <- function(download_key) {
  zip_path <- rgbif::occ_download_get(download_key, overwrite = TRUE)

  # unzip to a temp dir
  tmpdir <- tempfile("dwca_"); dir.create(tmpdir)
  utils::unzip(zip_path, exdir = tmpdir)

  # find the occurrence core file
  occ_file <- list.files(
    tmpdir,
    pattern = "^occurrence\\.(txt|csv)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(occ_file) == 0) stop("occurrence core not found in the archive")

  # robust read (treat everything as character first; respect quotes/newlines)
  df <- readr::read_delim(
    occ_file[1],
    delim = "\t",
    col_types = readr::cols(.default = readr::col_character()),
    na = c("", "NA"),
    quote = "\"",
    progress = FALSE,
    locale = readr::locale(encoding = "UTF-8")
  )
  message("Read ", nrow(df), " rows from occurrence core.")
  df
}


# 2) Import + standardize one archive
import_one_dwca <- function(download_key, island, country) {

  
  # --- NEW: use robust DWCA reader -------------------------------------------
  dat <- read_dwca_occurrence(download_key)
  
  dat |>
    dplyr::select(dplyr::any_of(c(
      "gbifID", "license", "institutionCode", "collectionCode", "occurrenceID",
      "catalogNumber", "recordNumber", "recordedBy", "eventDate", "habitat",
      "eventRemarks", "locality", "decimalLatitude", "decimalLongitude",
      "coordinateUncertaintyInMeters", "coordinatePrecision", "identifiedBy",
      "scientificName", "kingdom", "phylum", "class", "order", "family", "genus",
      "taxonRank", "taxonomicStatus", "elevation", "elevationAccuracy", "issue",
      "taxonKey", "acceptedTaxonKey", "speciesKey", "species", "acceptedScientificName",
      "verbatimScientificName", "iucnRedListCategory", "countryCode",
      "individualCount", "year", "basisOfRecord", "datasetName",
      "establishmentMeans", "references"
    ))) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) |>
    dplyr::mutate(
      decimalLatitude = suppressWarnings(as.numeric(decimalLatitude)),
      decimalLongitude = suppressWarnings(as.numeric(decimalLongitude)),
      coordinateUncertaintyInMeters = suppressWarnings(as.numeric(coordinateUncertaintyInMeters)),
      coordinatePrecision = suppressWarnings(as.numeric(coordinatePrecision)),
      elevation = suppressWarnings(as.numeric(elevation)),
      elevationAccuracy = suppressWarnings(as.numeric(elevationAccuracy)),
      individualCount = suppressWarnings(as.numeric(individualCount)),
      year = suppressWarnings(as.numeric(year))
    ) |>
    tibble::add_column(ISLAND = island, COUNTRY = country)
}

# 3) Import all regions with tagging
import_and_bind <- function(download_keys_named) {
  out <- map(names(download_keys_named), function(nm){
    mm <- REGION_META[[nm]]
    tryCatch(import_one_dwca(download_keys_named[[nm]], mm$island, mm$country),
             error = function(e){ message("Error importing ", nm, ": ", e$message); NULL })
  }) |> compact()
  bind_rows(out)
}

# 4) Light taxonomic and column standardization
standardize_for_modeling <- function(df){
  df |>
    separate(acceptedScientificName, into = c("gen", "sp", "authority"), sep = " ",
             extra = "merge", fill = "right") |>
    rename(collector = recordedBy, number = recordNumber, lat = decimalLatitude, long = decimalLongitude) |>
    select(family, gen, sp, authority, collector, number, lat, long,
           countryCode, year, coordinateUncertaintyInMeters, coordinatePrecision,
           taxonRank, ISLAND, COUNTRY)
}

# 5) Cleaning (CoordinateCleaner + clip to study shapefile + species-level filter)
clean_occurrences <- function(df, study_shp_path) {
  # Basic geo filter to your study countries
  df1 <- df |>
    filter(!is.na(lat), !is.na(long)) |>
    filter(countryCode %in% c("ID","PG","MY","PH","BN","TL"))
  
  df2 <- df1 |>
    mutate(species = if_else(!is.na(gen) & !is.na(sp), paste(gen, sp), NA_character_)) |>
    cc_dupl(lon="long", lat="lat") |>
    cc_zero(lon="long", lat="lat") |>
    cc_equ(lon="long", lat="lat") |>
    cc_val(lon="long", lat="lat") |>
    cc_cap(lon="long", lat="lat", buffer = 2000) |>
    cc_cen(lon="long", lat="lat", buffer = 2000) |>
    cc_gbif(lon="long", lat="lat", buffer = 2000) |>
    cc_inst(lon="long", lat="lat", buffer = 2000) |>
    select(-species)
  
  # --- NEW: read study area passed in as argument ---
  study_area <- sf::st_read(study_shp_path, quiet = TRUE)
  
  # Convert to sf points
  pts <- df2 |>
    dplyr::filter(!is.na(lat), !is.na(long)) |>
    sf::st_as_sf(coords = c("long", "lat"), crs = 4326)
  
  if (nrow(pts) == 0) {
    warning("No points after CoordinateCleaner. Returning empty tibble.")
    return(dplyr::slice(df2, 0))
  }
  
  # Clip to study area
  inside <- sf::st_intersection(pts, study_area)
  
  df2_5 <- inside %>%
    dplyr::mutate(
      long = sf::st_coordinates(.)[, 1],
      lat  = sf::st_coordinates(.)[, 2]
    ) %>%
    sf::st_drop_geometry()
  
  # Species-level & certain IDs
  df5 <- df2_5 |>
    filter(!is.na(gen), !is.na(sp)) |>
    filter(!grepl("\\bsp\\.?$|\\bcf\\.?\\b|\\baff\\.?\\b", sp, ignore.case = TRUE))
  
  df5
}

# 6) Find cleaned data, if it already exists
find_latest_cleaned_csv <- function(taxon_label, suffix = NULL) {
  base <- paste0("^gbifSEAsia_", stringr::str_to_title(taxon_label), "OccDataCleaned_\\d{8}")
  if (!is.null(suffix)) base <- paste0(base, "_", stringr::str_replace_all(suffix, "[^A-Za-z0-9_-]", ""))
  pat <- paste0(base, "\\.csv$")
  files <- list.files(pattern = pat)
  if (length(files) == 0) return(NULL)
  dates <- as.Date(stringr::str_extract(files, "\\d{8}"), "%Y%m%d")
  files[order(dates, decreasing = TRUE)][1]
}


# 7) Grid heatmap (same styling; returns ggplot object and writes SVG)
make_heatmap <- function(clean_df, out_prefix, region_pack = make_region_poly(), cellsize = 0.5) {
  regionPoly   <- region_pack$regionPoly
  borders_clip <- region_pack$borders_clip
  
  xy <- clean_df |> select(gen, sp, long, lat) |> tidyr::drop_na() |>
    mutate(id = row_number())
  
  bb <- sf::st_bbox(regionPoly)   # <- get bounding box from sf object
  r <- terra::rast(
    xmin = bb["xmin"], xmax = bb["xmax"],
    ymin = bb["ymin"], ymax = bb["ymax"],
    resolution = cellsize,
    crs = "EPSG:4326"
  )
  
  pts <- terra::vect(xy[, c("long","lat")], geom = c("long","lat"), crs = "EPSG:4326")
  rcount <- terra::rasterize(pts, r, fun = "count", background = 0)
  rcount <- terra::mask(rcount, terra::vect(regionPoly))
  
  # After you compute `rcount` and mask it
  names(rcount) <- "numCollections"  # give the layer a stable name
  
  # Convert to polygons AND keep the values as an attribute
  grid_v  <- terra::as.polygons(rcount, dissolve = FALSE, values = TRUE)
  
  # Drop NA cells (outside mask) to keep plot tidy
  grid_v  <- grid_v[!is.na(terra::values(grid_v)$numCollections), ]
  
  # To sf
  grid_sf <- sf::st_as_sf(grid_v)
  
  # Make sure type is integer
  grid_sf$numCollections <- as.integer(grid_sf$numCollections)
  
  
  # --- Build a one-time sea mask = bbox minus land (used to hide overhang) ---
  bb_sfc   <- sf::st_as_sfc(sf::st_bbox(regionPoly))
  land     <- sf::st_union(regionPoly)
  sea_mask <- suppressWarnings(sf::st_difference(bb_sfc, land))
  
  # Ocean & graticule styling
  ocean_fill <- "grey90"
  
  # 
  # --- Bin + palette (unchanged) ---
  bin_levels <- c("0","1–25","26–100","101–500","501–1000","1001–1500","1501–2000",">2000")

  grid_sf_binned <- grid_sf %>%
    dplyr::mutate(
      collection_bin = cut(
        numCollections,
        breaks = c(-Inf, 0, 25, 100, 500, 1000, 1500, 2000, Inf),
        labels = bin_levels,
        include.lowest = TRUE, right = TRUE
      ),
      collection_bin = factor(collection_bin, levels = bin_levels)
    )

  pal <- c(
    "0"          = "#FFFFFF",
    "1–25"       = "#EAF2FF",
    "26–100"     = "#D7E9FF",
    "101–500"    = "#C3DEFF",
    "501–1000"   = "#AECFFF",
    "1001–1500"  = "#FFF6BF",
    "1501–2000"  = "#FFE1B8",
    ">2000"      = "#FFB6C8"
  )


  # --- ALIGN GRATICULE WITH TICKS & PANEL GRID ----
  bb <- sf::st_bbox(regionPoly)
  
  # choose steps
  lon_step <- 10
  lat_step <- 5
  
  # snap to neat multiples so ticks fall on round degrees
  lon_min <- floor(bb["xmin"] / lon_step) * lon_step
  lon_max <- ceiling(bb["xmax"] / lon_step) * lon_step
  lat_min <- floor(bb["ymin"] / lat_step) * lat_step
  lat_max <- ceiling(bb["ymax"] / lat_step) * lat_step
  
  lon_breaks <- seq(lon_min, lon_max, by = lon_step)
  lat_breaks <- seq(lat_min, lat_max, by = lat_step)
  
  # build graticule exactly at those breaks
  graticule <- sf::st_graticule(
    x   = bb,
    crs = sf::st_crs(regionPoly),
    lon = lon_breaks,
    lat = lat_breaks
  )
  
  # keep only the parts of the graticule that lie over the sea
  grat_sea <- suppressWarnings(sf::st_intersection(graticule, sea_mask))
  
  # pretty axis labels (optional)
  lab_lon <- function(x) paste0(abs(x), "°", ifelse(x < 0, "W", ifelse(x > 0, "E", "")))
  lab_lat <- function(y) paste0(abs(y), "°", ifelse(y < 0, "S", ifelse(y > 0, "N", "")))
  
  # --- PLOT (layer order matters) ---
  p <- ggplot() +
    # heatmap cells
    geom_sf(data = grid_sf_binned, aes(fill = collection_bin), color = NA) +
    # sea mask (hides coastal overhang)
    geom_sf(data = sea_mask, fill = ocean_fill, color = NA) +
    # graticule *on top of sea* so lines are visible and do not cross land
    geom_sf(data = grat_sea, color = "white", linewidth = 0.5) +
    # color scale
    scale_fill_manual(
      name   = "Sample locations",
      values = pal,
      limits = bin_levels,
      drop   = FALSE,
      guide_legend(override.aes = list(fill = pal))
    ) +
    # make ticks align with our graticule
    scale_x_continuous(breaks = lon_breaks, labels = lab_lon, expand = c(0, 0)) +
    scale_y_continuous(breaks = lat_breaks, labels = lab_lat, expand = c(0, 0)) +
    # keep the plotted window exactly to the bbox so breaks land cleanly
    coord_sf(
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"]),
      expand = FALSE, clip = "on"
    ) +
    # outlines
    geom_sf(data = regionPoly, fill = NA, colour = "grey60", linewidth = 0.2) +
    geom_sf(data = borders_clip, color = "grey60", linewidth = 0.2) +
    theme(
      # use our own graticule lines; turn off the background grid
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = ocean_fill, colour = NA),
      legend.position        = "inside",
      legend.position.inside = c(0.96, 0.96),
      legend.justification   = c(1, 1),
      legend.background      = element_rect(fill = NA, colour = NA),
      legend.title           = element_text(size = 9, face = 'bold'),
      legend.margin          = margin(4,4,4,4),
      legend.key.width       = unit(0.4, "cm"),
      legend.key.height      = unit(0.4, "cm")
    )
  
    
  ggsave(
    filename = paste0(out_prefix, "_heatmap.svg"),
    plot = p, width = 8, height = 6, units = "in", dpi = 300
  )
}

# 7) Region density table (per COUNTRY by default; change to ISLAND if preferred)
region_density <- function(clean_df, region_col = "COUNTRY") {
  counts <- clean_df |>
    count(.data[[region_col]], name = "Observations") |>
    rename(Region = !!region_col)
  
  # Your areas table (keep as-is; adjust if you change region_col)
  areas <- data.frame(
    Region = c("Bali","Brunei","Indonesian Borneo","Indonesian Papua",
               "Java","Lesser Sunda Islands","Malaysian Borneo","Maluku",
               "Papua New Guinea","Phillipines","Sulawesi","Sumatra"),
    Area_km2 = c(5590.15,5765,534698.27,412214.61,132598.77,67128.38 + 14950,
                 198445.64,78897,462840,300000,174416.16,482286.55)
  )
  
  counts |>
    left_join(areas, by = "Region") |>
    mutate(Density = Observations / Area_km2)
}

# =============================== MASTER RUNNER ===============================#

# Master runner:
# - If you already have GBIF download keys, pass them in as named list to skip downloading.
# - scope = character vector of names from REGIONS (e.g., c("papua","papua_barat","png",...))
# - study_shp_path = "territory_selection.shp" (as per your code)

run_gbif_pipeline <- function(taxon_label,
                              taxon_key,
                              scope,
                              study_shp_path = "territory_selection.shp",
                              download_keys_override = NULL,
                              reuse_cleaned = TRUE,
                              cleaned_csv = NULL,
                              basis_of_record = VALID_BOR,
                              run_suffix = NULL,
                              max_concurrent = 3,
                              poll_every_sec = 15,
                              # NEW: mapping extent controls
                              map_scope = NULL,           # vector of names from REGIONS to define map extent
                              map_bbox = NULL,            # explicit bbox override (named c(xmin,xmax,ymin,ymax))
                              map_margin_deg = 0.5)     # padding for the bbox
  {
  date_stamp <- format(Sys.Date(), "%Y%m%d")
  ttl <- stringr::str_to_title(taxon_label)
  borT <- bor_tag(basis_of_record)
  suffix <- if (!is.null(run_suffix)) run_suffix else borT
  
  message("=== ", toupper(taxon_label), " / taxonKey=", taxon_key,
          " / ", paste0("BOR=", paste(basis_of_record, collapse=",")), " ===")
  
  # 0) Reuse cleaned CSV?
  if (!is.null(cleaned_csv) && file.exists(cleaned_csv)) {
    message("Reusing cleaned occurrences from: ", cleaned_csv)
    clean <- readr::read_csv(cleaned_csv, show_col_types = FALSE)
    std <- NULL; dl_keys <- NULL
  } else if (isTRUE(reuse_cleaned)) {
    latest <- find_latest_cleaned_csv(taxon_label, suffix = suffix)  # <-- updated
    if (!is.null(latest)) {
      message("Reusing latest cleaned CSV found in working directory: ", latest)
      clean <- readr::read_csv(latest, show_col_types = FALSE)
      std <- NULL; dl_keys <- NULL
    } else {
      clean <- NULL
    }
  } else {
    clean <- NULL
  }
  
  if (!is.null(clean)) {
    rp <- make_region_poly()
    plot_obj <- make_heatmap(clean,
                             out_prefix = paste0(taxon_label, "_digitised_occurrences_", suffix),
                             region_pack = rp)
    dens <- region_density(clean, region_col = "COUNTRY")
    dens_out <- paste0(taxon_label, "_occurrence_density_", date_stamp, "_", suffix, ".csv")
    readr::write_csv(dens, dens_out)
    message("Saved density table: ", dens_out)
    
    return(list(
      download_keys = dl_keys,
      standardized  = NULL,
      cleaned       = clean,
      heatmap       = plot_obj,
      density_table = dens
    ))
  }
  
  # 1) Download
  if (is.null(download_keys_override)) {
    dl_keys <- occ_download_by_regions(
      taxon_key,
      scope,
      basis_of_record = basis_of_record,
      max_concurrent = max_concurrent,
      poll_every_sec = poll_every_sec
    )
  } else {
    dl_keys <- download_keys_override
    stopifnot(all(names(dl_keys) %in% names(REGIONS)))
    message("Using provided download keys (skipping new downloads).")
  }
  
  # 2) Import & bind
  raw <- import_and_bind(dl_keys)
  
  # 3) Standardize
  std <- standardize_for_modeling(raw)
  std_out <- paste0("gbifSEAsia_", ttl, "OccDataStandardized_", date_stamp, "_", suffix, ".csv")
  readr::write_csv(std, std_out)
  message("Saved standardized table: ", std_out)
  
  # 4) Clean
  clean <- clean_occurrences(std, study_shp_path)
  clean_out <- paste0("gbifSEAsia_", ttl, "OccDataCleaned_", date_stamp, "_", suffix, ".csv")
  readr::write_csv(clean, clean_out)
  message("Saved cleaned table: ", clean_out)
  
  # 5) Plot
  map_scope_eff <- map_scope %||% scope
  map_codes <- codes_for_scope(map_scope_eff)
  message("Map codes: ", paste(map_codes, collapse = ", "))
  
  # Build unioned map geometry (will error loudly if something is wrong)
  reg_geom <- region_geom_from_codes(map_codes)
  
  # Quick sanity: print bbox we’re about to use
  message("Union bbox (pre-margin): ",
          paste(round(unname(sf::st_bbox(reg_geom)), 4), collapse = ", "))
  
  rp <- make_region_poly(
    region_geom    = reg_geom,
    region_bbox    = map_bbox,        # normally NULL so bbox derives from reg_geom
    map_margin_deg = map_margin_deg
  )
  
  message("Plot bbox (post-margin): ",
          paste(round(unname(rp$bbox), 4), collapse = ", "))
  
  plot_obj <- make_heatmap(clean,
                           out_prefix = paste0(taxon_label, "_digitised_occurrences_", suffix),
                           region_pack = rp)
  
  # 6) Density
  dens <- region_density(clean, region_col = "COUNTRY")
  dens_out <- paste0(taxon_label, "_occurrence_density_", date_stamp, "_", suffix, ".csv")
  readr::write_csv(dens, dens_out)
  message("Saved density table: ", dens_out)
  
  list(
    download_keys = dl_keys,
    standardized  = std,
    cleaned       = clean,
    heatmap       = plot_obj,
    density_table = dens
  )
}


# ============================== TAXON AND GEOGRAPHY ARGUMENTS ================#

mammal_keys <- list(
  papua_barat = "0024029-250920141307145",
  papua       = "0024027-250920141307145",
  png         = "0024057-250920141307145",
  kalimantan  = "0024055-250920141307145",
  phillipines = "0024030-250920141307145",
  sabah       = "0024038-250920141307145",
  sarawak     = "0024039-250920141307145",
  brunei      = "0024037-250920141307145",
  java        = "0027224-250920141307145",
  sumatra     = "0027225-250920141307145",
  sulawesi    = "0027229-250920141307145",
  lesser_sunda_islands = "0027230-250920141307145",
  bali        = "0027231-250920141307145",
  maluku      = "0027244-250920141307145"
)

frog_keys <- list(
  papua_barat = "0020372-250920141307145",
  papua = "0020368-250920141307145", 
  png = "0022860-250920141307145",
  kalimantan = "0022979-250920141307145",
  phillipines = "0020374-250920141307145",
  sabah = "0020403-250920141307145",
  sarawak = "0020405-250920141307145",
  brunei = "0022904-250920141307145",
  java = "0026479-250920141307145",
  sumatra = "0026481-250920141307145",
  bali = "0026498-250920141307145",
  maluku = "0026509-250920141307145",
  sulawesi = "0026482-250920141307145",
  lesser_sunda_islands = "0026497-250920141307145"
)

bird_keys <- list(
  papua_barat = "0024334-250920141307145",
  papua = "0024333-250920141307145", 
  png = "0024457-250920141307145",
  kalimantan = "0024396-250920141307145",
  phillipines = "0024335-250920141307145",
  sabah = "0024394-250920141307145",
  sarawak = "0024395-250920141307145",
  brunei = "0024387-250920141307145",
  java = "0028378-250920141307145",
  sumatra = "0028379-250920141307145",
  sulawesi = "0028380-250920141307145",
  lesser_sunda_islands = "0028510-250920141307145",
  bali = "0028511-250920141307145",
  maluku = "0028564-250920141307145"
)

squamate_keys <- list(
  papua_barat = "0030972-250920141307145",
  papua = "0030971-250920141307145", 
  png = "0031079-250920141307145",
  kalimantan = "0031078-250920141307145",
  phillipines = "0030973-250920141307145",
  sabah = "0030995-250920141307145",
  sarawak = "0030997-250920141307145",
  brunei = "0030992-250920141307145",
  java = "0031084-250920141307145",
  sumatra = "0031122-250920141307145",
  sulawesi = "0031632-250920141307145",
  lesser_sunda_islands = "0031183-250920141307145",
  bali = "0031184-250920141307145",
  maluku = "0031217-250920141307145"
)

plants_keys <- list(
  papua_barat = "0047000-250920141307145",
  papua = "0046996-250920141307145", 
  png = "0047318-250920141307145",
  kalimantan = "0047131-250920141307145",
  phillipines = "0047001-250920141307145",
  sabah = "0047007-250920141307145",
  sarawak = "0047029-250920141307145",
  brunei = "0047005-250920141307145",
  java = "0047132-250920141307145",
  sumatra = "0047133-250920141307145",
  sulawesi = "0047203-250920141307145",
  lesser_sunda_islands = "0047296-250920141307145",
  bali = "0047297-250920141307145",
  maluku = "0047317-250920141307145"
)




frog_res <- run_gbif_pipeline(
  taxon_label = "anura",
  taxon_key   = 952,
  scope       = names(frog_keys),                # data scope (unchanged)
  study_shp_path = "territory_selection.shp",
  download_keys_override = frog_keys,
  reuse_cleaned = FALSE,
  basis_of_record =  c(
    "OBSERVATION","MACHINE_OBSERVATION","HUMAN_OBSERVATION",
    "MATERIAL_SAMPLE","MATERIAL_CITATION","PRESERVED_SPECIMEN","OCCURRENCE"
  ),
  map_scope    = c("papua_barat","papua","png"), 
  map_margin_deg = 0.25
)



mammal_res <- run_gbif_pipeline(
  taxon_label = "mammalia",
  taxon_key   = 359,
  scope       = names(mammal_keys),
  study_shp_path = "territory_selection.shp",
  download_keys_override = mammal_keys,
  reuse_cleaned = FALSE        # default; can even omit
)

bird_res <- run_gbif_pipeline(
  taxon_label = "aves",
  taxon_key   = 212,
  scope       = names(bird_keys),
  study_shp_path = "territory_selection.shp",
  download_keys_override = bird_keys,
  reuse_cleaned = FALSE        # default; can even omit
)

squamates_res <- run_gbif_pipeline(
  taxon_label = "squamates",
  taxon_key   = 11592253,
  scope       = names(squamate_keys),
  study_shp_path = "territory_selection.shp",
  download_keys_override = squamate_keys,
  reuse_cleaned = FALSE        # default; can even omit
)


plants_res <- run_gbif_pipeline(
  taxon_label = "tracheophyta",
  taxon_key   = 7707728,
  scope       = names(plants_keys),
  study_shp_path = "territory_selection.shp",
  download_keys_override = NULL,
  reuse_cleaned = FALSE,
  basis_of_record = "PRESERVED_SPECIMEN"
)
