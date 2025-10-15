# One-time requirements

# ── Libraries (your list as-is) ─────────────────────────────────────────────────
library(terra); library(maps); library(mapdata); library(letsR); library(sf)
library(dplyr); library(ggplot2); library(scico); library(rnaturalearth)
library(purrr); library(smoothr); library(readr); library(tidyverse)
library(CoordinateCleaner); library(countrycode); library(ggrastr); library(rgbif)
library(magrittr); library(bit64); library(keyring); library(geodata)
library(stringr); library(tidyr); library(units)

## ── Per-island map extents ─────────────────────────────────────────────────
ISLAND_BBOX <- list(
  "New Guinea" = c(xmin = 129, xmax = 156,  ymin = -12,  ymax = 4.5),
  "Borneo"     = c(xmin = 108, xmax = 120.5, ymin =  -4.5, ymax = 7.8),
  "Madagascar" = c(xmin =  43, xmax =  51.5, ymin = -26.5, ymax = -11)
)


## ── Island-split writer ────────────────────────────────────────────────────
write_island_svgs <- function(clean_df, out_prefix, study_shp_path) {
  islands <- c("New Guinea", "Borneo", "Madagascar")
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  purrr::walk(islands, function(isle){
    sub <- dplyr::filter(clean_df, .data$ISLAND == isle)
    if (nrow(sub) == 0) {
      message("No points for ", isle, " — skipping.")
      return(invisible(NULL))
    }
    rp <- make_region_poly(bbox = ISLAND_BBOX[[isle]])
    make_heatmap(
      sub,
      out_prefix = paste0(out_prefix, "_", gsub("\\s+","", tolower(isle))), # suffix island name
      region_pack = rp,
      stamp = ts
    )
  })
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
  png         = "PNG",
  brunei      = "BRN",
  sabah       = "MYS.13_1",
  sarawak     = "MYS.14_1",
  madagascar  = "MDG"
)

# For tagging records with island & country labels in your outputs
REGION_META <- list(
  papua       = list(island="New Guinea", country="Indonesian Papua"),
  papua_barat = list(island="New Guinea", country="Indonesian Papua"),
  kalimantan  = list(island="Borneo", country="Indonesian Borneo"),
  brunei      = list(island="Borneo", country="Brunei"),
  sabah       = list(island="Borneo", country="Malaysian Borneo"),
  sarawak     = list(island="Borneo", country="Malaysian Borneo"),
  png         = list(island="New Guinea", country="Papua New Guinea"),
  madagascar  = list(island = "Madagascar", country="Madagascar")
)

# ---- helper: make a concise basisOfRecord tag for filenames ---------------
make_bor_tag <- function(basis_filter) {
  if (is.null(basis_filter) || length(basis_filter) == 0) return("ALLBOR")
  tag <- paste(sort(unique(toupper(basis_filter))), collapse = "+")
  paste0("BOR-", gsub("[^A-Z0-9_+]", "", tag))
}

# ---- helper: make a concise shapefile tag for filenames --------------------
make_shp_tag <- function(shp_path) {
  if (is.null(shp_path) || is.na(shp_path) || !nzchar(shp_path)) return("AOI-unknown")
  nm <- tools::file_path_sans_ext(basename(shp_path))
  nm <- gsub("\\s+", "_", nm)
  nm <- gsub("[^A-Za-z0-9_+\\-]", "", nm)  # <- hyphen placed last in class
  paste0("AOI-", nm)
}





# ── Plot region (same window you used) ────────────────────────────────────────

# Accepts:
# - numeric vector c(xmin, xmax, ymin, ymax), or
# - named list/vector with names xmin/xmax/ymin/ymax, or
# - NULL to use your current default
make_region_poly <- function(bbox = NULL) {
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  malaysia <- rnaturalearth::ne_countries(country = "Malaysia", returnclass = "sf")
  east_my <- st_crop(malaysia, c(xmin=108, xmax=131, ymin=-7, ymax=8))
  keep <- world |> dplyr::filter(name %in% c("Indonesia","Brunei","Philippines","Papua New Guinea","Timor-Leste"))
  region_keep <- dplyr::bind_rows(keep, east_my) |> st_make_valid() |> st_union()
  
  # --- NEW: normalize bbox input ---
  if (is.null(bbox)) {
    bb <- sf::st_bbox(c(xmin = 95, xmax = 156, ymin = -10.3, ymax = 22), crs = sf::st_crs(world))
  } else {
    # support numeric c(xmin, xmax, ymin, ymax) or named
    if (is.numeric(bbox) && length(bbox) == 4 && is.null(names(bbox))) {
      names(bbox) <- c("xmin","xmax","ymin","ymax")
    }
    stopifnot(all(c("xmin","xmax","ymin","ymax") %in% names(bbox)))
    bb <- sf::st_bbox(c(xmin = bbox[["xmin"]], xmax = bbox[["xmax"]],
                        ymin = bbox[["ymin"]], ymax = bbox[["ymax"]]),
                      crs = sf::st_crs(world))
  }
  
  regionPoly <- st_crop(region_keep, bb) |> st_make_valid() |> st_as_sf()
  borders <- rnaturalearth::ne_download(scale="large", type="admin_0_boundary_lines_land",
                                        category="cultural", returnclass="sf") |>
    st_make_valid()
  bbox_sfc <- sf::st_as_sfc(bb)
  borders_clip <- suppressWarnings(sf::st_intersection(borders, bbox_sfc))
  borders_clip <- borders_clip[sf::st_intersects(borders_clip, regionPoly, sparse = FALSE), ]
  list(regionPoly = st_make_valid(regionPoly), borders_clip = borders_clip, bbox = bb)
}

#=============================== CORE BUILDING BLOCKS =========================#

# 1) Download: by regions for a taxonKey
occ_download_by_regions <- function(taxon_key, scope_names, format = "DWCA") {
  stopifnot(all(scope_names %in% names(REGIONS)))
  keys <- map(scope_names, function(nm){
    region_codes <- REGIONS[[nm]]
    preds <- list(
      pred("taxonKey", taxon_key),
      pred_in("basisOfRecord", c("OBSERVATION","MACHINE_OBSERVATION","HUMAN_OBSERVATION",
                                 "MATERIAL_SAMPLE","MATERIAL_CITATION","PRESERVED_SPECIMEN","OCCURRENCE")),
      pred("hasGeospatialIssue", FALSE),
      pred("hasCoordinate", TRUE),
      pred("occurrenceStatus", "PRESENT")
    )
    if (length(region_codes) == 1) {
      preds <- append(preds, list(pred("gadm", region_codes)), after = 1)
    } else {
      preds <- append(preds, list(pred_in("gadm", region_codes)), after = 1)
    }
    do.call(occ_download, c(preds, format = format, user = GBIF_USER, pwd = GBIF_PWD, email = GBIF_EMAIL))
  })
  names(keys) <- scope_names
  cat("Started downloads for:", paste(scope_names, collapse = ", "), "\n")
  keys
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
           taxonRank, ISLAND, COUNTRY, basisOfRecord)   # << keep basisOfRecord
}

# ---- NEW: filter standardized data by basisOfRecord if requested -----------
filter_by_basis <- function(std_df, basis_filter = NULL){
  if (is.null(basis_filter) || length(basis_filter) == 0) return(std_df)
  std_df %>% dplyr::filter(.data$basisOfRecord %in% basis_filter)
}


# 5) Cleaning (CoordinateCleaner + clip to study shapefile + species-level filter)
clean_occurrences <- function(df, study_shp_path) {
  # Basic geo filter to your study countries
  df1 <- df |>
    filter(!is.na(lat), !is.na(long)) |>
    filter(countryCode %in% c("ID","PG","MY","PH","BN","TL","MG"))
  
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
find_latest_cleaned_csv <- function(taxon_label, bor_tag = "ALLBOR", shp_tag = "AOI-unknown") {
  base <- paste0("^gbifSEAsia_", stringr::str_to_title(taxon_label),
                 "OccDataCleaned_", bor_tag, "_", shp_tag, "_")
  # Accept YYYYMMDD or YYYYMMDD_HHMMSS
  pat <- paste0(base, "(\\d{8})(?:_(\\d{6}))?\\.csv$")
  files <- list.files(pattern = pat)
  if (length(files) == 0) return(NULL)
  # Extract date and optional time; order by date then time (missing time treated as 000000)
  m <- stringr::str_match(files, paste0(base, "(\\d{8})(?:_(\\d{6}))?\\.csv$"))
  d <- as.Date(m[,2], "%Y%m%d")
  t <- m[,3]; t[is.na(t)] <- "000000"
  ord <- order(d, t, decreasing = TRUE)
  files[ord][1]
}

# 7) Grid heatmap (same styling; returns ggplot object and writes SVG)
make_heatmap <- function(clean_df, out_prefix, region_pack = make_region_poly(), cellsize = 0.5, stamp = NULL) {
  # Unpack region data prepared by make_region_poly()
  regionPoly   <- region_pack$regionPoly
  borders_clip <- region_pack$borders_clip
  
  # ✅ Always honor the requested bbox if present; otherwise fall back to land bbox
  bb <- if (!is.null(region_pack$bbox)) region_pack$bbox else sf::st_bbox(regionPoly)
  
  # Points df
  xy <- clean_df |>
    dplyr::select(gen, sp, long, lat) |>
    tidyr::drop_na() |>
    dplyr::mutate(id = dplyr::row_number())
  
  # Raster grid based on the requested bbox
  r <- terra::rast(
    xmin = bb["xmin"], xmax = bb["xmax"],
    ymin = bb["ymin"], ymax = bb["ymax"],
    resolution = cellsize,
    crs = "EPSG:4326"
  )
  
  pts <- terra::vect(xy[, c("long","lat")], geom = c("long","lat"), crs = "EPSG:4326")
  rcount <- terra::rasterize(pts, r, fun = "count", background = 0)
  rcount <- terra::mask(rcount, terra::vect(regionPoly))
  names(rcount) <- "numCollections"
  
  # Raster cells -> polygons, keep values, drop outer NA cells
  grid_v  <- terra::as.polygons(rcount, dissolve = FALSE, values = TRUE)
  grid_v  <- grid_v[!is.na(terra::values(grid_v)$numCollections), ]
  grid_sf <- sf::st_as_sf(grid_v)
  grid_sf$numCollections <- as.integer(grid_sf$numCollections)
  
  # Build sea mask = bbox polygon minus land
  bb_sfc   <- sf::st_as_sfc(bb)
  land     <- sf::st_union(regionPoly)
  sea_mask <- suppressWarnings(sf::st_difference(bb_sfc, land))
  
  # Bins & palette
  bin_levels <- c("0","1","2-25", "26-100","101-500","501-1000","1001-2000",">2000")
  grid_sf_binned <- grid_sf |>
    dplyr::mutate(
      collection_bin = cut(
        numCollections,
        breaks = c(-Inf, 0, 1, 25, 100, 500, 1000, 2000, Inf),
        labels = bin_levels,
        include.lowest = TRUE, right = TRUE
      ),
      collection_bin = factor(collection_bin, levels = bin_levels)
    )
  
  pal <- c(
    "0"          = "#FFFFFF",
    "1"       = "#EAF2FF",
    "2-25"     = "#D7E9FF",
    "26-100"    = "#C3DEFF",
    "101-500"   = "#AECFFF",
    "501-1000"  = "#FFF6BF",
    "1001-2000"  = "#FFE1B8",
    ">2000"      = "#FFB6C8"
  )
  
  # Graticule aligned to nice breaks, based on the requested bbox
  lon_step <- 10
  lat_step <- 5
  lon_min <- floor(bb["xmin"] / lon_step) * lon_step
  lon_max <- ceiling(bb["xmax"] / lon_step) * lon_step
  lat_min <- floor(bb["ymin"] / lat_step) * lat_step
  lat_max <- ceiling(bb["ymax"] / lat_step) * lat_step
  lon_breaks <- seq(lon_min, lon_max, by = lon_step)
  lat_breaks <- seq(lat_min, lat_max, by = lat_step)
  
  graticule <- sf::st_graticule(
    x   = bb,
    crs = sf::st_crs(regionPoly),
    lon = lon_breaks,
    lat = lat_breaks
  )
  grat_sea  <- suppressWarnings(sf::st_intersection(graticule, sea_mask))
  
  # Axis labels
  lab_lon <- function(x) paste0(abs(x), "°", ifelse(x < 0, "W", ifelse(x > 0, "E", "")))
  lab_lat <- function(y) paste0(abs(y), "°", ifelse(y < 0, "S", ifelse(y > 0, "N", "")))
  
  ocean_fill <- "grey90"
  
  # Plot
  p <- ggplot() +
    # heatmap cells
    geom_sf(data = grid_sf_binned, aes(fill = collection_bin), color = NA) +
    # sea mask to hide coastal overhang
    geom_sf(data = sea_mask, fill = ocean_fill, color = NA) +
    # graticule lines (only over sea so they don't cross land)
    geom_sf(data = grat_sea, linewidth = 0.5, color = "white") +
    # manual palette & legend
    scale_fill_manual(
      name   = "Digitised occurrences",
      values = pal,
      limits = bin_levels,
      drop   = FALSE,
      guide  = guide_legend(override.aes = list(colour = NA))
    ) +
    # axes aligned with graticule
    scale_x_continuous(breaks = lon_breaks, labels = lab_lon, expand = c(0, 0)) +
    scale_y_continuous(breaks = lat_breaks, labels = lab_lat, expand = c(0, 0)) +
    # ✅ lock panel to requested bbox
    coord_sf(
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"]),
      expand = FALSE, clip = "on"
    ) +
    # outlines
    geom_sf(data = regionPoly, fill = NA, colour = "grey60", linewidth = 0.2) +
    geom_sf(data = borders_clip, color = "grey60", linewidth = 0.2) +
    # theme
    theme(
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
  
  # Write SVG
  suffix <- if (!is.null(stamp)) paste0("_", stamp) else ""
  ggsave(
    filename = paste0(out_prefix, "_heatmap", suffix, ".svg"),
    plot = p, width = 8, height = 6, units = "in", dpi = 300
  )
  
  # Also return the plot object for interactive sessions
  invisible(p)
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
                              basis_filter = NULL,
                              bbox = NULL,                 # << NEW
                              lon_lim = NULL,              # << NEW (c(min_lon, max_lon))
                              lat_lim = NULL) {            # << NEW (c(min_lat, max_lat))
  run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  ttl <- stringr::str_to_title(taxon_label)
  bor_tag <- make_bor_tag(basis_filter)
  shp_tag <- make_shp_tag(study_shp_path)
  
  # --- NEW: build a bbox if user gave lon/lat limits
  user_bbox <- NULL
  if (!is.null(bbox)) {
    # allow unnamed numeric vector
    if (is.numeric(bbox) && length(bbox) == 4 && is.null(names(bbox))) {
      names(bbox) <- c("xmin","xmax","ymin","ymax")
    }
    user_bbox <- bbox
  } else if (!is.null(lon_lim) && !is.null(lat_lim)) {
    stopifnot(length(lon_lim) == 2, length(lat_lim) == 2)
    user_bbox <- c(xmin = min(lon_lim), xmax = max(lon_lim),
                   ymin = min(lat_lim), ymax = max(lat_lim))
  }
  message("=== ", toupper(taxon_label), " / taxonKey=", taxon_key,
          " / ", bor_tag, " / ", shp_tag,
          if (!is.null(user_bbox)) paste0(" / BBOX=[", paste(user_bbox, collapse=", "), "]") else "")
  
  # ---------------- Reuse-cleaned branch (unchanged, but pass region_pack) ---
  if (!is.null(cleaned_csv) && file.exists(cleaned_csv)) {
    clean <- readr::read_csv(cleaned_csv, show_col_types = FALSE)
    # OLD:
    # rp <- make_region_poly(bbox = user_bbox)
    # plot_obj <- make_heatmap(
    #   clean,
    #   out_prefix = paste0(taxon_label, "_", bor_tag, "_", shp_tag, "_digitised_occurrences"),
    #   region_pack = rp, stamp = run_stamp
    # )
    
    # NEW:
    write_island_svgs(
      clean_df   = clean,
      out_prefix = paste0(taxon_label, "_", bor_tag, "_", shp_tag, "_digitised_occurrences"),
      study_shp_path = study_shp_path
    )
    
    dens <- region_density(clean, region_col = "COUNTRY")
    dens_out <- paste0(taxon_label, "_occurrence_density_", bor_tag, "_", shp_tag, "_", run_stamp, ".csv")
    readr::write_csv(dens, dens_out)
    return(list(download_keys = NULL, standardized = NULL, cleaned = clean,
                heatmap = NULL, density_table = dens))
  } else if (isTRUE(reuse_cleaned)) {
    latest <- find_latest_cleaned_csv(taxon_label, bor_tag = bor_tag, shp_tag = shp_tag)
    if (!is.null(latest)) {
      clean <- readr::read_csv(latest, show_col_types = FALSE)
      # OLD:
      # rp <- make_region_poly(bbox = user_bbox)
      # plot_obj <- make_heatmap(
      #   clean,
      #   out_prefix = paste0(taxon_label, "_", bor_tag, "_", shp_tag, "_digitised_occurrences"),
      #   region_pack = rp, stamp = run_stamp
      # )
      
      # NEW:
      write_island_svgs(
        clean_df   = clean,
        out_prefix = paste0(taxon_label, "_", bor_tag, "_", shp_tag, "_digitised_occurrences"),
        study_shp_path = study_shp_path
      )
      
      dens <- region_density(clean, region_col = "COUNTRY")
      dens_out <- paste0(taxon_label, "_occurrence_density_", bor_tag, "_", shp_tag, "_", run_stamp, ".csv")
      readr::write_csv(dens, dens_out)
      return(list(download_keys = NULL, standardized = NULL, cleaned = clean,
                  heatmap = NULL, density_table = dens))
    }
  }
  
  # ---------------- Download/import/clean as before --------------------------
  if (is.null(download_keys_override)) {
    dl_keys <- occ_download_by_regions(taxon_key, scope)  # unchanged
  } else {
    dl_keys <- download_keys_override
    stopifnot(all(names(dl_keys) %in% names(REGIONS)))
  }
  raw <- import_and_bind(dl_keys)
  std <- standardize_for_modeling(raw)
  std <- filter_by_basis(std, basis_filter = basis_filter)
  
  std_out <- paste0("gbifSEAsia_", ttl, "OccDataStandardized_", bor_tag, "_", shp_tag, "_", run_stamp, ".csv")
  readr::write_csv(std, std_out)
  
  clean <- clean_occurrences(std, study_shp_path)
  clean_out <- paste0("gbifSEAsia_", ttl, "OccDataCleaned_", bor_tag, "_", shp_tag, "_", run_stamp, ".csv")
  readr::write_csv(clean, clean_out)
  
  # --- Build region pack with the requested bbox and plot --------------------
  # OLD:
  # rp <- make_region_poly(bbox = user_bbox)
  # plot_obj <- make_heatmap(
  #   clean,
  #   out_prefix = paste0(taxon_label, "_", bor_tag, "_", shp_tag, "_digitised_occurrences"),
  #   region_pack = rp, stamp = run_stamp
  # )
  
  # NEW:
  write_island_svgs(
    clean_df   = clean,
    out_prefix = paste0(taxon_label, "_", bor_tag, "_", shp_tag, "_digitised_occurrences"),
    study_shp_path = study_shp_path
  )
  
  
  dens <- region_density(clean, region_col = "COUNTRY")
  dens_out <- paste0(taxon_label, "_occurrence_density_", bor_tag, "_", shp_tag, "_", run_stamp, ".csv")
  readr::write_csv(dens, dens_out)
  
  list(download_keys = dl_keys, standardized = std, cleaned = clean,
       heatmap = NULL, density_table = dens)
}


# ============================== TAXON AND GEOGRAPHY ARGUMENTS ================#

mammal_keys <- list(
  papua_barat = "0024029-250920141307145",
  papua       = "0024027-250920141307145",
  png         = "0024057-250920141307145",
  kalimantan  = "0024055-250920141307145",
  sabah       = "0024038-250920141307145",
  sarawak     = "0024039-250920141307145",
  brunei      = "0024037-250920141307145",
  madagascar  = "0025618-251009101135966"
)

frog_keys <- list(
  papua_barat = "0020372-250920141307145",
  papua = "0020368-250920141307145", 
  png = "0022860-250920141307145",
  kalimantan = "0022979-250920141307145",
  sabah = "0020403-250920141307145",
  sarawak = "0020405-250920141307145",
  brunei = "0022904-250920141307145",
  madagascar = "0025544-251009101135966"
)

bird_keys <- list(
  papua_barat = "0024334-250920141307145",
  papua = "0024333-250920141307145", 
  png = "0024457-250920141307145",
  kalimantan = "0024396-250920141307145",
  sabah = "0024394-250920141307145",
  sarawak = "0024395-250920141307145",
  brunei = "0024387-250920141307145",
  madagascar = "0025640-251009101135966"
)

squamate_keys <- list(
  papua_barat = "0030972-250920141307145",
  papua = "0030971-250920141307145", 
  png = "0031079-250920141307145",
  kalimantan = "0031078-250920141307145",
  sabah = "0030995-250920141307145",
  sarawak = "0030997-250920141307145",
  brunei = "0030992-250920141307145",
  madagascar = "0025668-251009101135966"
)

plants_keys <- list(
  papua_barat = "0047000-250920141307145",
  papua = "0046996-250920141307145", 
  png = "0047318-250920141307145",
  kalimantan = "0047131-250920141307145",
  sabah = "0047007-250920141307145",
  sarawak = "0047029-250920141307145",
  brunei = "0047005-250920141307145",
  madagascar = "0025675-251009101135966"
)




# Frogs (Anura = 952)
frog_res <- run_gbif_pipeline(
  taxon_label = "anura", 
  taxon_key = 952, 
  scope = names(frog_keys),
  study_shp_path = "NGMadBor.shp",
  download_keys_override = frog_keys, 
  reuse_cleaned = FALSE,
  #basis_filter = c("PRESERVED_SPECIMEN"),
  #bbox = c(xmin = 95, xmax = 156, ymin = -10.3, ymax = 22)
)

#xmin = 95, xmax = 156, ymin = -10.3, ymax = 22 - IndoPacific
#xmin = 129, xmax = 156, ymin = -12, ymax = 4.5 - NG

mammal_res <- run_gbif_pipeline(
  taxon_label = "mammalia",
  taxon_key   = 359,
  scope       = names(mammal_keys),
  study_shp_path = "indoPacificIslands.shp",
  download_keys_override = mammal_keys,
  reuse_cleaned = FALSE,
  #basis_filter = c("PRESERVED_SPECIMEN")
  bbox = c(xmin = 95, xmax = 156, ymin = -10.3, ymax = 22)
)

bird_res <- run_gbif_pipeline(
  taxon_label = "aves",
  taxon_key   = 212,
  scope       = names(bird_keys),
  study_shp_path = "indoPacificIslands.shp",
  download_keys_override = bird_keys,
  reuse_cleaned = FALSE,
  bbox = c(xmin = 95, xmax = 156, ymin = -10.3, ymax = 22)
)

squamates_res <- run_gbif_pipeline(
  taxon_label = "squamates",
  taxon_key   = 11592253,
  scope       = names(squamate_keys),
  study_shp_path = "indoPacificIslands.shp",
  download_keys_override = squamate_keys,
  reuse_cleaned = FALSE,
  #basis_filter = c("PRESERVED_SPECIMEN")
  bbox = c(xmin = 95, xmax = 156, ymin = -10.3, ymax = 22)
)

plants_res <- run_gbif_pipeline(
  taxon_label = "tracheophyta",
  taxon_key   = 7707728,
  scope       = names(plants_keys),
  study_shp_path = "indoPacificIslands.shp",
  download_keys_override = plants_keys,
  reuse_cleaned = FALSE,
  #basis_filter = c("PRESERVED_SPECIMEN")
  bbox = c(xmin = 95, xmax = 156, ymin = -10.3, ymax = 22)
)
