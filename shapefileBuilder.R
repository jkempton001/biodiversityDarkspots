# install.packages(c("sf","rnaturalearth","rnaturalearthdata","dplyr"))
library(sf)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(rnaturalearthhires)

#=========================== INDO-PACIFC ISLANDS ==============================#

# ------------------------------------------------------------------------------
# 1) Pull admin boundaries
# ------------------------------------------------------------------------------

# Countries at 1:10m (highest detail in Natural Earth via rnaturalearth)
ctry <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf") %>%
  st_make_valid()

# Keep: Indonesia, Brunei, Philippines, Papua New Guinea
sel_ctry <- ctry %>%
  filter(iso_a3 %in% c("IDN","BRN","PHL","PNG","TLS"))

# Malaysian Borneo = Sabah, Sarawak, (and Labuan which is on Borneo)
mys_states <- rnaturalearth::ne_states(country = "Malaysia", returnclass = "sf") %>%
  st_make_valid()

mys_borneo <- mys_states %>%
  filter(name %in% c("Sabah","Sarawak","W.P. Labuan")) %>%
  summarise(geometry = st_union(geometry)) %>%      # dissolve into 1 polygon
  mutate(src = "MY_Borneo")

# ------------------------------------------------------------------------------
# 2) Get the New Guinea island polygon (independent of admin borders)
# ------------------------------------------------------------------------------

# Natural Earth has a named island polygon for "New Guinea"
# (type = "geography_regions_polys" holds major islands)
islands <- rnaturalearth::ne_download(
  scale = "large",
  type  = "geography_regions_polys",
  category = "physical",
  returnclass = "sf"
) %>% st_make_valid()

new_guinea <- islands %>%
  filter(NAME == "NEW GUINEA") %>%
  summarise(geometry = st_union(geometry)) %>%
  mutate(src = "New_Guinea")

# ------------------------------------------------------------------------------
# 3) Combine everything and dissolve
# ------------------------------------------------------------------------------

combo <- bind_rows(
  sel_ctry %>% select(src = name_long, geometry),
  mys_borneo %>% select(src, geometry),
  new_guinea %>% select(src, geometry)
) %>%
  # Keep land-only geometries (ctry/islands are already land; this is just a safety)
  st_make_valid()

# Dissolve to a single (multi)polygon
territory <- combo %>%
  summarise(geometry = st_union(geometry)) %>%
  st_make_valid()

# ------------------------------------------------------------------------------
# 4) (Optional) Reproject to an equal-area CRS for later spatial ops
#    This won’t change the shape when you save to WGS84 again;
#    but it’s useful if you’ll do area-based calculations.
# ------------------------------------------------------------------------------
# Asia/Australia equal-area CRS example (World Cylindrical Equal Area, EPSG: 6933)
territory_eq <- st_transform(territory, 6933)

ggplot(territory) +
  geom_sf(fill = "white", color = "darkgreen") +
  theme_minimal() +
  labs(title = "Territory Mask (hi-res)")

# ------------------------------------------------------------------------------
# 5) Save as a shapefile (or, better, a GeoPackage)
# ------------------------------------------------------------------------------

# Shapefile (creates multiple files in this folder)
st_write(territory, "territory_selection.shp", delete_layer = TRUE)

# Strongly recommended alternative (single-file, modern format):
# st_write(territory, "territory_selection.gpkg", layer = "territory", delete_layer = TRUE)

#================================== NEW GUINEA ================================#

# ------------------------------------------------------------------------------
# 1) Countries: keep PNG
# ------------------------------------------------------------------------------
ctry <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf") %>%
  st_make_valid()

png_ctry <- ctry %>%
  filter(iso_a3 == "PNG") %>%
  select(iso_a3, name_long, geometry) %>%
  mutate(src = "Papua New Guinea")

# ------------------------------------------------------------------------------
# 2) Indonesian provinces: return to legacy 2-province scheme
# ------------------------------------------------------------------------------
idn_states <- rnaturalearth::ne_states(country = "Indonesia", returnclass = "sf") %>%
  st_make_valid()

pick_names <- function(df, nm) {
  cols <- intersect(c("name", "name_en"), names(df))
  df %>%
    dplyr::filter(dplyr::if_any(dplyr::all_of(cols), ~ .x %in% nm))
}


# Try direct old names first
papua_old       <- pick_names(idn_states, c("Papua"))
papua_barat_old <- pick_names(idn_states, c("Papua Barat", "West Papua"))

# If either missing, aggregate from the newer six-province split
if (nrow(papua_old) == 0 | nrow(papua_barat_old) == 0) {
  papua_group_names <- c("Papua", "Papua Tengah", "Central Papua",
                         "Papua Pegunungan", "Highland Papua",
                         "Papua Selatan", "South Papua")
  papua_barat_group_names <- c("Papua Barat", "West Papua",
                               "Papua Barat Daya", "Southwest Papua")
  
  papua_grp <- pick_names(idn_states, papua_group_names) %>%
    summarise(geometry = st_union(geometry), .groups = "drop") %>%
    mutate(name = "Papua")
  
  papua_barat_grp <- pick_names(idn_states, papua_barat_group_names) %>%
    summarise(geometry = st_union(geometry), .groups = "drop") %>%
    mutate(name = "Papua Barat")
  
  id_papua2 <- bind_rows(papua_grp, papua_barat_grp) %>%
    mutate(src = paste0("IDN - ", name)) %>%
    select(src, geometry) %>%
    st_make_valid()
} else {
  # Old scheme present — just keep as-is (dissolve to one feature per province)
  id_papua2 <- bind_rows(
    papua_old  %>% mutate(name = "Papua"),
    papua_barat_old %>% mutate(name = "Papua Barat")
  ) %>%
    group_by(name) %>%
    summarise(geometry = st_union(geometry), .groups = "drop") %>%
    mutate(src = paste0("IDN - ", name)) %>%
    select(src, geometry) %>%
    st_make_valid()
}

# ------------------------------------------------------------------------------
# 3) Combine PNG + legacy two Papuan provinces
# ------------------------------------------------------------------------------
combo <- bind_rows(
  png_ctry %>% select(src, geometry),
  id_papua2 %>% select(src, geometry)
) %>% st_make_valid()

# Optional: dissolve to single mask (comment out to keep three features)
territory_mask <- combo %>%
  summarise(geometry = st_union(geometry)) %>%
  st_make_valid()

# ------------------------------------------------------------------------------
# 4) Plot
# ------------------------------------------------------------------------------
ggplot() +
  geom_sf(data = combo, aes(fill = src), color = "gray25", linewidth = 0.25) +
  theme_minimal() +
  labs(title = "Papua New Guinea + Legacy Indonesian Papua Provinces")

# ------------------------------------------------------------------------------
# 5) Save
# ------------------------------------------------------------------------------
# Recommended: per-feature GeoPackage (PNG + Papua + Papua Barat)
#st_write(combo, "png_papua_legacy.gpkg", layer = "png_papua_legacy", delete_layer = TRUE)

# If you want a single dissolved mask instead:
# st_write(territory_mask, "png_papua_legacy_mask.gpkg", layer = "mask", delete_layer = TRUE)

# Shapefile option if required:
st_write(combo, "PNGIDP.shp", delete_layer = TRUE)

## ==================== NG, BORNEO, and MADGASCAR ============================##

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
st_valid <- function(g) suppressWarnings(st_make_valid(g))

# Try to select admin-1 rows by ISO-3166-2 code first, else by name
pick_admin1 <- function(df, iso_codes = character(), names_fallback = character()) {
  # Potential code fields found in Natural Earth / GADM-like exports
  code_cols <- intersect(c("iso_3166_2", "postal", "postal_code", "code_hasc", "gn_a1_code"), names(df))
  name_cols <- intersect(c("name", "name_en", "name_long", "name_alt"), names(df))
  
  sel <- tibble()
  if (length(code_cols) && length(iso_codes)) {
    sel <- df %>% filter(reduce(code_cols, ~ .x | .y %in% iso_codes, .init = (.[[code_cols[1]]] %in% iso_codes)))
  }
  if (nrow(sel) == 0 && length(names_fallback) && length(name_cols)) {
    sel <- df %>% filter(reduce(name_cols, ~ .x | .y %in% names_fallback, .init = (.[[name_cols[1]]] %in% names_fallback)))
  }
  st_valid(sel)
}

# ------------------------------------------------------------
# Load countries and states (Natural Earth, 1:10m = "large")
# ------------------------------------------------------------
ctry <- ne_countries(scale = "large", returnclass = "sf") |> st_valid()
# We'll need states/provinces for Indonesia & Malaysia
idn_states <- ne_states(country = "Indonesia", returnclass = "sf") |> st_valid()
mys_states <- ne_states(country = "Malaysia",  returnclass = "sf") |> st_valid()

# ------------------------------------------------------------
# Country-level pieces (whole-island countries)
# - BRN (Brunei) contributes to Borneo
# - PNG (Papua New Guinea) contributes to New Guinea
# - MDG (Madagascar) is the entire island
# ------------------------------------------------------------
brn <- ctry %>% filter(iso_a3 == "BRN") %>% select(src = name_long, geometry)
png <- ctry %>% filter(iso_a3 == "PNG") %>% select(src = name_long, geometry)
mdg <- ctry %>% filter(iso_a3 == "MDG") %>% select(src = name_long, geometry)

# ------------------------------------------------------------
# Province-level pieces
# Borneo (Kalimantan Indonesia + Malaysia Borneo states)
# ------------------------------------------------------------
# Indonesia (Kalimantan) — ISO 3166-2 preferred
id_kalim_iso <- c("ID-KB","ID-KT","ID-KS","ID-KI","ID-KU")
id_kalim_names <- c("Kalimantan Barat","Kalimantan Tengah","Kalimantan Selatan",
                    "Kalimantan Timur","Kalimantan Utara")
id_kalim <- pick_admin1(idn_states, iso_codes = id_kalim_iso, names_fallback = id_kalim_names) %>%
  mutate(src = paste0("IDN - ", coalesce(name, name_en)))

# Malaysia (Borneo) — Sabah, Sarawak, Labuan
my_borneo_iso <- c("MY-12","MY-13","MY-15")  # Sarawak, Sabah, Labuan
my_borneo_names <- c("Sarawak","Sabah","W.P. Labuan","Labuan")
my_borneo <- pick_admin1(mys_states, iso_codes = my_borneo_iso, names_fallback = my_borneo_names) %>%
  mutate(src = paste0("MYS - ", coalesce(name, name_en)))

# ------------------------------------------------------------
# Province-level pieces
# New Guinea (Indonesian Papua provinces)
# Support both the legacy 2-province scheme and the newer 5-province split
# ------------------------------------------------------------
# Legacy two:
id_papua_legacy_iso   <- c("ID-PA","ID-PB")                  # Papua, Papua Barat
id_papua_legacy_names <- c("Papua","Papua Barat","West Papua")

# Newer split (as ISO codes and common English/Indo names)
id_papua_new_iso <- c("ID-PT","ID-PS","ID-PP","ID-PD","ID-PB") # Central, South, Highland, Southwest, West Papua
id_papua_new_names <- c("Papua Tengah","Central Papua",
                        "Papua Selatan","South Papua",
                        "Papua Pegunungan","Highland Papua",
                        "Papua Barat Daya","Southwest Papua",
                        "Papua Barat","West Papua",
                        "Papua")  # include "Papua" as some datasets still carry it

# Try to get whichever set exists in the data
id_papua <- pick_admin1(idn_states,
                        iso_codes = c(id_papua_new_iso, id_papua_legacy_iso),
                        names_fallback = c(id_papua_new_names, id_papua_legacy_names)) %>%
  mutate(src = paste0("IDN - ", coalesce(name, name_en)))

# ------------------------------------------------------------
# Assemble islands
# ------------------------------------------------------------
# Borneo = Brunei + Malaysian Borneo states + Indonesian Kalimantan provinces
borneo_parts <- bind_rows(
  brn %>% select(src, geometry),
  my_borneo %>% select(src, geometry),
  id_kalim %>% select(src, geometry)
) %>% st_valid()

borneo <- borneo_parts %>%
  summarise(island = "BORNEO", geometry = st_union(geometry), .groups = "drop") %>%
  st_valid()

# New Guinea = PNG + Indonesian Papua provinces
new_guinea_parts <- bind_rows(
  png %>% select(src, geometry),
  id_papua %>% select(src, geometry)
) %>% st_valid()

new_guinea <- new_guinea_parts %>%
  summarise(island = "NEW GUINEA", geometry = st_union(geometry), .groups = "drop") %>%
  st_valid()

# Madagascar = whole country
madagascar <- mdg %>%
  summarise(island = "MADAGASCAR", geometry = st_union(geometry), .groups = "drop") %>%
  st_valid()

# Combine to 3 features
three_islands <- bind_rows(borneo, new_guinea, madagascar) %>%
  st_set_crs(4326) %>%
  st_valid()

# ------------------------------------------------------------
# Quick visual check
# ------------------------------------------------------------
ggplot(three_islands) +
  geom_sf(aes(fill = island), color = "gray25", linewidth = 0.25, show.legend = FALSE) +
  theme_minimal() +
  labs(title = "Borneo, New Guinea, Madagascar (built from country + province codes)")

# ------------------------------------------------------------
# Write out
# ------------------------------------------------------------
st_write(three_islands, "NGMadBor.shp", delete_layer = TRUE)
# Alternative modern container:
# st_write(three_islands, "three_islands_by_codes.gpkg", layer = "islands", delete_layer = TRUE)

## ============================= Rainforest regions ==========================##

## ==================== AMAZON, CONGO BASIN, INDO-MALAY ======================##

# ========================= Rainforest regions (hi-res) ======================= #

st_valid <- function(g) suppressWarnings(st_make_valid(g))

# ------------------------------------------------------------------------------
# 1) Countries at 1:10m (highest detail available from Natural Earth)
# ------------------------------------------------------------------------------
# Use numeric scale = 10 to be explicit; requires rnaturalearthhires to be loaded
ctry <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf") %>%
  st_valid()

# ------------------------------------------------------------------------------
# 2) Region membership via ISO3 (your lists)
# ------------------------------------------------------------------------------
AMAZON_ISO3     <- c("BRA","PER","COL","ECU","BOL","VEN","GUY","SUR","GUF")
CONGOBASIN_ISO3 <- c("COD","COG","GAB","CMR","CAF","GNQ","AGO","RWA","BDI")
INDOMALAY_ISO3  <- c("IDN","MYS","BRN","SGP","PHL","THA","VNM","KHM","LAO","MMR","TLS","PNG")

# ------------------------------------------------------------------------------
# 3) (Optional) Restrict Angola to Cabinda or northern piece only
#     This keeps outlines sharp and avoids including Angola’s far south.
#     We try an admin-level pick for "Cabinda"; if not found, we crop with a bbox.
# ------------------------------------------------------------------------------
maybe_cabinda <- function(x) {
  # try admin name first (works with NE states)
  ang_states <- tryCatch(
    rnaturalearth::ne_states(country = "Angola", returnclass = "sf") %>% st_valid(),
    error = function(e) NULL
  )
  if (!is.null(ang_states) && any(grepl("^Cabinda$", ang_states$name))) {
    cab <- ang_states %>% filter(name == "Cabinda") %>% st_union() %>% st_make_valid()
    return(cab)
  } else {
    # coarse north-only fallback (tweak as you like)
    bb <- st_as_sfc(st_bbox(c(xmin = 10.5, xmax = 24.0, ymin = -8.0, ymax = -4.0), crs = 4326))
    return(bb)
  }
}

# ------------------------------------------------------------------------------
# 4) Build each region (union of hi-res countries)
#     Only topology fix is st_make_valid() (no snapping/rounding).
# ------------------------------------------------------------------------------
amazon <- ctry %>%
  filter(iso_a3 %in% AMAZON_ISO3) %>%
  summarise(region = "AMAZON", geometry = st_union(geometry), .groups = "drop") %>%
  st_valid()

congo_raw <- ctry %>%
  filter(iso_a3 %in% CONGOBASIN_ISO3) %>%
  st_valid()

# Optionally trim Angola so only Cabinda/northern part contributes
if ("AGO" %in% CONGOBASIN_ISO3 && any(congo_raw$iso_a3 == "AGO")) {
  ang_mask <- maybe_cabinda(congo_raw)
  congo_raw <- congo_raw %>%
    mutate(geometry = ifelse(
      iso_a3 == "AGO",
      st_intersection(geometry, ang_mask),
      geometry
    )) %>%
    st_valid()
}

congo_basin <- congo_raw %>%
  summarise(region = "CONGO BASIN", geometry = st_union(geometry), .groups = "drop") %>%
  st_valid()

indomalay <- ctry %>%
  filter(iso_a3 %in% INDOMALAY_ISO3) %>%
  summarise(region = "INDO-MALAY", geometry = st_union(geometry), .groups = "drop") %>%
  st_valid()

# ------------------------------------------------------------------------------
# 5) Combine to 3 features
# ------------------------------------------------------------------------------
three_rainforests <- bind_rows(amazon, congo_basin, indomalay) %>%
  st_set_crs(4326) %>%
  st_valid()

# ------------------------------------------------------------------------------
# 6) Quick visual check
# ------------------------------------------------------------------------------
ggplot(three_rainforests) +
  geom_sf(aes(fill = region), color = "gray25", linewidth = 0.25) +
  theme_minimal() +
  labs(title = "Major Tropical Rainforest Regions (country-based unions, 1:10m)",
       fill = "Region")

# ------------------------------------------------------------------------------
# 7) Write out
# ------------------------------------------------------------------------------
st_write(three_rainforests, "rainforestRegions.shp", delete_layer = TRUE)
# Or:
# st_write(three_rainforests, "rainforest_regions.gpkg", layer = "regions", delete_layer = TRUE)
