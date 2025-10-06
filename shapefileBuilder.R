# install.packages(c("sf","rnaturalearth","rnaturalearthdata","dplyr"))
library(sf)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)

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

