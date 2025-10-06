# install.packages(c("sf","rnaturalearth","rnaturalearthdata","dplyr"))
library(sf)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)

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
