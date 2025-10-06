# ---- Quick viewer for .shp territory borders ----

library(sf)
library(dplyr)
library(ggplot2)
library(tmap)
library(mapview)

# EDIT THIS: path to your shapefile (.shp)
shp_path <- "territory_selection.shp"


# 1) Read the shapefile
#   - stringsAsFactors=FALSE: keep attributes as characters
#   - quiet=TRUE: reduce console noise
borders <- st_read(shp_path, stringsAsFactors = FALSE, quiet = TRUE)

# 2) Basic checks
cat("Geometry type:", unique(st_geometry_type(borders)), "\n")
cat("Features:", nrow(borders), "\n")
print(st_crs(borders))

# 3) Fix common geometry issues and standardize CRS to WGS84 (EPSG:4326)
borders <- st_make_valid(borders)
if (is.na(st_crs(borders))) {
  warning("Layer has no CRS defined. If you know it, set it, e.g.: st_crs(borders) <- 3857")
}
borders_4326 <- tryCatch(st_transform(borders, 4326), error = function(e) borders)

# 4) Quick static plot (ggplot2)
#    If you have a name column (e.g., "name" or "admin"), set label_col accordingly.
label_col <- intersect(names(borders_4326), c("name", "admin", "NAME", "NAME_EN"))[1]

p <- ggplot(borders_4326) +
  geom_sf(fill = NA, linewidth = 0.4) +
  labs(title = "Territory Borders", subtitle = shp_path) +
  theme_minimal()

print(p)

# 5) Interactive map in RStudio Viewer (mapview)
#    - Automatically picks a sensible basemap
#    - Hover to see attributes; set 'zcol' to color by a field if you like
z_field <- if (!is.null(label_col)) label_col else NA
mapviewOptions(fgb = TRUE)  # fast rendering
mv <- mapview(borders_4326, zcol = z_field, layer.name = "Territory Borders")
mv  # appears in the Viewer pane

# 6) (Optional) Interactive tmap with basemap tiles
tmap_mode("view")
tm <- tm_shape(borders_4326) +
  tm_polygons(alpha = 0, border.col = "black") +
  tm_view(view.legend.position = c("left", "bottom"))
tm
