# ---- Quick viewer for .shp territory borders ----

library(sf)
library(dplyr)
library(ggplot2)

# EDIT THIS: path to your shapefile (.shp)
shp_path <- "PNGIDP.shp"


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
