# extract-crown-metrics.R
# -----------------------
# Extracts LiDAR-derived structural metrics for tree crowns in a single
# survey tile. Called from 03a_lidar.py via subprocess for each ns-f tile.
#
# Usage (from command line):
#   Rscript extract-crown-metrics.R <ns> <f>
#
# Arguments:
#   ns  : survey direction, either 'north' or 'south'
#   f   : tile identifier (e.g. '338000_4051000')
#
# Private inputs required (not shared):
#   - LiDAR point cloud tiles (.laz)
#   - Tile neighbor lookup shapefiles
#   - Crown polygons: data/labels/level_2/all_trees_match-unmatched-wfire-n27751.shp
#
# Outputs (written to LIDAR_VARS_DIR):
#   {ns}-{f}-10-ladder.shp
#   {ns}-{f}-{d1}-{d2}-rumple.shp   for d2 in 30, 60, 90, 120
#   {ns}-{f}-{d1}-{d2}-cdensity.shp for d2 in 30, 60, 90, 120

rm(list = ls(globalenv()))

library(lidRmetrics)
library(lidR)
library(sf)
library(dplyr)

options(warn = -1)

# ---------------------------------------------------------------------------
# Private paths -- update if re-running
# ---------------------------------------------------------------------------

CROWN_FILE    <- 'data/intermediate/trees_with_crowns.shp'          # <-- update
LAZ_DIR       <- '/path/to/laz_chm'                                  # <-- not shared
TILES_DIR     <- '/path/to/tiles_and_neighbors'                      # <-- not shared
LIDAR_VARS_DIR <- '/path/to/database/lidarvars'                      # <-- not shared

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
ns <- args[1]
f  <- args[2]

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

function_sum <- function(z) length(z)

cdensity <- function(z) sum(z > 4) / length(z)  # 4 m canopy height cutoff

create_rings <- function(crowns, distance1, distance2) {
  # Creates donut-shaped buffer rings (outer - inner) around each crown
  polygon_list <- vector('list', nrow(crowns))
  for (i in seq_len(nrow(crowns))) {
    crown        <- crowns[i, ]
    inner_buffer <- st_buffer(crown, distance1)
    outer_buffer <- st_buffer(crown, distance2)
    polygon_list[[i]] <- st_difference(outer_buffer, inner_buffer)
  }
  rings <- do.call(rbind, polygon_list)
  rings[, c('tree_id', 'geometry')]
}


extract_ladder_metrics <- function(buffers, las) {
  # Extracts ladder fuel metrics for each crown buffer:
  #   zlt4    : returns below 4 m
  #   metrics2: returns between 4 m and zq95/2 (ladder fuels)
  #   metrics3: returns above zq95/2
  #   metrics4: returns above zq95
  results_list <- vector('list', nrow(buffers))

  for (i in seq_len(nrow(buffers))) {
    print(i)
    buffer <- buffers[i, ]

    las_clip  <- clip_roi(las, buffer)
    crown     <- st_buffer(buffer, -10)
    las_inner <- clip_roi(las_clip, crown)
    zq95      <- quantile(las_inner$Z, 0.95, na.rm = TRUE)
    zq50      <- zq95 / 2
    tree_id   <- buffer$tree_id

    # Returns below 4 m
    las_under  <- las_clip[las_clip$Z < 4]
    metrics1   <- polygon_metrics(las_under, func = ~function_sum(Z), geometry = buffer)
    colnames(metrics1) <- c('zlt4', 'geometry')

    # Returns between 4 m and zq50 (ladder fuels)
    metrics2 <- 0
    las_ladder <- las_clip[las_clip$Z > 4 & las_clip$Z < zq50]
    if (nrow(las_ladder@data) > 0) {
      m2       <- polygon_metrics(las_ladder, func = ~function_sum(Z), geometry = buffer)
      metrics2 <- m2$V1
    }

    # Returns above zq50
    las_gt50  <- las_clip[las_clip$Z > zq50]
    metrics3  <- polygon_metrics(las_gt50, func = ~function_sum(Z), geometry = buffer)$V1

    # Returns above zq95
    las_gt95  <- las_clip[las_clip$Z > zq95]
    metrics4  <- polygon_metrics(las_gt95, func = ~function_sum(Z), geometry = buffer)$V1

    crown_metrics          <- cbind(metrics1, metrics2, metrics3, metrics4)
    crown_metrics$tree_id  <- tree_id
    crown_metrics$zq95     <- zq95
    results_list[[i]]      <- crown_metrics
  }

  do.call(rbind, results_list)
}


extract_other_metrics <- function(rings, las, func) {
  # Applies an arbitrary lidR metric function to each ring polygon
  results_list <- vector('list', nrow(rings))
  for (i in seq_len(nrow(rings))) {
    print(i)
    ring              <- rings[i, ]
    metrics           <- polygon_metrics(las, func, geometry = ring)
    results_list[[i]] <- cbind(ring, metrics)
  }
  do.call(rbind, results_list)
}

# ---------------------------------------------------------------------------
# Main extraction function
# ---------------------------------------------------------------------------

get_neighborhood_metrics <- function(ns, f) {

  # Load crowns for this tile
  crowns <- st_read(CROWN_FILE, quiet = TRUE)
  crowns <- crowns[crowns$ns == ns & crowns$f == f, c('tree_id', 'geometry')]

  # Load tile neighbor list and read + mosaic LiDAR tiles
  tiles_file <- file.path(TILES_DIR, paste0(ns, '-', f, '.shp'))
  lidar_tiles <- st_read(tiles_file, quiet = TRUE)
  tiles_to_read <- as.list(lidar_tiles$file_chm)

  center_las_path <- file.path(LAZ_DIR, paste0(f, '-', ns, '.laz'))
  las <- readLAS(center_las_path, select = 'xyzric')

  tiles_to_read <- setdiff(tiles_to_read, center_las_path)
  tile_geom     <- st_buffer(lidar_tiles[lidar_tiles$f == f, ], dist = 120)

  if (length(tiles_to_read) > 0) {
    cat('Reading neighboring tiles\n')
    for (las_file in tiles_to_read) {
      las_outer <- readLAS(las_file, select = 'xyzric')
      las_outer <- clip_roi(las_outer, tile_geom)
      las       <- rbind(las, las_outer)
    }
  }

  # Ladder fuels (10 m buffer around crown)
  cat('Extracting ladder metrics\n')
  outname_ladder <- file.path(LIDAR_VARS_DIR, sprintf('%s-%s-10-ladder.shp', ns, f))
  if (!file.exists(outname_ladder)) {
    buffers        <- st_buffer(crowns, 10)
    ladder_metrics <- extract_ladder_metrics(buffers, las)
    st_write(ladder_metrics, outname_ladder, append = FALSE, quiet = TRUE)
  }

  # Rumple and canopy density at four buffer distances
  ring_list <- list(c(0, 30), c(0, 60), c(0, 90), c(0, 120))

  cat('Extracting rumple metrics\n')
  for (rd in ring_list) {
    d1 <- rd[1]; d2 <- rd[2]
    outname <- file.path(LIDAR_VARS_DIR, sprintf('%s-%s-%d-%d-rumple.shp', ns, f, d1, d2))
    if (!file.exists(outname)) {
      rings <- create_rings(crowns, d1, d2)
      st_write(extract_other_metrics(rings, las, func = .metrics_rumple),
               outname, append = FALSE, quiet = TRUE)
    }
  }

  cat('Extracting canopy density metrics\n')
  las_first <- las[las$ReturnNumber == 1]
  for (rd in ring_list) {
    d1 <- rd[1]; d2 <- rd[2]
    outname <- file.path(LIDAR_VARS_DIR, sprintf('%s-%s-%d-%d-cdensity.shp', ns, f, d1, d2))
    if (!file.exists(outname)) {
      rings <- create_rings(crowns, d1, d2)
      st_write(extract_other_metrics(rings, las_first, func = ~cdensity(Z)),
               outname, append = FALSE, quiet = TRUE)
    }
  }
}

get_neighborhood_metrics(ns, f)