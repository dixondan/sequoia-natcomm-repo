"""
03b_topo.py
-----------
!! NOT REPRODUCIBLE — requires private raster data !!

Extracts topographic variables (solar radiation, TPI, slope) for each tree
by buffering the crown centroid at four radii (30, 60, 90, 120 m) and
computing the mean raster value within each buffer.

The output is pre-computed and provided in data/intermediate/. You do not
need to run this script to reproduce the modeling results.

Private inputs required
-----------------------
RSUN_RASTER  : solar radiation raster (rsun)
TPI_RASTER   : topographic position index raster
SLOPE_RASTER : slope raster

Shared input
------------
data/intermediate/trees_with_crowns.shp    from 01_prepare_crowns.py

Outputs (pre-computed, provided in repo)
----------------------------------------
data/intermediate/all-trees-topo-vars-30.csv
data/intermediate/all-trees-topo-vars-60.csv
data/intermediate/all-trees-topo-vars-90.csv
data/intermediate/all-trees-topo-vars-120.csv
"""

import os
import json
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio import mask

# ---------------------------------------------------------------------------
# Private paths -- update if re-running
# ---------------------------------------------------------------------------

RSUN_RASTER  = '/path/to/north-south-rsun.tif'   # <-- not shared
TPI_RASTER   = '/path/to/tpi_seki-3m.tif'        # <-- not shared
SLOPE_RASTER = '/path/to/seki-slope-10m.tif'     # <-- not shared

CROWNS_PATH = 'data/intermediate/trees_with_crowns.shp'
OUT_DIR     = 'data/intermediate'

os.makedirs(OUT_DIR, exist_ok=True)


# ---------------------------------------------------------------------------
# Extract topo variables at 30, 60, 90, 120 m buffers
# ---------------------------------------------------------------------------

df = gpd.read_file(CROWNS_PATH)
tree_ids = df['tree_id'].tolist()

for size in [30, 60, 90, 120]:
    l = []
    for tree_id in tree_ids:
        print(len(l))
        dft = df[df['tree_id'] == tree_id]
        dft_buf = dft['geometry'].centroid.buffer(size)

        crown_json = dft_buf.to_json()
        geom2clip = [feature['geometry'] for feature in json.loads(crown_json)['features']]

        with rasterio.open(TPI_RASTER) as src:
            arr_tpi, trans = mask.mask(src, geom2clip, crop=True)
            arr_tpi = arr_tpi[arr_tpi > -9999]
            tpi = arr_tpi.mean()

        with rasterio.open(RSUN_RASTER) as src:
            arr_solar, trans = mask.mask(src, geom2clip, crop=True)
            arr_solar = arr_solar[arr_solar > -9999]
            solar = arr_solar.mean()

        with rasterio.open(SLOPE_RASTER) as src:
            arr_slope, trans = mask.mask(src, geom2clip, crop=True)
            arr_slope = arr_slope[arr_slope > -9999]
            slope = arr_slope.mean()

        l.append((tree_id, solar, tpi, slope))

    dfout = pd.DataFrame(l, columns=['tree_id', 'solar', 'tpi', 'slope'])
    out_path = os.path.join(OUT_DIR, f'all-trees-topo-vars-{size}.csv')
    dfout.to_csv(out_path, index=False)
    print(f'Done ({size} m buffer): {len(dfout)} trees saved to {out_path}')