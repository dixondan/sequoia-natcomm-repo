"""
03a_lidar.py
------------
!! NOT REPRODUCIBLE — requires private LiDAR tile data !!

Extracts LiDAR-derived structural variables (rumple, canopy density, ladder
fuels) for each tree crown, using per-tile shapefiles produced by the R script
Rscripts/extract-crown-metrics.R. One set of variables is extracted at four
height ranges (0-30, 0-60, 0-90, 0-120 m).

The output is pre-computed and provided in data/intermediate/. You do not
need to run this script to reproduce the modeling results.

Private inputs required
-----------------------
LIDAR_VARS_DIR  : directory of per-tile LiDAR variable shapefiles
                  (rumple, cdensity, ladder), one set per ns-f tile
RSCRIPT         : path to Rscript executable
R_METRICS_SCRIPT: Rscripts/extract-crown-metrics.R

Shared input
------------
data/intermediate/trees_with_crowns.shp    from 01_prepare_crowns.py

Output (pre-computed, provided in repo)
-------
data/intermediate/all-trees-lidar-vars.csv
"""

import os
import subprocess
import pandas as pd
import geopandas as gpd

# ---------------------------------------------------------------------------
# Private paths -- update if re-running
# ---------------------------------------------------------------------------

LIDAR_VARS_DIR   = '/path/to/database/lidarvars'           # <-- not shared
RSCRIPT          = '/path/to/Rscript'                      # <-- not shared
R_METRICS_SCRIPT = 'Rscripts/extract-crown-metrics.R'

CROWNS_PATH = 'data/intermediate/trees_with_crowns.shp'
OUT_PATH    = 'data/intermediate/all-trees-lidar-vars.csv'

os.makedirs('data/intermediate', exist_ok=True)


# ---------------------------------------------------------------------------
# Section 1: Extract crown metrics per tile via R script
# ---------------------------------------------------------------------------

def process_metrics(ns, f):
    outshp = os.path.join(LIDAR_VARS_DIR, f'{ns}-{f}-0-120-cdensity.shp')
    if not os.path.exists(outshp):
        print(f'Processing crown metrics: {ns}-{f}')
        subprocess.run([RSCRIPT, R_METRICS_SCRIPT, ns, f])

df = gpd.read_file(CROWNS_PATH)
df['ns_f'] = df['ns'] + '-' + df['f']
items = list(df['ns_f'].unique())

for i, item in enumerate(items):
    print(i, item)
    ns, f = item.split('-')
    process_metrics(ns, f)


# ---------------------------------------------------------------------------
# Section 2: Compile lidar variables across all tiles
# ---------------------------------------------------------------------------

df = gpd.read_file(CROWNS_PATH)
df['ns_f'] = df['ns'] + '-' + df['f']
items = list(df['ns_f'].unique())

l2 = []
for item in items:
    distances = [(0, 30), (0, 60), (0, 90), (0, 120)]
    l = []
    for d1, d2 in distances:
        ns, f = item.split('-')

        # rumple
        crowns_rumple = gpd.read_file(
            os.path.join(LIDAR_VARS_DIR, f'{ns}-{f}-{d1}-{d2}-rumple.shp'))
        rumple_col = f'rumple_{d1}_{d2}'
        crowns_rumple[rumple_col] = crowns_rumple['rumple']
        crowns_rumple = crowns_rumple[['tree_id', rumple_col]]

        # canopy density
        crowns_cdensity = gpd.read_file(
            os.path.join(LIDAR_VARS_DIR, f'{ns}-{f}-{d1}-{d2}-cdensity.shp'))
        cdensity_col = f'cdensity_{d1}_{d2}'
        crowns_cdensity[cdensity_col] = crowns_cdensity['V1']
        crowns_cdensity = crowns_cdensity[['tree_id', cdensity_col]]

        dft = pd.merge(crowns_rumple, crowns_cdensity, on='tree_id')
        l.append(dft)

    merged_df_fuels = l[0]
    for dft in l[1:]:
        merged_df_fuels = merged_df_fuels.merge(dft, on='tree_id')

    # ladder fuels
    ladder_df = gpd.read_file(
        os.path.join(LIDAR_VARS_DIR, f'{ns}-{f}-10-ladder.shp'))
    merged_df_fuels = merged_df_fuels.merge(ladder_df, on='tree_id')

    l2.append(merged_df_fuels)

dfout = pd.concat(l2)
dfout.to_csv(OUT_PATH, index=False)
print(f'Done: {len(dfout)} trees saved to {OUT_PATH}')