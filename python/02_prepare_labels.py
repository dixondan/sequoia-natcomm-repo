"""
02_prepare_labels.py
--------------------
Extracts pixel-level mortality and survival predictions from the KNP and
Castle fire prediction rasters for each tree crown polygon. Outputs a CSV
with per-tree mortality probability used as the response variable in modeling.

Inputs
------
data/intermediate/trees_with_crowns.shp     from 01_prepare_crowns.py
data/raw/pred-knp_2020-2022.tif
data/raw/pred-castle_2019-2021.tif

Output
------
data/intermediate/tree-mortality-outcomes.csv
"""

import os
import json
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio import mask

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

CROWNS_PATH   = 'data/intermediate/trees_with_crowns.shp'
RASTER_KNP    = 'data/raw/pred-knp_2020-2022.tif'
RASTER_CASTLE = 'data/raw/pred-castle_2019-2021.tif'
OUT_PATH      = 'data/intermediate/tree-mortality-outcomes.csv'

os.makedirs('data/intermediate', exist_ok=True)

# ---------------------------------------------------------------------------
# Extract mortality predictions from rasters
# ---------------------------------------------------------------------------

df = gpd.read_file(CROWNS_PATH)
df = df.to_crs('EPSG:5070')
tree_ids = df['tree_id'].tolist()

l = []
for tree_id in tree_ids:
    dft = df[df['tree_id'] == tree_id]
    crown_json = dft.to_json()
    geom2clip = [feature['geometry'] for feature in json.loads(crown_json)['features']]
    fire = dft['FIRE_NAME'].unique()[0]

    if fire == 'KNP':
        raster = RASTER_KNP
    elif fire == 'CASTLE':
        raster = RASTER_CASTLE

    with rasterio.open(raster) as src:
        arr, trans = mask.mask(src, geom2clip, crop=True)
        sum_0 = np.sum(arr, axis=0)
        arr_mask = np.zeros((sum_0.shape[0], sum_0.shape[1]), dtype=np.float32)
        arr_mask[sum_0 == 0] = np.nan
        arr_masked = arr + arr_mask

        # Pixel-level predicted classes
        pred_classes = arr_masked[5, :, :]
        nonwoody      = np.sum(np.isin(pred_classes, [1]))
        shrub_survival = np.sum(np.isin(pred_classes, [2]))
        shrub_mort    = np.sum(np.isin(pred_classes, [3]))
        tree_survival  = np.sum(np.isin(pred_classes, [4]))
        tree_mort      = np.sum(np.isin(pred_classes, [5]))

        # Probability extraction — masked to woody pixels only (class >= 4)
        nonwoody_mask = arr_masked[5, :, :] >= 4

        survival_prob = arr_masked[3]
        survival_prob = np.where(nonwoody_mask, survival_prob, np.nan)
        survival_prob = np.median(survival_prob[~np.isnan(survival_prob)])

        mortality_prob = arr_masked[4]
        mortality_prob = np.where(nonwoody_mask, mortality_prob, np.nan)
        mortality_prob = np.median(mortality_prob[~np.isnan(mortality_prob)])

        mort_surv_prob = survival_prob + mortality_prob
        survival_prob_norm = np.round((survival_prob / mort_surv_prob) * 100, 2)
        mortality_prob_norm = np.round((mortality_prob / mort_surv_prob) * 100, 2)

        l.append((tree_id, fire, nonwoody, tree_survival, shrub_survival, tree_mort, shrub_mort,
                  survival_prob, mortality_prob, mort_surv_prob, survival_prob_norm, mortality_prob_norm))

    print(np.round(len(l)/len(df['tree_id'].unique())*100, 1), '% complete')

result = pd.DataFrame(l, columns=[
    'tree_id', 'fire', 'nonwoody', 'tree_survival', 'shrub_survival',
    'tree_mort', 'shrub_mort', 'survival_prob', 'mortality_prob',
    'mort_surv_prob', 'survival_prob_norm', 'mortality_prob_norm'
])

result.to_csv(OUT_PATH, index=False)
print(f'Done: {len(result)} trees saved to {OUT_PATH}')