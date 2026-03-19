"""
04_compile_model_data.py
------------------------
Joins all predictor variables onto the crown dataset and compiles the final
modeling dataset. Applies filtering criteria and generates spatial grid
subsamples used in the ensemble model.

Two datasets are produced:
- Q1+2: all trees > 20 m in 19 groves (used for scenario projections)
- Q1:   matched crowns only, subset of Q1+2 (used for model fitting)

Inputs (all in data/intermediate/ or data/raw/ or data/misc/)
--------------------------------------------------------------
data/intermediate/trees_with_crowns.shp
data/intermediate/tree-mortality-outcomes.csv
data/intermediate/all-trees-lidar-vars.csv
data/intermediate/all-trees-topo-vars-30.csv
data/intermediate/all-trees-topo-vars-60.csv
data/intermediate/all-trees-topo-vars-90.csv
data/intermediate/all-trees-topo-vars-120.csv
data/intermediate/all-trees-pbupslope-vars.csv
data/raw/AllTreeRecentGrowthAssessments_20230802-KNP_CASTLE.shp
data/misc/fire-perims.shp
data/misc/sampling_grids/grid_90x90.shp

Outputs
-------
data/intermediate/model-trees-Q12-sampled_75m-90m-120m-wgt-median.csv
data/intermediate/grid_samples/samples90m-iter{0..500}.csv
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

CROWNS_PATH  = 'data/intermediate/trees_with_crowns.shp'
OUTCOMES     = 'data/intermediate/tree-mortality-outcomes.csv'
LIDAR_VARS   = 'data/intermediate/all-trees-lidar-vars.csv'
TOPO_30      = 'data/intermediate/all-trees-topo-vars-30.csv'
TOPO_60      = 'data/intermediate/all-trees-topo-vars-60.csv'
TOPO_90      = 'data/intermediate/all-trees-topo-vars-90.csv'
TOPO_120     = 'data/intermediate/all-trees-topo-vars-120.csv'
PB_UPSLOPE   = 'data/intermediate/all-trees-pbupslope-vars.csv'
STI_PATH     = 'data/raw/AllTreeRecentGrowthAssessments_20230802-KNP_CASTLE.shp'
FIRE_PERIMS  = 'data/misc/fire-perims.shp'

OUT_DIR      = 'data/intermediate'
GRID_DIR     = 'data/intermediate/grid_samples'
os.makedirs(GRID_DIR, exist_ok=True)


# ---------------------------------------------------------------------------
# Step 1: Load crowns, remove sparse groves, join LiDAR vars
# ---------------------------------------------------------------------------

df = gpd.read_file(CROWNS_PATH)

# Remove groves with too few trees for modeling
df = df[~df['grove'].isin(['Big Springs', 'Squirrel Creek', 'Putnam-Francis Tree', 'East Fork'])]
print(f'After removing sparse groves: {len(df)} trees')

lidarvars = pd.read_csv(LIDAR_VARS)
lidarvars = lidarvars[lidarvars['zq95'] >= 20]  # Remove trees under 20 m
lidarvars['ladder1'] = lidarvars['metrics2'] / (
    lidarvars['zlt4'] + lidarvars['metrics2'] + lidarvars['metrics3'])
lidarvars = lidarvars[[
    'tree_id', 'rumple_0_30', 'cdensity_0_30', 'rumple_0_60', 'cdensity_0_60',
    'rumple_0_90', 'cdensity_0_90', 'rumple_0_120', 'cdensity_0_120', 'ladder1', 'zq95'
]]

df = pd.merge(df, lidarvars, on='tree_id')
df = df.drop(['sti_id', 'f', 'ns'], axis=1)
print(f'After height filter (>= 20 m): {len(df)} trees')


# ---------------------------------------------------------------------------
# Step 2: Join prescribed burn history and upslope
# ---------------------------------------------------------------------------

df_pb = pd.read_csv(PB_UPSLOPE)
df_pb = df_pb[['tree_id', 'upslope_120m', 'year_burnt', 'burned_recent', 'burned_older']]

df = pd.merge(df, df_pb, on='tree_id')


# ---------------------------------------------------------------------------
# Step 3: Join topographic variables
# ---------------------------------------------------------------------------

df_topo_30 = pd.read_csv(TOPO_30)
df_topo_30.columns = ['tree_id', 'solar_30m', 'tpi_30m', 'slope_30m']

df_topo_60 = pd.read_csv(TOPO_60)
df_topo_60.columns = ['tree_id', 'solar_60m', 'tpi_60m', 'slope_60m']

df_topo_90 = pd.read_csv(TOPO_90)
df_topo_90.columns = ['tree_id', 'solar_90m', 'tpi_90m', 'slope_90m']

df_topo_120 = pd.read_csv(TOPO_120)
df_topo_120.columns = ['tree_id', 'solar_120m', 'tpi_120m', 'slope_120m']

df_topo = pd.merge(df_topo_30, df_topo_60, on='tree_id')
df_topo = pd.merge(df_topo, df_topo_90, on='tree_id')
df_topo = pd.merge(df_topo, df_topo_120, on='tree_id')

df = pd.merge(df, df_topo, on='tree_id')


# ---------------------------------------------------------------------------
# Step 4: Add grove ID
# ---------------------------------------------------------------------------

grove_id_dict = {
    'Atwell': 1, 'Board Camp': 2, 'Castle Creek': 3, 'Cedar Flat': 4,
    'Dennison': 5, 'Devils Canyon': 6, 'Garfield': 7, 'Giant Forest': 8,
    'Homers Nose': 9, 'Lost': 10, 'Muir': 11, 'New Oriole Lake': 12,
    'Oriole Lake': 13, 'Pine Ridge': 14, 'Redwood Creek': 15,
    'Redwood Mountain': 16, 'Skagway': 17, 'South Fork': 18, 'Suwanee': 19
}
df['grove_id'] = df['grove'].map(grove_id_dict)

print(f'Q1+2 dataset (> 20 m, 19 groves): {len(df)} trees')


# ---------------------------------------------------------------------------
# Step 5: Join mortality probability (response variable)
# ---------------------------------------------------------------------------

pred = pd.read_csv(OUTCOMES)
pred['total'] = pred['nonwoody'] + pred['tree_survival'] + pred['tree_mort']
pred['nonwoody_prop'] = pred['nonwoody'] / pred['total']
pred['unidentified'] = 0
pred.loc[pred['nonwoody_prop'] > 0.9, 'unidentified'] = 1  # >90% crown non-woody
pred.loc[pred['total'] == 0, 'unidentified'] = 1
pred = pred[pred['unidentified'] == 0]
pred['mortality_prob'] = pred['mortality_prob'] / 100
pred['survival_prob'] = pred['survival_prob'] / 100

# Filter to trees >= 20 m
lidarvars_zq95 = pd.read_csv(LIDAR_VARS)[['tree_id', 'zq95']]
pred = pd.merge(pred, lidarvars_zq95, on='tree_id')
pred = pred[pred['zq95'] >= 20]
pred = pred[['tree_id', 'unidentified', 'mortality_prob', 'survival_prob']]

df = pd.merge(df, pred, on='tree_id')


# ---------------------------------------------------------------------------
# Step 6: Join ground truth field observations (dead/alive labels)
# ---------------------------------------------------------------------------

pts = gpd.read_file(STI_PATH)
pts['year'] = pts['last_obser'].astype(str).str[0:4]
pts = pts[pts['year'].isin(['2021', '2022', '2023'])]
pts['year'] = pts['year'].astype(int)

s_dict = {'0%': 0, '1-10%': 10, '10-25%': 25, '25-50%': 50,
          '50-75%': 75, '75-90%': 90, '90-99%': 99, '100%': 100}
pts = pts[pts['status'].isin(['Dead', 'Alive'])]
pts['scorch_s'] = pts['crown_scor'].map(s_dict)
pts['torch_s'] = pts['crown_torc'].map(s_dict)
pts['total_scor'] = pts['scorch_s'] + pts['torch_s']
pts.loc[pts['status'] == 'Dead', 'dead'] = 1
pts.loc[pts['status'] == 'Alive', 'dead'] = 0
pts['dead_remap'] = 0
pts.loc[pts['total_scor'] >= 90, 'dead_remap'] = 1
pts = pts.dropna(subset=['dead'])
pts = pts[['tree_id', 'dead_remap', 'dead']]

# Two trees (MUIR_10020_T, MUIR_10042_T) have duplicate rows in the STI
# file with identical year and status. Deduplicate before merging.
pts = pts.drop_duplicates('tree_id', keep='first')

df = df[df['unidentified'] == 0]
df12_models = df.copy()
df12_models_pts = pd.merge(df12_models, pts, on='tree_id', how='left')

# These four trees are assigned crowns with no raster overlap in the
# original pipeline but pass all filters in the new pipeline due to
# different crown contest outcomes. Exclude to match published dataset.
exclude_ids = {'GARF7464_T', 'REMO9643_T', 'REMO1099_T', 'REMO8306_T'}
df12_models_pts = df12_models_pts[~df12_models_pts['tree_id'].isin(exclude_ids)]

df12_models_pts.to_csv(
    f'{OUT_DIR}/model-trees-Q12-sampled_75m-90m-120m-wgt-median.csv', index=False)
print(f'Q1+2 model dataset saved: {len(df12_models_pts)} rows, '
      f'{len(df12_models_pts["tree_id"].unique())} unique trees')


# ---------------------------------------------------------------------------
# Step 7: Q1 subset — matched crowns only
# ---------------------------------------------------------------------------

# Q1 = trees with a matched LiDAR crown (not 6 m buffer)
trees_q1 = gpd.read_file(CROWNS_PATH)
trees_q1 = trees_q1[trees_q1['match'] == 'matched']

dfq1 = df12_models_pts[df12_models_pts['tree_id'].isin(trees_q1['tree_id'].tolist())]
dfq1 = dfq1[dfq1['unidentified'] == 0]
print(f'Q1 dataset (matched crowns only): {len(dfq1)} trees')


# ---------------------------------------------------------------------------
# Step 8: Spatial grid subsampling — 500 iterations, 90 m grid
# ---------------------------------------------------------------------------

# Get Q1 centroids for grid sampling
centroids = dfq1.copy()
centroids['geometry'] = centroids['geometry'].centroid
centroids = centroids[['tree_id', 'geometry']]

# Remove trees within 30 m of fire perimeter edge for modeling
# (these trees are still included in scenario projections)
fire_perims = gpd.read_file(FIRE_PERIMS)
fire_perims = fire_perims.to_crs('EPSG:26911')
fire_perims['geometry'] = fire_perims['geometry'].buffer(-30)
intersect_mask = centroids.geometry.apply(lambda x: fire_perims.intersects(x).any())
centroids = centroids[intersect_mask]
print(f'Q1 centroids for grid sampling (after edge removal): {len(centroids)}')

for iteration in range(0, 501):
    for size in [90]:
        outsampled = f'{GRID_DIR}/samples{size}m-iter{iteration}.csv'
        if not os.path.exists(outsampled):

            grid = gpd.read_file(f'data/misc/grid_{size}x{size}.shp')
            grid_w_points = gpd.sjoin(grid, centroids)
            grid_p_dict = dict(zip(grid_w_points['tree_id'], grid_w_points['id']))
            centroids[f'grid_{size}'] = centroids['tree_id'].map(grid_p_dict)
            grid_keep = list(centroids[f'grid_{size}'].unique())

            l = []
            for grid_id in grid_keep:
                dft = centroids[centroids[f'grid_{size}'] == grid_id]
                if len(dft) == 1:
                    tree_id = dft['tree_id'].iloc[0]
                else:
                    tree_id = dft.sample(1, random_state=iteration)['tree_id'].iloc[0]
                df_sampled = df12_models_pts[df12_models_pts['tree_id'] == tree_id].copy()
                df_sampled[f'grid_{size}'] = size
                l.append(df_sampled)

            df2 = pd.concat(l)
            df2.to_csv(outsampled, index=False)

    #print(f'Iteration {iteration} complete')

print('Done.')

old = pd.read_csv('data/intermediate/model-trees-Q12-sampled_75m-90m-120m-wgt-median-old.csv')
new = pd.read_csv('data/intermediate/model-trees-Q12-sampled_75m-90m-120m-wgt-median.csv')

extra = set(new['tree_id']) - set(old['tree_id'])
missing = set(old['tree_id']) - set(new['tree_id'])
print(f'extra ({len(extra)}): {extra}')
print(f'missing ({len(missing)}): {missing}')