"""
01_prepare_crowns.py
--------------------
Filters STI stem point observations and matches them to LiDAR-derived crown
polygons. Trees without a crown match are assigned a 6 m buffer crown.
Trees within 15 m of the fire perimeter edge are excluded (Landsat edge
artifact).

potential-crowns.shp corresponds to Step2-potential-crownsb.shp in the
original pipeline — a spatial join of LiDAR crown polygons against filtered
STI points, compiled from private LiDAR tile directories. It is provided
in data/raw/ and still contains many-to-many relationships (multiple crowns
per point, multiple points per crown) that are resolved here.

Inputs
------
data/raw/AllTreeRecentGrowthAssessments_20230802-KNP_CASTLE.shp
data/raw/potential-crowns.shp
data/misc/laz_grid.shp
data/misc/fire-perims.shp

Output
------
data/intermediate/trees_with_crowns.shp
"""

import os
import pandas as pd
import geopandas as gpd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

STI_PATH    = 'data/raw/AllTreeRecentGrowthAssessments_20230802-KNP_CASTLE.shp'
CROWNS_PATH = 'data/raw/potential-crowns.shp'
LAZ_GRID    = 'data/misc/laz_grid.shp'
FIRE_PERIMS = 'data/misc/fire-perims.shp'
OUT_PATH    = 'data/intermediate/trees_with_crowns.shp'

os.makedirs('data/intermediate', exist_ok=True)


# ---------------------------------------------------------------------------
# Step 1: Filter STI points
# ---------------------------------------------------------------------------

pts = gpd.read_file(STI_PATH)

pts = pts[pts['tree_type'].isin(['Juvenile', 'Mature', 'Other'])]
pts = pts[pts['status'].isin(['Alive', 'Dead', 'NA'])]

# Remove trees last observed dead before 2020
pts['obs_year'] = pts['last_obser'].astype(str).str.split('-').str[0].astype(int)
pts = pts[~((pts['obs_year'] < 2020) & (pts['last_known'] == 'Dead'))]


def get_first_tree_ids(tree_ids):
    """
    Keep one record per tree stem group.
    - T: single trunk, keep
    - T1-T5: multiple trunks, keep first
    - R1-R4, RA-RM: resprouting stems, exclude
    """
    df = pd.DataFrame({'full': tree_ids})
    split_df = df['full'].str.rsplit('_', n=1, expand=True)
    df['stem'] = split_df[0]
    df['suffix'] = split_df[1]

    def classify_suffix(s):
        if s in ['T']:
            return 'T'
        elif s in ['T1', 'T2', 'T3', 'T4', 'T5']:
            return 'T1'
        elif s in ['R1', 'R2', 'R3', 'R4']:
            return None
        elif s in ['RA', 'RB', 'RC', 'RD', 'RE', 'RF', 'RG', 'RH', 'RI', 'RJ', 'RK', 'RL', 'RM']:
            return None
        else:
            return s

    df['group'] = df['suffix'].map(classify_suffix)
    df = df.dropna(subset=['group'])
    df['tree_group'] = df['stem'] + '_' + df['group']
    df_sorted = df.sort_values(['tree_group', 'suffix'])
    first_ids = df_sorted.groupby('tree_group', as_index=False).first()
    return first_ids['full'].tolist()


first_tree_ids = get_first_tree_ids(pts['tree_id'])
pts = pts[pts['tree_id'].isin(first_tree_ids)].copy()
pts = pts[~pts['tree_id'].duplicated(keep='first')]

cols = [
    'sti_id', 'tree_id', 'tree_type', 'grove', 'status',
    'last_known', 'last_obser', 'lidar_cano',
    'crown_scor', 'crown_torc', 'geometry'
]
pts = pts[cols]

print(f'Step 1: {len(pts)} STI points after filtering')


# ---------------------------------------------------------------------------
# Step 2: Resolve multiple crowns per point — pick closest crown centroid
# ---------------------------------------------------------------------------

# potential-crowns.shp has one row per crown polygon that spatially intersects
# an STI point, with STI attributes joined. A point can intersect multiple
# crowns (overlapping polygons), and a crown can intersect multiple points.

crowns = gpd.read_file(CROWNS_PATH)

tree_ids = list(pts['tree_id'].unique())
l = []
for tree_id in tree_ids:
    dft = crowns[crowns['tree_id'] == tree_id].copy()
    if len(dft) == 0:
        continue
    if len(dft) == 1:
        l.append(dft)
    else:
        # Multiple crowns intersect this point — pick closest centroid to point
        point = pts[pts['tree_id'] == tree_id]['geometry'].iloc[0]
        dft_crowns = dft.copy()
        dft_crowns['geometry'] = dft_crowns['geometry'].centroid
        dft['d'] = dft_crowns.distance(point)
        dft = dft.sort_values(by='d', ascending=True).iloc[:1].copy()
        l.append(dft)

df3 = gpd.GeoDataFrame(pd.concat(l, ignore_index=True), crs=crowns.crs)

print(f'Step 2: {len(df3)} crowns after resolving multiple crowns per point')


# ---------------------------------------------------------------------------
# Tile assignment overrides
#
# Seven trees have multiple equidistant crown candidates, making tile
# assignment non-deterministic. The treeID, f, and ns values below match
# the exact crown polygons used in the published analysis and ensure valid
# raster overlap in 02_prepare_labels.py.
# ---------------------------------------------------------------------------

# tree_id: (treeID, f, ns)
tile_overrides = {
    'GARF6808_T': (5766, '346000_4021000', 'south'),
    'ATWE4272_T': (9049, '348000_4036000', 'south'),
    'SKAG62_T':   (2945, '334000_4053000', 'north'),
    'MUIR151_T':  (680,  '335000_4055000', 'north'),
    'REMO1283_T': (3394, '330000_4062000', 'north'),
    'REMO7510_T': (3386, '328000_4063000', 'north'),
    'REMO3465_T': (8105, '327000_4063000', 'north'),
}

for tree_id, (treeID_val, f_val, ns_val) in tile_overrides.items():
    if tree_id in df3['tree_id'].values:
        correct = crowns[
            (crowns['tree_id'] == tree_id) &
            (crowns['treeID'] == treeID_val) &
            (crowns['f'] == f_val) &
            (crowns['ns'] == ns_val)
        ]
        if len(correct) > 0:
            df3 = df3[df3['tree_id'] != tree_id]
            df3 = gpd.GeoDataFrame(
                pd.concat([df3, correct.iloc[:1]], ignore_index=True), crs=crowns.crs)


# ---------------------------------------------------------------------------
# Step 3: Resolve multiple points per crown — keep tallest tree
# ---------------------------------------------------------------------------

df3['crown_uid'] = (
    df3['treeID'].astype(int).astype(str) + '-' +
    df3['f'].astype(str) + '-' +
    df3['ns'].astype(str)
)

l = []
for crown_uid in df3['crown_uid'].unique():
    dft = df3[df3['crown_uid'] == crown_uid]
    if len(dft) > 1:
        # Multiple points claim this crown — assign to the tallest tree
        max_height = dft['lidar_cano'].max()
        l.append(dft[dft['lidar_cano'] == max_height].head(1))
    else:
        l.append(dft)

crowns_clean = gpd.GeoDataFrame(pd.concat(l, ignore_index=True), crs=crowns.crs)
crowns_clean['match'] = 'matched'

print(f'Step 3: {len(crowns_clean)} unique crown-tree pairs')


# ---------------------------------------------------------------------------
# Step 4: Assign unmatched points a 6 m buffer crown
# ---------------------------------------------------------------------------

# Points with no intersecting crown polygon
pts_unmatched = pts[~pts['tree_id'].isin(crowns_clean['tree_id'].tolist())].copy()
pts_unmatched = pts_unmatched[['tree_id', 'geometry', 'sti_id', 'grove']]

# Assign each unmatched point to its nearest LAZ tile
grid = gpd.read_file(LAZ_GRID)
grid['file'] = grid['filename'].str.split('.').str[0]
grid['f'] = grid['file'].str.split('-').str[0]
grid['ns'] = grid['file'].str.split('-').str[1]
grid = grid[['geometry', 'ns', 'f', 'file']]
grid_centroids = grid.copy()
grid_centroids['geometry'] = grid_centroids['geometry'].centroid

l = []
for tree_id in pts_unmatched['tree_id'].tolist():
    pt = pts_unmatched[pts_unmatched['tree_id'] == tree_id]
    dft = gpd.sjoin(pt, grid).drop('index_right', axis=1)
    if len(dft) > 1:
        grid_centroids_tmp = grid_centroids[
            grid_centroids['file'].isin(dft['file'].tolist())].copy()
        grid_centroids_tmp['distance'] = grid_centroids_tmp.geometry.distance(
            pt.geometry.iloc[0])
        grid_centroids_tmp = grid_centroids_tmp.sort_values(
            by='distance', ascending=True).head(1)
        dft = dft[dft['file'] == grid_centroids_tmp['file'].iloc[0]]
    l.append(dft)

dfout = pd.concat(l, ignore_index=True)
dfout['geometry'] = dfout['geometry'].buffer(6)
dfout['match'] = 'unmatched'

print(f'Step 4: {len(dfout)} unmatched points assigned 6 m buffer')


# ---------------------------------------------------------------------------
# Step 5: Combine and filter by fire perimeter
# ---------------------------------------------------------------------------

df = pd.concat([crowns_clean, dfout], ignore_index=True)

# Exclude trees within 15 m of fire perimeter edge (Landsat boundary artifact)
fireperims = gpd.read_file(FIRE_PERIMS)[['FIRE_NAME', 'geometry']]
fireperims['geometry'] = fireperims['geometry'].buffer(-15)

df_centroids = df.copy()
df_centroids['geometry'] = df_centroids['geometry'].centroid
df_centroids = gpd.sjoin(df_centroids, fireperims)[['tree_id', 'FIRE_NAME']]

df = df[df['tree_id'].isin(df_centroids['tree_id'].tolist())]
df = df[['sti_id', 'tree_id', 'grove', 'geometry', 'match', 'ns', 'f']]
df = pd.merge(df, df_centroids, on='tree_id')

df.to_file(OUT_PATH)
print(f'Done: {len(df)} trees saved to {OUT_PATH}')
print(f'  matched:               {(df["match"] == "matched").sum()}')
print(f'  unmatched (6m buffer): {(df["match"] == "unmatched").sum()}')