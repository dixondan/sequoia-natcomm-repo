"""
03c_pb_upslope.py
-----------------
!! NOT REPRODUCIBLE — requires private data !!

Extracts two predictor variables for each tree:

1. Prescribed burn history: spatial join of tree centroids against the SEKI
   burn history polygon dataset. Each tree is assigned the most recent burn
   year (1990-2021, excluding the KNP and SQF wildfires). Binary flags are
   created for recent (2011-2021) and older (2001-2010) burns.

2. Upslope classification: determines whether each tree was on the upslope
   side of the fire front at the time of burning, based on fire spread
   direction and local terrain aspect. Fire spread direction requires daily
   fire perimeter data and a manual editing step (see notes below).

Note on upslope classification
------------------------------
Fire spread direction (avg_az) was computed from daily KNP and Castle fire
perimeter data, then manually reviewed and corrected for trees where the
automated method failed (e.g. no fire front found within 1 km buffer).
The manually edited file (points-avg-az-post-edit.shp) is required to
produce the final upslope classification and is not shared. The final
upslope values are included in the pre-computed output CSV.

Private inputs required
-----------------------
BURN_HISTORY_SHP : SEKI prescribed burn history polygons
ASPECT_RASTER    : terrain aspect raster
FIRE_DIR_KNP     : daily KNP fire perimeter shapefile
FIRE_DIR_CASTLE  : daily Castle fire perimeter shapefile
UPSLOPE_EDIT_SHP : manually edited fire azimuth points (post-edit)
UPSLOPE_CSV      : pre-computed upslope_120m values (tree_id-pb-upslope.csv)

Shared input
------------
data/intermediate/trees_with_crowns.shp    from 01_prepare_crowns.py

Output (pre-computed, provided in repo)
-------
data/intermediate/all-trees-pbupslope-vars.csv
"""

import os
import json
import math
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio import mask
from shapely.geometry import Point, LineString

# ---------------------------------------------------------------------------
# Private paths -- update if re-running
# ---------------------------------------------------------------------------

BURN_HISTORY_SHP = '/path/to/seki_burn_history_2023.shp'              # <-- not shared
ASPECT_RASTER    = '/path/to/north-south-aspect.tif'                   # <-- not shared
FIRE_DIR_KNP     = '/path/to/knp-daily.shp'                           # <-- not shared
FIRE_DIR_CASTLE  = '/path/to/castle-daily.shp'                        # <-- not shared
UPSLOPE_EDIT_SHP = '/path/to/points-avg-az-post-edit.shp'             # <-- not shared
UPSLOPE_CSV      = '/path/to/tree_id-pb-upslope.csv'                  # <-- not shared

CROWNS_PATH = 'data/intermediate/trees_with_crowns.shp'
OUT_PATH    = 'data/intermediate/all-trees-pbupslope-vars.csv'

os.makedirs('data/intermediate', exist_ok=True)


# ---------------------------------------------------------------------------
# Section 1: Prescribed burn history
# ---------------------------------------------------------------------------

df_burn = gpd.read_file(CROWNS_PATH)
df_burn['geometry'] = df_burn['geometry'].centroid

burn_history = gpd.read_file(BURN_HISTORY_SHP)
burn_history = burn_history.to_crs('EPSG:26911')
burn_history = burn_history[burn_history['GIS_Acres'] > 1]
burn_history = burn_history[burn_history['YEARNO'] < 2021]
burn_history = burn_history[burn_history['YEARNO'] > 1990]
burn_history = burn_history[~burn_history['FireName'].isin(['KNP Complex', 'SQF Complex'])]
burn_history['geometry'] = burn_history.geometry.buffer(30)

# Assign most recent burn year to each tree
joined = gpd.sjoin(df_burn, burn_history[['YEARNO', 'geometry']], how='left', predicate='intersects')
most_recent_burn = joined.groupby(joined.index)['YEARNO'].max()
df_burn['year_burnt'] = most_recent_burn

print(f"Burned before 2011: {len(df_burn[df_burn['year_burnt'] < 2011])}")
print(f"Burned 2011-2021:   {len(df_burn[df_burn['year_burnt'] >= 2011])}")
print(f"Never burned:       {len(df_burn[df_burn['year_burnt'].isna()])}")

# Binary burn flags
df_burn['burned_recent'] = np.where(
    (df_burn['year_burnt'] >= 2011) & (df_burn['year_burnt'] <= 2021), 1, 0)
df_burn['burned_older'] = np.where(
    (df_burn['year_burnt'] >= 2001) & (df_burn['year_burnt'] < 2011), 1, 0)
df_burn['burned_recent'] = df_burn['burned_recent'].fillna(0)
df_burn['burned_older'] = df_burn['burned_older'].fillna(0)
df_burn = df_burn[['tree_id', 'year_burnt', 'burned_recent', 'burned_older']]


# ---------------------------------------------------------------------------
# Section 2: Fire spread direction (azimuth)
# Note: output requires manual review before upslope classification
# ---------------------------------------------------------------------------

def calculate_azimuth_and_distance(coord, target_point):
    x1, y1 = coord
    x2, y2 = target_point.lon, target_point.lat
    dx = x2 - x1
    dy = y2 - y1
    angle = math.degrees(math.atan2(dy, dx))
    azimuth = (450 - angle) % 360
    distance = math.sqrt(dx ** 2 + dy ** 2)
    return azimuth, distance


def average_azimuth(azimuths):
    radians = [math.radians(a) for a in azimuths]
    x = sum(math.cos(r) for r in radians) / len(radians)
    y = sum(math.sin(r) for r in radians) / len(radians)
    avg_angle = math.atan2(y, x)
    return (math.degrees(avg_angle) + 360) % 360


def classify_direction(tree_id, points, fires, max_buffer):
    point = points[points['tree_id'] == tree_id]
    point['lon'] = point['geometry'].x
    point['lat'] = point['geometry'].y

    fire_clip = fires.clip(point)
    burn_date = fire_clip[fire_clip.contains(point.iloc[0].geometry)]['time'].min()
    burn_date_index = fires[fires['time'] == burn_date].index[0]
    pre_date = fires.iloc[burn_date_index - 1]

    point_buffer = point.copy()
    point_buffer['geometry'] = point_buffer['geometry'].buffer(max_buffer)
    pre_fire = fires[fires['time'] == pre_date['time']]
    pre_fire['geometry'] = pre_fire['geometry'].boundary
    fire_fire_clip = pre_fire.clip(point_buffer)

    if len(fire_fire_clip) == 0:
        print('no front found')
    else:
        length = fire_fire_clip['geometry'].length.iloc[0]
        spacing = 30
        distances = list(range(0, int(length), spacing))
        distances.append(int(length))
        new_coords = [fire_fire_clip['geometry'].iloc[0].interpolate(d).coords[0] for d in distances]
        new_fire_front = LineString(new_coords)

        data = []
        for coord in new_fire_front.coords:
            azimuth, distance = calculate_azimuth_and_distance(coord, point)
            data.append({'geometry': Point(coord), 'azimuth': azimuth, 'distance': distance})

        result = gpd.GeoDataFrame(data, crs='EPSG:26911')
        result['tree_id'] = tree_id

        average = average_azimuth(result['azimuth'].tolist())
        if len(result) > 10:
            return average


# NOTE: the fire azimuth output requires manual review and correction before
# upslope classification can proceed. The manually edited file is
# UPSLOPE_EDIT_SHP (points-avg-az-post-edit.shp).

points = gpd.read_file(CROWNS_PATH)
points['geometry'] = points['geometry'].centroid
points = points[['tree_id', 'geometry']]

knp = gpd.read_file(FIRE_DIR_KNP)
castle = gpd.read_file(FIRE_DIR_CASTLE)
fires = pd.concat([knp, castle])
fires['time'] = pd.to_datetime(fires['t'])
fires = fires.sort_values(by='time').reset_index(drop=True)

l = []
for tree_id in points['tree_id'].tolist():
    print(tree_id)
    r = classify_direction(tree_id, points, fires, max_buffer=1000)
    l.append((tree_id, r))

res_out = pd.DataFrame(l, columns=['tree_id', 'avg_az'])
# Save for manual review before proceeding to upslope classification
# res_out.to_file('/path/to/points-avg-az-pre-edit.shp')


# ---------------------------------------------------------------------------
# Section 3: Upslope classification from aspect raster
# (requires manually edited azimuth file)
# ---------------------------------------------------------------------------

df = gpd.read_file(UPSLOPE_EDIT_SHP)
tree_ids = df['tree_id'].tolist()

l = []
for tree_id in tree_ids:
    dft = df[df['tree_id'] == tree_id]
    for buffer_size in [10, 50, 120]:
        dft['geometry'] = dft['geometry'].buffer(buffer_size)
        crown_json = dft.to_json()
        geom2clip = [feature['geometry'] for feature in json.loads(crown_json)['features']]

        with rasterio.open(ASPECT_RASTER) as src:
            arr, trans = mask.mask(src, geom2clip, crop=True)
            arr = arr[arr != -9999]
            radians = [math.radians(a) for a in arr.reshape(-1, 1)[:, 0]]
            x = sum(math.cos(r) for r in radians) / len(radians)
            y = sum(math.sin(r) for r in radians) / len(radians)
            avg_angle = math.atan2(y, x)
            avg_azimuth = (math.degrees(avg_angle) + 360) % 360
            l.append((tree_id, buffer_size, avg_azimuth))
            print(len(l))

terrain = pd.DataFrame(l, columns=['tree_id', 'buffer_size', 'avg_azimuth'])
terrain = terrain.pivot(index='tree_id', columns='buffer_size', values='avg_azimuth').reset_index()
terrain.columns = ['tree_id', '10m', '50m', '120m']

direction = gpd.read_file(UPSLOPE_EDIT_SHP)


def classify_upslope_fire(row):
    fire_min = (row['avg_az-f'] - 60) % 360
    fire_max = (row['avg_az-f'] + 60) % 360
    aspect_min = (row['120m'] - 60) % 360
    aspect_max = (row['120m'] + 60) % 360

    def ranges_overlap(r1_min, r1_max, r2_min, r2_max):
        range1 = [(r1_min, r1_max)] if r1_min <= r1_max else [(r1_min, 360), (0, r1_max)]
        range2 = [(r2_min, r2_max)] if r2_min <= r2_max else [(r2_min, 360), (0, r2_max)]
        for r1_start, r1_end in range1:
            for r2_start, r2_end in range2:
                if not (r1_end < r2_start or r2_end < r1_start):
                    return True
        return False

    overlap = ranges_overlap(fire_min, fire_max, aspect_min, aspect_max)
    return 0 if overlap else 1  # 1 = upslope fire, 0 = not upslope


df_upslope = pd.merge(direction, terrain, on='tree_id')
df_upslope['upslope'] = df_upslope.apply(classify_upslope_fire, axis=1)


# ---------------------------------------------------------------------------
# Section 4: Compile pb + upslope into final output
# ---------------------------------------------------------------------------

# upslope_120m was manually reviewed — load from pre-computed CSV
upslope = pd.read_csv(UPSLOPE_CSV)
upslope = upslope[['tree_id', 'upslope_120m']]

df = gpd.read_file(CROWNS_PATH)
df = df[['tree_id', 'grove']]
df = pd.merge(df, upslope, on='tree_id', how='left')
df = pd.merge(df, df_burn, on='tree_id')

df.to_csv(OUT_PATH, index=False)
print(f'Done: {len(df)} trees saved to {OUT_PATH}')