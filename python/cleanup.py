"""
cleanup_sti_columns.py
----------------------
NOT FOR SHARING — run locally before publishing the data repo.

Removes unused columns from the STI shapefile and updates potential-crowns.shp
to match. Only columns used in the analysis pipeline are kept.

Inputs
------
data/raw/AllTreeRecentGrowthAssessments_20230802-KNP_CASTLE.shp
data/raw/potential-crowns.shp

Outputs (overwrites in place)
------------------------------
data/raw/AllTreeRecentGrowthAssessments_20230802-KNP_CASTLE.shp
data/raw/potential-crowns.shp
"""

import geopandas as gpd

STI_PATH    = 'data/raw/AllTreeRecentGrowthAssessments_20230802-KNP_CASTLE.shp'
CROWNS_PATH = 'data/raw/potential-crowns.shp'

# Columns used in the analysis pipeline:
#   01_prepare_crowns.py : tree_id, tree_type, status, last_obser, last_known,
#                          lidar_cano, grove, sti_id, geometry
#   04_compile_model_data.py : tree_id, last_obser, status, crown_scor, crown_torc
STI_KEEP = [
    'tree_id', 'sti_id', 'tree_type', 'grove', 'status',
    'last_known', 'last_obser', 'lidar_cano',
    'crown_scor', 'crown_torc', 'geometry'
]

# ---------------------------------------------------------------------------
# Clean STI shapefile
# ---------------------------------------------------------------------------

sti = gpd.read_file(STI_PATH)
print(f'STI before: {len(sti.columns)} columns')

missing = [c for c in STI_KEEP if c not in sti.columns]
if missing:
    raise ValueError(f'Expected columns not found in STI: {missing}')

sti = sti[STI_KEEP]
sti.to_file(STI_PATH)
print(f'STI after:  {len(sti.columns)} columns')


# ---------------------------------------------------------------------------
# Clean potential-crowns shapefile
# (has STI columns joined in — keep only what's needed downstream)
# ---------------------------------------------------------------------------

crowns = gpd.read_file(CROWNS_PATH)
print(f'\nCrowns before: {len(crowns.columns)} columns')

# Columns needed from crowns in 01_prepare_crowns.py
CROWNS_KEEP = [
    'treeID', 'zmax', 'zq95', 'f', 'ns',
    'tree_id', 'lidar_cano', 'grove', 'sti_id', 'geometry'
]

missing = [c for c in CROWNS_KEEP if c not in crowns.columns]
if missing:
    print(f'Warning: columns not found in crowns (will skip): {missing}')
    CROWNS_KEEP = [c for c in CROWNS_KEEP if c in crowns.columns]

crowns = crowns[CROWNS_KEEP]
crowns.to_file(CROWNS_PATH)
print(f'Crowns after: {len(crowns.columns)} columns')

print('\nDone.')