import os, gc
import rasterio as rio
from rasterio.mask import mask
import numpy as np
import geopandas as gpd


CLASSES = [
    'water', 'trees', 'grass', 'wetland', 'crops',
    'shrub', 'built', 'bareland', 'snowice']
CLASS_VALUES = {c: i for i, c in enumerate(CLASSES)}
DW_CLASSES = CLASS_VALUES
ESA_CLASSES = {
    'water': 80, 'trees': 10, 'grass': 30, 'wetland': 90, 'crops': 40,
    'shrub': 20, 'built': 50, 'bareland': 60, 'snowice': 70
}
ESRI_CLASSES = {
    'water': 1, 'trees': 2, 'grass': 11, 'wetland': 4, 'crops': 5,
    'shrub': 6, 'built': 7, 'bareland': 8, 'snowice': 9
}

# OSM_PATH = '/home/user/osm_rus_2022/RU-60797a9e-20221123-ru-shape/data'
# REGIONS_PATH = os.path.join(OSM_PATH, 'boundary-polygon.shp')
# REGIONS_NAME_PROPERTY = 'NAME_RU'
REGIONS_PATH = 'work_rus_regions_2023_with_buffer_1km.gpkg'
REGIONS_NAME_PROPERTY = 'name'
region_name = 'Краснодарский край'
year = 2021
intersect_classes = ['crops', 'trees']
regions_gdf = gpd.read_file(REGIONS_PATH).to_crs('EPSG:3857')
region_geometry = regions_gdf.loc[regions_gdf[REGIONS_NAME_PROPERTY] == region_name].iloc[0].geometry

work_path = os.path.join('/home/user/usage_project/data', region_name)
save_path = os.path.join(work_path, f'test_masks/{year}')

esa_path = os.path.join(work_path, f'esa/esa_{year}.tif')
dw_path = os.path.join(work_path, f'dw/dw_{year}.tif')
osm_path = os.path.join(work_path, f'osm/osm_mask.tif')
if not os.path.exists(save_path):
    os.makedirs(save_path)

with rio.open(esa_path, 'r') as rds:
    raster_profile = rds.profile.copy()
    esa_raw_data, raster_transform = mask(rds, [region_geometry], crop=True, nodata=255)
esa_raw_data = esa_raw_data[0].astype(np.uint8)
esa_data = np.full(esa_raw_data.shape, 100, np.uint8)
for c, v in ESA_CLASSES.items():
    esa_data[esa_raw_data == v] = CLASS_VALUES[c]
del esa_raw_data
gc.collect()

with rio.open(dw_path, 'r') as rds:
    raster_profile = rds.profile.copy()
    dw_raw_data, raster_transform = mask(rds, [region_geometry], crop=True, nodata=255)
dw_raw_data = dw_raw_data[0].astype(np.uint8)
dw_data = np.full(dw_raw_data.shape, 100, np.uint8)
for c, v in DW_CLASSES.items():
    dw_data[dw_raw_data == v] = CLASS_VALUES[c]
del dw_raw_data
gc.collect()

with rio.open(osm_path, 'r') as rds:
    osm_raw_data, _ = mask(rds, [region_geometry], crop=True, nodata=255)
osm_raw_data = osm_raw_data[0].astype(np.uint8)
osm_data = np.full(osm_raw_data.shape, 100, np.uint8)
for c, v in DW_CLASSES.items():
    osm_data[osm_raw_data == v] = CLASS_VALUES[c]
del osm_raw_data
gc.collect()

min_height = min(esa_data.shape[0], dw_data.shape[0], osm_data.shape[0])
min_width = min(esa_data.shape[1], dw_data.shape[1], osm_data.shape[1])
raster_profile.update({
    'transform': raster_transform,
    'width': min_width,
    'height': min_height,
    'compress': 'lzw'
})
esa_data = esa_data[:min_height, :min_width]
dw_data = dw_data[:min_height, :min_width]
osm_data = osm_data[:min_height, :min_width]

esa_trees = np.where(
    (esa_data == CLASS_VALUES['trees']),
    esa_data, np.full(dw_data.shape, 100, np.uint8))
built_over_esa_trees = np.where(
    (osm_data == CLASS_VALUES['built']) | (dw_data == CLASS_VALUES['built']),
    np.full(dw_data.shape, CLASS_VALUES['built'], np.uint8), esa_trees)
with rio.open(os.path.join(save_path, f'built_over_esa_trees.tif'), 'w', **raster_profile) as wds:
    wds.write(built_over_esa_trees, 1)
del built_over_esa_trees
gc.collect()
print('built over esa trees')
