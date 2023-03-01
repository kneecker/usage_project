import os, gc
import rasterio as rio
import numpy as np
import geopandas as gpd
from numpy import ma
from rasterio.mask import mask as rio_mask
from rasterio.features import shapes as rio_shapes
from rasterio.features import rasterize
from shapely.geometry import shape as shapely_shape


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
OSM_CLASSES = {
    'water': 0, 'trees': 1, 'grass': 2, 'wetland': 3, 'crops': 4,
    'shrub': 5, 'built': 6, 'bareland': 7, 'snowice': 8,
    'vineyard': 32
}

REGIONS_PATH = 'work_rus_regions_2023_with_buffer_1km.gpkg'
REGIONS_NAME_PROPERTY = 'name'
region_name = 'Краснодарский край'
years = [2020, 2021]
regions_gdf = gpd.read_file(REGIONS_PATH).to_crs('EPSG:3857')
region_geometry = regions_gdf.loc[regions_gdf[REGIONS_NAME_PROPERTY] == region_name].iloc[0].geometry

work_path = os.path.join('/home/user/usage_project/data', region_name)
save_path = os.path.join(work_path, 'new_masks')
if not os.path.exists(work_path):
    os.makedirs(work_path)
for year in years:
    esa_path = os.path.join(work_path, f'esa/esa_{year}.tif')
    dw_path = os.path.join(work_path, f'dw/dw_{year}.tif')
    esri_path = os.path.join(work_path, f'esri/esri_{year}.tif')
    year_save_path = os.path.join(save_path, f'{year}')
    if not os.path.exists(year_save_path):
        os.makedirs(year_save_path)
    heights = []
    widths = []
    with rio.open(dw_path, 'r') as rds:
        raster_profile = rds.profile.copy()
        dw_raw_data, raster_transform = rio_mask(rds, [region_geometry], crop=True, nodata=255)
    dw_raw_data = dw_raw_data[0].astype(np.uint8)
    dw_data = np.full(dw_raw_data.shape, 100, np.uint8)
    for c, v in DW_CLASSES.items():
        dw_data[dw_raw_data == v] = CLASS_VALUES[c]
    del dw_raw_data
    gc.collect()
    heights.append(dw_data.shape[0])
    widths.append(dw_data.shape[1])

    with rio.open(esri_path, 'r') as rds:
        esri_raw_data, _ = rio_mask(rds, [region_geometry], crop=True, nodata=255)
    esri_raw_data = esri_raw_data[0].astype(np.uint8)
    esri_data = np.full(esri_raw_data.shape, 100, np.uint8)
    for c, v in ESRI_CLASSES.items():
        esri_data[esri_raw_data == v] = CLASS_VALUES[c]
    del esri_raw_data
    gc.collect()
    heights.append(esri_data.shape[0])
    widths.append(esri_data.shape[1])

    with rio.open(esa_path, 'r') as rds:
        esa_raw_data, _ = rio_mask(rds, [region_geometry], crop=True, nodata=255)
    esa_raw_data = esa_raw_data[0].astype(np.uint8)
    esa_data = np.full(esa_raw_data.shape, 100, np.uint8)
    for c, v in ESA_CLASSES.items():
        esa_data[esa_raw_data == v] = CLASS_VALUES[c]
    del esa_raw_data
    gc.collect()
    heights.append(esa_data.shape[0])
    widths.append(esa_data.shape[1])

    min_height = min(heights)
    min_width = min(widths)
    raster_profile.update({
        'transform': raster_transform,
        'width': min_width,
        'height': min_height,
        'compress': 'lzw'
    })
    dw_data = dw_data[:min_height, :min_width]
    esri_data = esri_data[:min_height, :min_width]
    esa_data = esa_data[:min_height, :min_width]

    intersect_all = np.where((esa_data == CLASS_VALUES['crops']) & (dw_data == CLASS_VALUES['crops']) &
                             (esri_data == CLASS_VALUES['crops']),
                             dw_data, np.full(dw_data.shape, 100, np.uint8))
    with rio.open(os.path.join(year_save_path, f'crops_intersect_all.tif'), 'w', **raster_profile) as wds:
        wds.write(intersect_all, 1)
    # del intersect_all
    # gc.collect()
    print('intersect all')

    intersect_esa_dw = np.where(
        (esa_data == CLASS_VALUES['crops']) & (dw_data == CLASS_VALUES['crops']),
        esa_data, np.full(dw_data.shape, 100, np.uint8))
    with rio.open(os.path.join(year_save_path, f'crops_intersect_esa_dw.tif'), 'w', **raster_profile) as wds:
        wds.write(intersect_esa_dw, 1)
    # del intersect_esa_dw
    # gc.collect()
    print('intersect esa dw')

    intersect_esa_esri = np.where(
        (esa_data == CLASS_VALUES['crops']) & (esri_data == CLASS_VALUES['crops']),
        esa_data, np.full(dw_data.shape, 100, np.uint8))
    with rio.open(os.path.join(year_save_path, f'crops_intersect_esa_esri.tif'), 'w', **raster_profile) as wds:
        wds.write(intersect_esa_esri, 1)
    # del intersect_esa_esri
    # gc.collect()
    print('intersect esa esri')
    if 'esa_data' in globals():
        del esa_data
    if 'dw_data' in globals():
        del dw_data
    if 'esri_data' in globals():
        del esri_data
    gc.collect()
    year_crops = np.where((intersect_all == CLASS_VALUES['crops']) | (intersect_esa_dw == CLASS_VALUES['crops']) |
                          (intersect_esa_esri == CLASS_VALUES['crops']),
                          np.full(intersect_all.shape, CLASS_VALUES['crops'], np.uint8),
                          np.full(intersect_all.shape, 100, np.uint8))
    if 'full_crops' in globals():
        if year_crops.shape[0] != full_crops.shape[0]:
            min_height = min(year_crops.shape[0], full_crops.shape[0])
        if year_crops.shape[1] != full_crops.shape[1]:
            min_width = min(year_crops.shape[1], full_crops.shape[1])
        full_crops = full_crops[:min_height, :min_width]
        year_crops = year_crops[:min_height, :min_width]
        full_crops[year_crops == CLASS_VALUES['crops']] = CLASS_VALUES['crops']
    else:
        full_crops = year_crops.copy()
    with rio.open(os.path.join(save_path, f'1_crops_all_{year}.tif'), 'w', **raster_profile) as wds:
        wds.write(year_crops, 1)
    del year_crops
    gc.collect()
with rio.open(os.path.join(save_path, f'2_intersect_2020_2021.tif'), 'w', **raster_profile) as wds:
    wds.write(full_crops, 1)

esa_path = os.path.join(work_path, f'esa/esa_2021.tif')
with rio.open(esa_path, 'r') as rds:
    esa_raw_data, _ = rio_mask(rds, [region_geometry], crop=True, nodata=255)
esa_raw_data = esa_raw_data[0, :min_height, :min_width].astype(np.uint8)
esa_data = np.full(esa_raw_data.shape, 100, np.uint8)
for c, v in ESA_CLASSES.items():
    esa_data[esa_raw_data == v] = CLASS_VALUES[c]
del esa_raw_data
gc.collect()
after_esa_trees = np.where(esa_data == CLASS_VALUES['trees'],
                           esa_data, full_crops)
with rio.open(os.path.join(save_path, f'3_after_esa_trees.tif'), 'w', **raster_profile) as wds:
    wds.write(after_esa_trees, 1)
del esa_data, full_crops
gc.collect()


dw_path = os.path.join(work_path, f'dw/dw_2022.tif')
with rio.open(dw_path, 'r') as rds:
    dw_raw_data, _ = rio_mask(rds, [region_geometry], crop=True, nodata=255)
dw_raw_data = dw_raw_data[0, :min_height, :min_width].astype(np.uint8)
dw_data = np.full(dw_raw_data.shape, 100, np.uint8)
for c, v in DW_CLASSES.items():
    dw_data[dw_raw_data == v] = CLASS_VALUES[c]
del dw_raw_data
gc.collect()
after_dw_built = np.where(dw_data == CLASS_VALUES['built'],
                          dw_data, after_esa_trees)
with rio.open(os.path.join(save_path, f'4_after_dw_built.tif'), 'w', **raster_profile) as wds:
    wds.write(after_dw_built, 1)
del dw_data, after_esa_trees
gc.collect()

osm_path = os.path.join(work_path, f'osm_new/osm_mask_no_settlement.tif')  # osm_new/osm_mask_no_settlement.tif
with rio.open(osm_path, 'r') as rds:
    osm_raw_data, _ = rio_mask(rds, [region_geometry], crop=True, nodata=255)
osm_raw_data = osm_raw_data[0, :min_height, :min_width].astype(np.uint8)
osm_data = np.full(osm_raw_data.shape, 100, np.uint8)
for c, v in OSM_CLASSES.items():
    if c in CLASS_VALUES.keys():
        osm_data[osm_raw_data == v] = CLASS_VALUES[c]
    else:
        osm_data[osm_raw_data == v] = v
del osm_raw_data
gc.collect()
after_osm = np.where(~np.isin(osm_data, (100, 255)),
                     osm_data, after_dw_built)
with rio.open(os.path.join(save_path, f'5_after_osm.tif'), 'w', **raster_profile) as wds:
    wds.write(after_osm, 1)
del osm_data, after_dw_built
gc.collect()

max_areas = [1, 2]
small_polys = {m: [] for m in max_areas}
for shp, value in rio_shapes(after_osm, after_osm == CLASS_VALUES['crops'],  # ~np.isin(after_osm, [100, 255])
                             connectivity=4, transform=raster_transform):
    poly = shapely_shape(shp)
    for max_area in max_areas:
        if poly.area * 0.0001 < max_area:
            small_polys[max_area].append(poly)
for max_area in max_areas:
    small_mask = rasterize(small_polys[max_area], out_shape=after_osm.shape, fill=255,
                           default_value=200, transform=raster_transform, dtype=np.uint8)
    raster_arr = np.where(small_mask == 200, small_mask, after_osm)
    masked_arr = ma.masked_array(raster_arr, raster_arr == 200)
    iter_num = 0
    prev_count = masked_arr.mask.sum()
    while np.any(masked_arr.mask):
        for shift in (-1, 1):
            for axis in (0, 1):
                shifted_arr = np.roll(masked_arr, shift=shift, axis=axis)
                shifted_arr[shifted_arr == 255] = ma.masked
                idx = ~shifted_arr.mask * masked_arr.mask
                masked_arr[idx] = shifted_arr[idx]
        cur_count = masked_arr.mask.sum()
        if prev_count <= cur_count:
            break
        else:
            prev_count = cur_count
        print('iter', iter_num, cur_count)
        iter_num += 1
    raster_arr = ma.getdata(masked_arr)
    raster_arr[raster_arr == 200] = 100
    with rio.open(os.path.join(save_path, f'6_after_osm_no_small_{max_area}ha.tif'), 'w', **raster_profile) as wds:
        wds.write(raster_arr, 1)
    del small_mask, raster_arr, masked_arr
