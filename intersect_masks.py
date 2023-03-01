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
region_name = 'Новгородская область'
years = [2020, 2022]
intersect_classes = ['crops', 'trees']
regions_gdf = gpd.read_file(REGIONS_PATH).to_crs('EPSG:3857')
region_geometry = regions_gdf.loc[regions_gdf[REGIONS_NAME_PROPERTY] == region_name].iloc[0].geometry

work_path = os.path.join('/home/user/usage_project/data', region_name)
for year in years:
    esa_path = os.path.join(work_path, f'esa/esa_{year}.tif')
    dw_path = os.path.join(work_path, f'dw/dw_{year}.tif')
    esri_path = os.path.join(work_path, f'esri/esri_{year}.tif')
    save_path = os.path.join(work_path, f'test_masks/{year}')
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    if not (2023 > year > 2016):
        exit()
    heights = []
    widths = []
    if year > 2014:
        with rio.open(dw_path, 'r') as rds:
            raster_profile = rds.profile.copy()
            dw_raw_data, raster_transform = mask(rds, [region_geometry], crop=True, nodata=255)
        dw_raw_data = dw_raw_data[0].astype(np.uint8)
        dw_data = np.full(dw_raw_data.shape, 100, np.uint8)
        for c, v in DW_CLASSES.items():
            dw_data[dw_raw_data == v] = CLASS_VALUES[c]
        del dw_raw_data
        gc.collect()
        heights.append(dw_data.shape[0])
        widths.append(dw_data.shape[1])

    if year > 2016:
        with rio.open(esri_path, 'r') as rds:
            esri_raw_data, _ = mask(rds, [region_geometry], crop=True, nodata=255)
        esri_raw_data = esri_raw_data[0].astype(np.uint8)
        esri_data = np.full(esri_raw_data.shape, 100, np.uint8)
        for c, v in ESRI_CLASSES.items():
            esri_data[esri_raw_data == v] = CLASS_VALUES[c]
        del esri_raw_data
        gc.collect()
        heights.append(esri_data.shape[0])
        widths.append(esri_data.shape[1])


    if 2022 > year > 2019:
        with rio.open(esa_path, 'r') as rds:
            esa_raw_data, _ = mask(rds, [region_geometry], crop=True, nodata=255)
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
    if year > 2014:
        dw_data = dw_data[:min_height, :min_width]
    if year > 2016:
        esri_data = esri_data[:min_height, :min_width]
    if 2022 > year > 2019:
        esa_data = esa_data[:min_height, :min_width]

    for intersect_class in intersect_classes:
        print(intersect_class)
        if 2022 > year > 2019:
            intersect_all = np.where((esa_data == CLASS_VALUES[intersect_class]) & (dw_data == CLASS_VALUES[intersect_class]) &
                                     (esri_data == CLASS_VALUES[intersect_class]),
                                     dw_data, np.full(dw_data.shape, 100, np.uint8))
            with rio.open(os.path.join(save_path, f'{intersect_class}_intersect_all.tif'), 'w', **raster_profile) as wds:
                wds.write(intersect_all, 1)
            del intersect_all
            gc.collect()
            print('intersect all')

            intersect_esa_dw = np.where(
                (esa_data == CLASS_VALUES[intersect_class]) & (dw_data == CLASS_VALUES[intersect_class]),
                esa_data, np.full(dw_data.shape, 100, np.uint8))
            with rio.open(os.path.join(save_path, f'{intersect_class}_intersect_esa_dw.tif'), 'w', **raster_profile) as wds:
                wds.write(intersect_esa_dw, 1)
            del intersect_esa_dw
            gc.collect()
            print('intersect esa dw')

            intersect_esa_esri = np.where(
                (esa_data == CLASS_VALUES[intersect_class]) & (esri_data == CLASS_VALUES[intersect_class]),
                esa_data, np.full(dw_data.shape, 100, np.uint8))
            with rio.open(os.path.join(save_path, f'{intersect_class}_intersect_esa_esri.tif'), 'w', **raster_profile) as wds:
                wds.write(intersect_esa_esri, 1)
            del intersect_esa_esri
            gc.collect()
            print('intersect esa esri')

        if 2023 > year > 2016:
            intersect_dw_esri = np.where(
                (dw_data == CLASS_VALUES[intersect_class]) & (esri_data == CLASS_VALUES[intersect_class]),
                dw_data, np.full(dw_data.shape, 100, np.uint8))
            with rio.open(os.path.join(save_path, f'{intersect_class}_intersect_dw_esri.tif'), 'w', **raster_profile) as wds:
                wds.write(intersect_dw_esri, 1)
            del intersect_dw_esri
            gc.collect()
            print('intersect dw esri')
    if 'esa_data' in globals():
        del esa_data
    if 'dw_data' in globals():
        del dw_data
    if 'esri_data' in globals():
        del esri_data
    gc.collect()
