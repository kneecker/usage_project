import os, shutil
from pyproj import CRS
import geopandas as gpd
from shapely.geometry import box
import rasterio as rio
from rasterio.merge import merge
from rasterio.mask import mask
from rasterio.warp import calculate_default_transform, reproject


# OSM_PATH = '/home/user/osm_rus_2022/RU-60797a9e-20221123-ru-shape/data'
# REGIONS_PATH = os.path.join(OSM_PATH, 'boundary-polygon.shp')
# REGIONS_NAME_PROPERTY = 'NAME_RU'
# REGIONS_PATH = 'work_rus_regions_2023_with_buffer_1km.gpkg'
# REGIONS_NAME_PROPERTY = 'name'
REGIONS_PATH = 'work_se_regions_2023.gpkg'
REGIONS_NAME_PROPERTY = 'Name'


def merge_reproject_crop(load_path, save_path, geometry=None,
                         need_crs='EPSG:3857', need_res=10):
    raster_paths = [os.path.join(load_path, f) for f in os.listdir(load_path)
                    if f.endswith('.tif') or f.endswith('.jp2')]
    for raster_path in raster_paths.copy():
        print(raster_path, end=' ')
        with rio.open(raster_path, 'r') as rds:
            raster_extent = gpd.GeoSeries([box(*rds.bounds)], crs=rds.crs).to_crs(need_crs)
            if geometry is not None:
                if not raster_extent.iloc[0].intersects(geometry):
                    del raster_paths[raster_paths.index(raster_path)]
                    print('не пересекает геометрию')
                    continue
            need_reproject = \
                False if rds.crs == CRS(need_crs) and rds.res[0] == need_res and rds.res[1] == need_res else True
        if need_reproject:
            folder, file = os.path.split(raster_path)
            raw_raster_path = os.path.join(folder, f'{os.path.splitext(file)[0]}_raw.tif')
            os.rename(raster_path, raw_raster_path)
            with rio.open(raw_raster_path, 'r') as rds:
                new_transform, new_width, new_height = calculate_default_transform(
                    rds.crs, need_crs, rds.width, rds.height, *rds.bounds, resolution=10
                )
                new_profile = rds.profile.copy()
                new_profile.update({
                    'driver': 'GTiff',
                    'transform': new_transform,
                    'width': new_width,
                    'height': new_height,
                    'crs': need_crs,
                    'compress': 'lzw'
                })
                with rio.open(raster_path, 'w', **new_profile) as wds:
                    for b in range(1, wds.count + 1):
                        reproject(
                            rio.band(rds, b),
                            rio.band(wds, b),
                            src_transform=rds.transform,
                            src_crs=rds.crs,
                            dst_transform=new_transform,
                            dst_crs=need_crs,
                            num_threads=4
                        )
            os.remove(raw_raster_path)
        print('Готово')
    if geometry is not None:
        for raster_path in raster_paths:
            folder, file = os.path.split(raster_path)
            raster_reprojected_path = os.path.join(folder, f'{os.path.splitext(file)[0]}_reprojected.tif')
            os.rename(raster_path, raster_reprojected_path)
            with rio.open(raster_reprojected_path, 'r') as rds:
                new_profile = rds.profile.copy()
                cropped_data, new_transform = mask(rds, [geometry], crop=True)
                new_profile.update({
                    'transform': new_transform,
                    'width': cropped_data.shape[2],
                    'height': cropped_data.shape[1]
                })
            with rio.open(raster_path, 'w', **new_profile) as wds:
                wds.write(cropped_data)
        print('начало объединения')
        if len(raster_paths) > 1:
            merge(raster_paths, nodata=new_profile['nodata'], dst_path=save_path)
        else:
            shutil.copy(raster_paths[0], save_path)
        for raster_path in raster_paths:
            folder, file = os.path.split(raster_path)
            raster_reprojected_path = os.path.join(folder, f'{os.path.splitext(file)[0]}_reprojected.tif')
            os.remove(raster_path)
            os.rename(raster_reprojected_path, raster_path)
        print('готово')


if __name__ == '__main__':
    region_name = 'Херсонская область'
    mask_type = 'esri'
    # years = list(range(2022, 2016, -1))
    # mask_type = 'esa'
    # years = list(range(2021, 2019, -1))  # esa
    # mask_type = 'dw'
    # years = list(range(2022, 2014, -1))  # dw
    years = list(range(2021, 2019, -1))
    for year in years:
        print(year)
        load_path = f'/home/user/usage_project/data/{mask_type}/{year}'
        # load_path = f'/home/user/usage_project/data/{region_name}/{mask_type}/{mask_type}_{year}'
        save_path = f'/home/user/usage_project/data/{region_name}/{mask_type}/{mask_type}_{year}.tif'
        if not os.path.exists(os.path.split(save_path)[0]):
            os.makedirs(os.path.split(save_path)[0])
        regions_gdf = gpd.read_file(REGIONS_PATH).to_crs('EPSG:3857')
        region_geometry = regions_gdf.loc[regions_gdf[REGIONS_NAME_PROPERTY] == region_name].iloc[0].geometry
        merge_reproject_crop(load_path, save_path, region_geometry)
