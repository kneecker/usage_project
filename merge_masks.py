import os, gc
import rasterio as rio
import numpy as np
from numpy import ma
import geopandas as gpd
from rasterio.mask import mask as rio_mask
from rasterio.features import shapes as rio_shapes
from rasterio.features import rasterize
from shapely.geometry import shape as shapely_shape

nodata = 255
CLASSES = [
    'water', 'trees', 'grass', 'wetland', 'crops',
    'shrub', 'built', 'bareland', 'snowice']

REGIONS_PATH = 'work_rus_regions_2023_with_buffer_1km.gpkg'
REGIONS_NAME_PROPERTY = 'name'
region_name = 'Адыгея'
regions_gdf = gpd.read_file(REGIONS_PATH).to_crs('EPSG:3857')
region_geometry = regions_gdf.loc[regions_gdf[REGIONS_NAME_PROPERTY] == region_name].iloc[0].geometry

work_path = os.path.join('/home/user/usage_project/data', region_name)
# year_masks = {}
# for year in (2020, 2021):
#     year_path = os.path.join(work_path, f'test_masks/{year}')
#     paths = [os.path.join(year_path, f'{f}.tif')
#              for f in ['crops_intersect_all', 'crops_intersect_esa_dw', 'crops_intersect_esa_esri']]
#     for p, path in enumerate(paths):
#         print(path)
#         with rio.open(path, 'r') as rds:
#             p_mask, p_transform = rio_mask(rds, region_geometry, nodata=nodata, crop=True)
#             p_mask = p_mask[0].astype(np.uint8)
#             if p == 0:
#                 raster_profile = rds.profile.copy()
#                 raster_transform = p_transform
#                 raster_profile.update(transform=raster_transform, width=p_mask.shape[1], height=p_mask.shape[0])
#         p_mask = np.where(np.isin(p_mask, [CLASSES.index('crops'), nodata]),
#                           p_mask, np.full(p_mask.shape, 100, np.uint8))
#         if p == 0:
#             year_masks[year] = p_mask.copy()
#         else:
#             year_masks[year] = np.where(p_mask == CLASSES.index('crops'), p_mask, year_masks[year])
#     save_path = os.path.join(year_path, f'merge_masks_with_small.tif')
#     with rio.open(save_path, 'w', **raster_profile) as wds:
#         wds.write(year_masks[year], 1)
#     max_areas = [0.5, 1, 2, 3]
#     small_polys = {a: [] for a in max_areas}
#     for shp, value in rio_shapes(year_masks[year], ~np.isin(year_masks[year], [100, nodata]),
#                                  connectivity=4, transform=raster_transform):
#         poly = shapely_shape(shp)
#         for max_area in max_areas:
#             if poly.area * 0.0001 < max_area:
#                 small_polys[max_area].append(poly)
#     for max_area in max_areas:
#         print(max_area)
#         small_mask = rasterize(small_polys[max_area], out_shape=year_masks[year].shape, fill=nodata,
#                                default_value=200, transform=raster_transform, dtype=np.uint8)
#         raster_arr = np.where(small_mask == 200, small_mask, year_masks[year])
#         masked_arr = ma.masked_array(raster_arr, raster_arr == 200)
#         iter_num = 0
#         prev_count = masked_arr.mask.sum()
#         while np.any(masked_arr.mask):
#             for shift in (-1, 1):
#                 for axis in (0, 1):
#                     shifted_arr = np.roll(masked_arr, shift=shift, axis=axis)
#                     shifted_arr[shifted_arr == nodata] = ma.masked
#                     idx = ~shifted_arr.mask * masked_arr.mask
#                     masked_arr[idx] = shifted_arr[idx]
#             cur_count = masked_arr.mask.sum()
#             if prev_count <= cur_count:
#                 break
#             else:
#                 prev_count = cur_count
#             print('iter', iter_num, cur_count)
#             iter_num += 1
#         raster_arr = ma.getdata(masked_arr)
#         raster_arr[raster_arr == 200] = 100
#         save_path = os.path.join(year_path, f'merge_masks_no_small_{max_area}ha.tif')
#         with rio.open(save_path, 'w', **raster_profile) as wds:
#             wds.write(raster_arr, 1)


mask_2020_path = '/home/user/usage_project/data/Краснодарский край/test_masks/2020/merged_masks/merge_masks_no_small_2ha.tif'
mask_2021_path = '/home/user/usage_project/data/Краснодарский край/test_masks/2021/merged_masks/merge_masks_no_small_2ha.tif'
merged_save_path = os.path.join(work_path, 'test_masks', 'merge_2020_2021_no_small_2ha.tif')
with rio.open(mask_2020_path, 'r') as rds:
    mask_2020, raster_transform = rio_mask(rds, region_geometry, nodata=nodata, crop=True)
    mask_2020 = mask_2020[0].astype(np.uint8)
    raster_profile = rds.profile.copy()
    raster_profile.update(transform=raster_transform, width=mask_2020.shape[1], height=mask_2020.shape[0])
with rio.open(mask_2021_path, 'r') as rds:
    mask_2021, _ = rio_mask(rds, region_geometry, nodata=nodata, crop=True)
    mask_2021 = mask_2021[0].astype(np.uint8)
merged_mask = np.where((mask_2020 == CLASSES.index('crops')) | (mask_2021 == CLASSES.index('crops')),
                       np.full(mask_2020.shape, CLASSES.index('crops'), np.uint8),
                       np.full(mask_2020.shape, 100, np.uint8))
with rio.open(merged_save_path, 'w', **raster_profile) as wds:
    wds.write(merged_mask, 1)
