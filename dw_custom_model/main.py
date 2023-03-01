import os
import tensorflow as tf
import numpy as np
import ee
from shapely.geometry import mapping, box
import geopandas as gpd
import matplotlib.pyplot as plt
import rasterio as rio
from rasterio.warp import calculate_default_transform


MAX_PIXELS = 500
NEED_RES = 10


def make_grid(geometry, deltas=(MAX_PIXELS * NEED_RES, MAX_PIXELS * NEED_RES), only_intersecting=True,
              buffer=None, input_crs='EPSG:4326', out_crs='EPSG:4326', cells_crs='EPSG:3857',
              crop_by_geometry=False):
    if cells_crs is None:
        cells_crs = out_crs
    geometry_series = gpd.GeoSeries([geometry], crs=input_crs).to_crs(cells_crs)
    geometry = geometry_series.iloc[0]
    min_x, min_y, max_x, max_y = geometry.bounds
    if buffer is not None:
        min_x -= buffer
        min_y -= buffer
        max_x += buffer
        max_y += buffer
    grid_polys = []
    for cur_x in np.arange(min_x, max_x, deltas[0]):
        for cur_y in np.arange(min_y, max_y, deltas[1]):
            poly = box(cur_x, cur_y, cur_x + deltas[0], cur_y + deltas[1])
            if only_intersecting:
                if geometry.intersects(poly):
                    if crop_by_geometry:
                        poly = poly.intersection(geometry)
                    grid_polys.append(poly)
            else:
                grid_polys.append(poly)
    grid_features = [
        {'type': 'Feature', 'properties': {'Name': f'G_{g_n}'}, 'geometry': mapping(g_p)}
        for g_n, g_p in enumerate(grid_polys)
    ]
    grid_gdf = gpd.GeoDataFrame.from_features(grid_features, crs=cells_crs).to_crs(out_crs)
    grids = grid_gdf.geometry.tolist()
    return grids


def get_s2_by_aoi(aoi, date_from, date_to, process=True, bands=None):
    if bands is None:
        bands = ('B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B11', 'B12')
    s2 = ee.ImageCollection('COPERNICUS/S2')
    s2_collection = s2.filterBounds(aoi).filterDate(date_from, date_to).limit(3, 'CLOUDY_PIXEL_PERCENTAGE')

    s2_image = s2_collection.select(*bands).first()
    s2_image = s2_image.toFloat().resample('bilinear').reproject(
        'EPSG:3857', scale=10
    )
    s2_image_sample = s2_image.toArray().sampleRectangle(aoi)

    image = np.array(s2_image_sample.getInfo()['properties']['array'])
    print(image.shape)
    if process:
        original_image = image.copy()
        NORM_PERCENTILES = np.array([
            [1.7417268007636313, 2.023298706048351],
            [1.7261204997060209, 2.038905204308012],
            [1.6798346251414997, 2.179592821212937],
            [1.7734969472909623, 2.2890068333026603],
            [2.289154079164943, 2.6171674549378166],
            [2.382939712192371, 2.773418590375327],
            [2.3828939530384052, 2.7578332604178284],
            [2.1952484264967844, 2.789092484314204],
            [1.554812948247501, 2.4140534947492487]
        ])

        image = np.log(image * 0.005 + 1)
        image = (image - NORM_PERCENTILES[:, 0]) / NORM_PERCENTILES[:, 1]
        image = np.exp(image * 5 - 1)
        image = image / (image + 1)
        return original_image, image
    else:
        return image


if __name__ == '__main__':
    regions_name_property = 'name'
    regions_path = '/home/user/usage_project/work_rus_regions_2023_with_buffer_1km.gpkg'
    regions_gdf = gpd.read_file(regions_path)
    ee.Initialize()
    region_name = 'Адыгея'
    date_from = '2020-06-01'
    date_to = '2020-06-30'
    save_path = os.path.join('ee_test_data', region_name)
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    region_geometry = regions_gdf.loc[regions_gdf[regions_name_property] == region_name].iloc[0].geometry
    grids = make_grid(region_geometry)
    for g, grid in enumerate(grids):
        print('grid', g)
        grid_path = os.path.join(save_path, f'{g}.tif')
        if os.path.exists(grid_path):
            print('уже существует')
            continue
        aoi = ee.Geometry(mapping(grid))
        original_image, image = get_s2_by_aoi(aoi, date_from, date_to)
        grid_transform, grid_width, grid_height = calculate_default_transform('EPSG:4326', 'EPSG:3857',
                                                                              image.shape[1], image.shape[0],
                                                                              *grid.bounds, resolution=NEED_RES)
        grid_profile = {
            'driver': 'GTiff',
            'crs': 'EPSG:3857',
            'width': grid_width,
            'height': grid_height,
            'count': image.shape[-1],
            'transform': grid_transform,
            'compress': 'lzw',
            'dtype': 'float32'
        }
        with rio.open(grid_path, 'w', **grid_profile) as wds:
            for b in range(image.shape[-1]):
                wds.write(image[..., b], b + 1)
        f, axarr = plt.subplots(1, 2)
        f.set_size_inches(20, 10)

        axarr[0].imshow(original_image[:, :, [2, 1, 0]] / 3000)
        axarr[0].axes.get_xaxis().set_ticks([])
        axarr[0].axes.get_yaxis().set_ticks([])
        axarr[0].axes.set_xlabel('Original Image')

        axarr[1].imshow(image[:, :, [2, 1, 0]])
        axarr[1].axes.get_xaxis().set_ticks([])
        axarr[1].axes.get_yaxis().set_ticks([])
        axarr[1].axes.set_xlabel('Normalized Image')

        plt.show()
        exit()
