import os, io
import requests
import ee
import numpy as np
import zipfile
from shapely.geometry import box, mapping
import geopandas as gpd
from transliterate import translit


NEED_RES = 10
MAX_PIXELS = 5000


def make_grid(geometry, deltas=(MAX_PIXELS * NEED_RES, MAX_PIXELS * NEED_RES), only_intersecting=True,
              buffer=None, input_crs='EPSG:3857', out_crs='EPSG:3857', cells_crs='EPSG:3857',
              crop_by_geometry=True):
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


def download_esa_by_aoi(aoi, year, save_path, aoi_name, use_gdrive=False):
    aoi_name = translit(str(aoi_name), 'ru', reversed=True).replace(' ', '_').replace('\'', '')
    if year in (2020, 2021):
        if year == 2021:
            esa = ee.ImageCollection('ESA/WorldCover/v200').first()
            esa_name = 'esa_2021'
        else:
            esa = ee.ImageCollection('ESA/WorldCover/v100').first()
            esa_name = 'esa_2020'
        esa_crop = esa.clip(aoi)
        download_params = {
            'scale': NEED_RES,
            'crs': 'EPSG:3857',
            'region': aoi
        }
        if use_gdrive:
            download_params.update({
                'image': esa_crop,
                'description': aoi_name + '_' + esa_name,
                'fileNamePrefix': aoi_name + '_' + esa_name,
                'maxPixels': 1e13
            })
            task = ee.batch.Export.image.toDrive(**download_params)
            task.start()
        else:
            download_url = esa_crop.getDownloadURL(download_params)
            r = requests.get(download_url, headers={'User-agent': 'Mozilla/5.0'}, allow_redirects=True)
            if r.status_code == 200:
                with zipfile.ZipFile(io.BytesIO(r.content)) as zf:
                    zf.extract('2020.Map.tif', save_path)
                os.rename(os.path.join(save_path, '2020.Map.tif'),
                          os.path.join(save_path, f'{esa_name}_{aoi_name}.tif'))


def download_dw_by_aoi(aoi, date_from, date_to, save_path, aoi_name, use_gdrive=False):
    aoi_name = translit(str(aoi_name), 'ru', reversed=True).replace(' ', '_').replace('\'', '')
    year = int(date_from[:4])
    if year > 2014:
        dw_name = f'dw_{year}'
        aoi_filter = ee.Filter.And(
            ee.Filter.bounds(aoi),
            ee.Filter.date(date_from, date_to)
        )
        dw = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1').filter(aoi_filter)
        dw_crop = dw.select('label').mode().clip(aoi)
        download_params = {
            'scale': NEED_RES,
            'crs': 'EPSG:3857',
            'region': aoi
        }
        if use_gdrive:
            download_params.update({
                'image': dw_crop,
                'description': aoi_name + '_' + dw_name,
                'fileNamePrefix': aoi_name + '_' + dw_name,
                'maxPixels': 1e13
            })
            task = ee.batch.Export.image.toDrive(**download_params)
            task.start()
        else:
            download_url = dw_crop.getDownloadURL(download_params)
            r = requests.get(download_url, headers={'User-agent': 'Mozilla/5.0'}, allow_redirects=True)
            if r.status_code == 200:
                with zipfile.ZipFile(io.BytesIO(r.content)) as zf:
                    zf.extract('download.label.tif', save_path)
                os.rename(os.path.join(save_path, 'download.label.tif'),
                          os.path.join(save_path, f'{dw_name}_{aoi_name}.tif'))


if __name__ == '__main__':
    regions_path = 'work_se_regions_2023.gpkg'
    regions_name_property = 'Name'
    # regions_path = 'work_rus_regions_2023_with_buffer_1km.gpkg'
    # regions_name_property = 'name'
    regions_gdf = gpd.read_file(regions_path)
    # ee.Authenticate()
    ee.Initialize()

    use_gdrive = True
    download_esa = True
    download_dw = True
    region_name = 'Херсонская область'
    region_path = os.path.join('masks', region_name)
    # years = list(range(2022, 2014, -1))
    years = list(range(2021, 2019, -1))
    for year in years:
        date_from = f'{year}-04-01'
        date_to = f'{year}-11-01'
        if use_gdrive:
            region_geometry = regions_gdf.loc[regions_gdf[regions_name_property] == region_name].iloc[0].geometry
            aoi = ee.Geometry(mapping(region_geometry))
            if download_esa:
                download_esa_by_aoi(aoi, year, region_path, region_name, use_gdrive)
                print(f'download esa {year} to gdrive')
            if download_dw:
                download_dw_by_aoi(aoi, date_from, date_to, region_path, region_name, use_gdrive)
                print(f'download dw {year} to gdrive')
        else:
            if not os.path.exists(region_path):
                os.makedirs(region_path)
            regions_gdf = regions_gdf.to_crs('EPSG:3857')
            region_geometry = regions_gdf.loc[regions_gdf[regions_name_property] == region_name].iloc[0].geometry
            width = (region_geometry.bounds[2] - region_geometry.bounds[0]) / NEED_RES
            height = (region_geometry.bounds[3] - region_geometry.bounds[1]) / NEED_RES
            region_grids = make_grid(region_geometry, out_crs='EPSG:4326')
            aois = [ee.Geometry(mapping(grid)) for i, grid in enumerate(region_grids)]
            for i, aoi in enumerate(aois):
                if download_esa:
                    if i == 0:
                        print(f'download esa {year} to local')
                    download_esa_by_aoi(aoi, year, region_path, i, use_gdrive)
                if download_dw:
                    if i == 0:
                        print(f'download dw {year} to local')
                    download_dw_by_aoi(aoi, date_from, date_to, region_path, i, use_gdrive)
                percent = int(i / len(aois) * 100)
                if percent % 5 == 0:
                    print(f'{percent}%', end='..')
