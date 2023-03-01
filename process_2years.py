import os, shutil
import gc
import numpy as np
from pyproj import CRS
import geopandas as gpd
from shapely.validation import make_valid
import rasterio as rio
from rasterio.mask import mask
from rasterio.merge import merge
from rasterio.features import rasterize
from rasterio.warp import calculate_default_transform, reproject


NODATA = 255
PLACEHOLDER = 100
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
OSM_PATH = '/home/user/osm_rus_2022/RU-60797a9e-20221123-ru-shape/data'
OSM_LAYERS = {
    'polygon': [
        'settlement', 'airport', 'railway-platform', 'parking',
        'water', 'landuse', 'surface', 'vegetation',
        'poi', 'building'
    ],
    'line': [
        'railway', 'water', 'highway'
    ]
}
OSM_CLASS_ORDER = [
    'wetland', 'water', 'vineyard',  # 'trees', 'shrub',
    'bareland', 'built'   # , 'crops'
]
# REGIONS_PATH = os.path.join(OSM_PATH, 'boundary-polygon.shp')
# REGIONS_NAME_PROPERTY = 'NAME_RU'
# REGIONS_PATH = 'work_rus_regions_2023_with_buffer_1km.gpkg'
# REGIONS_NAME_PROPERTY = 'name'
REGIONS_PATH = 'work_se_regions_2023.gpkg'
REGIONS_NAME_PROPERTY = 'Name'


def merge_reproject_crop(load_path, save_path, geometry=None,
                         need_crs='EPSG:3857', need_res=10, delete_raw=False):
    raster_paths = [os.path.join(load_path, f) for f in os.listdir(load_path)
                    if f.endswith('.tif') or f.endswith('.jp2')]
    for raster_path in raster_paths:
        with rio.open(raster_path, 'r') as rds:
            if rds.crs == CRS(need_crs) and rds.res[0] == need_res and rds.res[1] == need_res:
                need_reproject = False
            else:
                need_reproject = True
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
                    'crs': new_profile.crs,
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
    merge(raster_paths, nodata=new_profile['nodata'], dst_path=save_path)
    if geometry is not None:
        folder, file = os.path.split(save_path)
        raw_save_path = os.path.join(folder, f'{os.path.splitext(file)[0]}_raw.tif')
        os.rename(save_path, raw_save_path)
        with rio.open(raw_save_path, 'r') as rds:
            new_profile = rds.profile.copy()
            cropped_data, new_transform = mask(rds, geometry, crop=True)
            new_profile.update({
                'transform': new_transform,
                'width': cropped_data.shape[2],
                'height': cropped_data.shape[1]
            })
        with rio.open(save_path, 'w', **new_profile) as wds:
            wds.write(cropped_data)
    if delete_raw:
        for raster_path in raster_paths:
            os.remove(raster_path)
        if len(os.listdir(load_path)) == 0:
            os.rmdir(load_path)


def intersect_esa_dw_esri_2years(work_path, region_name, years,
                                 regions_path=REGIONS_PATH, regions_name_property=REGIONS_NAME_PROPERTY,
                                 classes=CLASSES, class_values=CLASS_VALUES,
                                 nodata=NODATA, placeholder=PLACEHOLDER):
    if not os.path.exists(work_path):
        os.makedirs(work_path)
    regions_gdf = gpd.read_file(regions_path).to_crs('EPSG:3857')
    region_geometry = regions_gdf.loc[regions_gdf[regions_name_property] == region_name].iloc[0].geometry
    dw = {}
    esa = {}
    esri = {}
    res = {}
    for year in years:
        print(year)
        save_path = os.path.join(work_path, f'{result_name}_{year}')
        if not os.path.exists(save_path):
            os.mkdir(save_path)
        intersect_save_path = os.path.join(save_path, 'intersect_esa_dw_esri.tif')
        dw_path = os.path.join(work_path, 'dw', f'dw_{year}.tif')
        esa_path = os.path.join(work_path, 'esa', f'esa_{year}.tif')
        esri_path = os.path.join(work_path, 'esri', f'esri_{year}.tif')
        # dw_path = os.path.join(work_path, 'dw', str(year))
        # if not os.path.exists(dw_path):
        #     dw_path = os.path.join(work_path, 'dw', f'dw_{year}.tif')
        # esa_path = os.path.join(work_path, 'esa', str(year))
        # if not os.path.exists(esa_path):
        #     esa_path = os.path.join(work_path, 'esa', f'esa_{year}.tif')
        # esri_path = os.path.join(work_path, 'esri', str(year))
        # if not os.path.exists(esri_path):
        #     esri_path = os.path.join(work_path, 'esri', f'esri_{year}.tif')
        # if os.path.isdir(dw_path):
        #     dw_folder = dw_path
        #     dw_path = os.path.join(os.path.split(dw_folder)[0], f'dw_{os.path.split(dw_folder)[1]}.tif')
        #     merge([os.path.join(dw_folder, f) for f in os.listdir(dw_folder)], dst_path=dw_path)
        #     shutil.rmtree(dw_folder)
        #     print('merge DW')
        # if os.path.isdir(esa_path):
        #     esa_folder = esa_path
        #     esa_path = os.path.join(os.path.split(esa_folder)[0], f'esa_{os.path.split(esa_folder)[1]}.tif')
        #     merge([os.path.join(esa_folder, f) for f in os.listdir(esa_folder)], dst_path=esa_path)
        #     shutil.rmtree(esa_folder)
        #     print('merge ESA')
        # if os.path.isdir(esri_path):
        #     esri_folder = esri_path
        #     esri_path = os.path.join(os.path.split(esri_folder)[0], f'esri_{os.path.split(esri_folder)[1]}.tif')
        #     merge([os.path.join(esri_folder, f) for f in os.listdir(esri_folder)], dst_path=esri_path)
        #     shutil.rmtree(esri_folder)
        #     print('merge ESRI')
        with rio.open(dw_path, 'r') as rds:
            dw_data, new_transform = mask(rds, [region_geometry], nodata=nodata, crop=True)
            raster_profile = rds.profile.copy()
        dw_class_data = np.full(dw_data.shape[1:], nodata, np.uint8)
        for i, c in enumerate(classes):
            dw_class_data[dw_data[0] == DW_CLASSES[c]] = i
        del dw_data
        print('dw loaded')
        with rio.open(esa_path, 'r') as rds:
            esa_data, _ = mask(rds, [region_geometry], nodata=nodata, crop=True)
        esa_class_data = np.full(esa_data.shape[1:], nodata, np.uint8)
        for i, c in enumerate(classes):
            esa_class_data[esa_data[0] == ESA_CLASSES[c]] = i
        del esa_data
        print('esa loaded')
        with rio.open(esri_path, 'r') as rds:
            esri_data, _ = mask(rds, [region_geometry], nodata=nodata, crop=True)
        esri_class_data = np.full(esri_data.shape[1:], nodata, np.uint8)
        for i, c in enumerate(classes):
            esri_class_data[esri_data[0] == ESRI_CLASSES[c]] = i
        del esri_data
        print('esri loaded')
        gc.collect()
        min_height = min(dw_class_data.shape[0], esa_class_data.shape[0], esri_class_data.shape[0])
        min_width = min(dw_class_data.shape[1], esa_class_data.shape[1], esri_class_data.shape[1])
        dw[year] = dw_class_data[:min_height, :min_width].copy()
        esa[year] = esa_class_data[:min_height, :min_width].copy()
        esri[year] = esri_class_data[:min_height, :min_width].copy()
        del dw_class_data, esa_class_data, esri_class_data
        gc.collect()
        raster_profile.update({
            'transform': new_transform,
            'height': min_height,
            'width': min_width,
            'compress': 'lzw',
            'count': 1,
            'nodata': nodata
        })
        stages = ['intersect_esa_dw', 'esa_trees', 'dw_esri_trees', 'esa_dw_wetland',
                  'esa_dw_woody_wetland', 'esa_dw_esri_crops', 'all_built']
        for stage_num, stage_name in enumerate(stages):
            if stage_name == 'intersect_esa_dw':
                res_data = np.where(dw[year] == esa[year],
                                    dw[year], np.full(dw[year].shape, placeholder, np.uint8))
            elif stage_name == 'esa_trees':
                res_data = np.where(
                    esa[year] == class_values['trees'],
                    esa[year], res_data
                )
            elif stage_name == 'dw_esri_trees':
                res_data = np.where(
                    (dw[year] == class_values['trees']) & (esri[year] == class_values['trees']) &
                    (res_data == placeholder),
                    esri[year], res_data
                )
            elif stage_name == 'esa_dw_wetland':
                res_data = np.where(
                    ((esa[year] == class_values['water']) | (dw[year] == class_values['water'])) &
                    ((esa[year] == class_values['wetland']) | (dw[year] == class_values['wetland'])) &
                    (res_data == placeholder),
                    np.full(res_data.shape, class_values['wetland'], np.uint8), res_data
                )
            elif stage_name == 'esa_dw_woody_wetland':
                res_data = np.where(
                    ((esa[year] == class_values['water']) | (esa[year] == class_values['wetland'])) &
                    (dw[year] == class_values['trees']),
                    np.full(res_data.shape, class_values['wetland'], np.uint8), res_data
                )
            elif stage_name == 'esa_dw_esri_crops':
                res_data = np.where(
                    (esri[year] == class_values['crops']) &
                    ((esa[year] == class_values['crops']) | (dw[year] == class_values['crops'])) &
                    (res_data == placeholder),
                    esri[year], res_data
                )
            elif stage_name == 'all_built':
                res_data = np.where(
                    (dw[year] == class_values['built']) | (esa[year] == class_values['built']) |
                    (esri[year] == class_values['built']),
                    np.full(dw[year].shape, class_values['built'], np.uint8), res_data
                )
            with rio.open(os.path.join(save_path, f'{stage_num + 1}_{stage_name}.tif'), 'w', **raster_profile) as wds:
                wds.write(res_data, 1)
            print('-->', stage_num + 1, stage_name)
        with rio.open(intersect_save_path, 'w', **raster_profile) as wds:
            wds.write(res_data, 1)
        res[year] = res_data.copy()
        del res_data
        del dw[year], esa[year], esri[year]
        gc.collect()
    print('intersect')
    save_path = os.path.join(work_path, f'{result_name}_{"_".join(years)}')
    if not os.path.exists(save_path):
        os.mkdir(save_path)
    result_save_path = os.path.join(save_path, 'intersect_esa_dw_esri_2years.tif')
    min_width = min(res[years[0]].shape[1], res[years[1]].shape[1])
    min_height = min(res[years[0]].shape[0], res[years[1]].shape[0])
    res[years[0]] = res[years[0]][:min_height, :min_width]
    res[years[1]] = res[years[1]][:min_height, :min_width]
    res[years[1]] = np.where(res[years[0]] == class_values['crops'],
                             res[years[0]], res[years[1]])
    del res[years[0]]
    for year in years:
        if year in dw.keys():
            del dw[year]
        if year in esa.keys():
            del esa[year]
        if year in esri.keys():
            del esri[year]
    gc.collect()
    with rio.open(result_save_path, 'w', **raster_profile) as wds:
        wds.write(res[years[1]], 1)
    print('ready')


def osm_conditions(layer_gdf, layer_name, layer_type='polygon'):
    class_layers = {}
    if layer_name in ('settlement', 'building', 'highway', 'railway', 'airport', 'railway-platform', 'parking'):
        class_layers['built'] = layer_gdf
    elif layer_name == 'water':
        if layer_type == 'polygon':
            class_layers['water'] = layer_gdf.loc[(layer_gdf['NATURAL'] != 'wetland') &
                                                  (layer_gdf['WETLAND'].isnull())
                                                  ]
            class_layers['wetland'] = layer_gdf.loc[(layer_gdf['NATURAL'] == 'wetland') &
                                                    (layer_gdf['WETLAND'].notnull())
                                                    ]
        elif layer_type == 'line':
            class_layers['water'] = layer_gdf.loc[layer_gdf['WATERWAY'] != 'drain']
    elif layer_name == 'landuse':
        class_layers['built'] = layer_gdf.loc[layer_gdf['LANDUSE'].isin(
            ('traffic_island', 'residential', 'education', 'educational', 'garage', 'garages', 'railway',
             'autodrome', 'cemetery', 'apiary', 'border_control', 'retail', 'industrial', 'school',
             'churchyard', 'institutional', 'commercial', 'religious', 'construction', 'depot', 'public', 'highway',
             'parking', 'common', 'storage', 'recreation_ground', 'shed', 'animal_keeping', 'police',
             'healthcare', 'harbour', 'hamlet', 'allotments', 'greenhouse_horticulture'))  # 'shooting_range',
        ]
        class_layers['crops'] = layer_gdf.loc[layer_gdf['LANDUSE'].isin(
            ('farmland', 'farmyard', 'field', 'plantation'))
        ]
        class_layers['trees'] = layer_gdf.loc[layer_gdf['LANDUSE'].isin(
            ('wood', 'orchard', 'plant_nursery', 'national_reserve', 'conservation'))
        ]
        class_layers['grass'] = layer_gdf.loc[layer_gdf['LANDUSE'].isin(
            (
            'logging', 'grassland', 'grass', 'greenfield', 'meadow', 'island', 'village_green', 'flowerbed', 'pasture'))
        ]
        class_layers['bareland'] = layer_gdf.loc[layer_gdf['LANDUSE'].isin(
            ('quarry', 'desert', 'peat_cutting', 'brownfield', 'landfill'))
        ]
        class_layers['shrub'] = layer_gdf.loc[layer_gdf['LANDUSE'].isin(
            ('scrub', 'mud'))
        ]
        class_layers['vineyard'] = layer_gdf.loc[layer_gdf['LANDUSE'] == 'vineyard'
        ]
    elif layer_name == 'surface':
        class_layers['grass'] = layer_gdf.loc[layer_gdf['NATURAL'].isin(
            ('fell', 'grassland'))
        ]
        class_layers['shrub'] = layer_gdf.loc[layer_gdf['NATURAL'] == 'scrub']
        class_layers['bareland'] = layer_gdf[layer_gdf['NATURAL'].isin(
            ('beach', 'sand', 'scree', 'resort'))
        ]
    elif layer_name == 'vegetation':
        class_layers['trees'] = layer_gdf.loc[(layer_gdf['NATURAL'].isin(
            ('forest', 'heath', 'tree', 'tree_group', 'tree_row', 'wood')
        )) | (layer_gdf['LANDUSE'].isin(
            ('conservation', 'forest', 'tree', 'plant_nursery', 'orchard')
        ))
        ]
        class_layers['grass'] = layer_gdf.loc[(layer_gdf['NATURAL'].isin(
            ('grassland', 'grass', 'greenfield', 'meadow', 'village_green', 'quarry', 'logging')
        )) & (layer_gdf['LANDUSE'] != 'farmland')
        ]
        class_layers['shrub'] = layer_gdf.loc[(layer_gdf['NATURAL'].isin(
            ('bush', 'cliff', 'coastline', 'fell', 'scrub')
        ))  # | (layer_gdf['LANDUSE'] == 'vineyard')
        ]
        class_layers['vineyard'] = layer_gdf.loc[layer_gdf['LANDUSE'] == 'vineyard'
        ]
    elif layer_name == 'poi':
        if layer_type == 'polygon':
            class_layers['built'] = layer_gdf.loc[(layer_gdf['LEISURE'].isnull()) |
                                                  ((layer_gdf['MAN_MADE'].notnull()) | (
                                                      layer_gdf['AMENITY'].notnull()) |
                                                   (layer_gdf['OFFICE'].notnull()) | (layer_gdf['SHOP'].notnull()) |
                                                   (layer_gdf['SPORT'].notnull()))
                                                  ]
    return class_layers


def osm_to_vector(work_path, region_name,
                  regions_path=REGIONS_PATH, regions_name_property=REGIONS_NAME_PROPERTY,
                  osm_path=OSM_PATH, osm_layers=OSM_LAYERS, conditions_function=osm_conditions,
                  osm_name='osm'):
    print('osm conditions')
    save_path = os.path.join(work_path, osm_name)
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    regions_gdf = gpd.read_file(regions_path).to_crs('EPSG:4326')
    region_gdf = regions_gdf.loc[regions_gdf[regions_name_property] == region_name]
    region_gdf.to_file(os.path.join(save_path, 'boundary-polygon.gpkg'))
    region_geometry = region_gdf.iloc[0].geometry
    for layer_type, layer_names in osm_layers.items():
        for layer_name in layer_names:
            layer_gdf = gpd.read_file(os.path.join(osm_path, f'{layer_name}-{layer_type}.shp'),
                                      bbox=region_geometry.bounds)  # mask = region_geometry
            print(layer_name, end=' ')
            layer_gdf.geometry = layer_gdf.geometry.map(lambda x: x if x.is_valid else make_valid(x))
            layer_gdf = layer_gdf.loc[layer_gdf.geometry.intersects(region_geometry)]
            print(len(layer_gdf))
            if len(layer_gdf) == 0:
                continue
            if layer_type == 'line':
                layer_gdf.to_crs('EPSG:3857', inplace=True)
                layer_gdf.geometry = layer_gdf.geometry.buffer(5)
                layer_gdf.to_crs('EPSG:4326', inplace=True)
            class_layers = conditions_function(layer_gdf, layer_name, layer_type)
            for class_name, class_layer in class_layers.items():
                class_path = os.path.join(save_path, class_name)
                if not os.path.exists(class_path):
                    os.mkdir(class_path)
                class_layer_path = os.path.join(class_path, f'{layer_name}-{layer_type}.gpkg')
                class_layer.to_file(class_layer_path)
                print('-->', class_name, len(class_layer))
    print('ready')


def osm_to_raster(work_path, region_name, use_settlement=False,
                  regions_path=REGIONS_PATH, regions_name_property=REGIONS_NAME_PROPERTY,
                  nodata=NODATA, placeholder=PLACEHOLDER, classes=OSM_CLASSES, class_order=OSM_CLASS_ORDER,
                  osm_name='osm'):
    print('osm rasterize')
    regions_gdf = gpd.read_file(regions_path).to_crs('EPSG:3857')
    region_geometry = regions_gdf.loc[
        regions_gdf[regions_name_property] == region_name].iloc[0].geometry
    template_mask_path = [os.path.join(work_path, 'dw', f'dw_{y}.tif') for y in range(2017, 2023)
                          if os.path.exists(os.path.join(work_path, 'dw', f'dw_{y}.tif'))][-1]
    osm_vector_path = os.path.join(work_path, osm_name)
    if use_settlement:
        osm_raster_path = os.path.join(work_path, osm_name, 'osm_mask.tif')
    else:
        osm_raster_path = os.path.join(work_path, osm_name, 'osm_mask_no_settlement.tif')
    with rio.open(template_mask_path, 'r') as rds:
        mask_data, raster_transform = mask(rds, [region_geometry], nodata=nodata, crop=True)
        raster_profile = rds.profile.copy()
    osm_mask = np.where(mask_data[0] == nodata,
                        np.full(mask_data.shape[1:], nodata, np.uint8),
                        np.full(mask_data.shape[1:], placeholder, np.uint8))
    del mask_data
    for c in class_order:
        print(c)
        class_folder = os.path.join(osm_vector_path, c)
        if not os.path.exists(class_folder):
            continue
        class_index = classes[c]
        class_paths = [os.path.join(class_folder, f) for f in os.listdir(class_folder)]
        for class_path in class_paths:
            print(class_path)
            cur_gdf = gpd.read_file(class_path).to_crs('EPSG:3857')
            if c == 'built':
                if 'settlement' in class_path and not use_settlement:
                    continue
            cur_features = cur_gdf.geometry.tolist()
            if len(cur_features) > 0:
                cur_mask = rasterize(cur_features, out_shape=osm_mask.shape, fill=placeholder,
                                     transform=raster_transform, default_value=class_index,
                                     dtype=np.uint8)
                osm_mask = np.where(osm_mask == placeholder, cur_mask, osm_mask)
                del cur_mask

    with rio.open(osm_raster_path, 'w', **raster_profile) as wds:
        wds.write(osm_mask, 1)
    print('ready')


def intersect_masks(work_path, years, nodata=PLACEHOLDER, use_settlement=False, osm_name='osm'):
    print('intersect mask and osm')
    if use_settlement:
        osm_path = os.path.join(work_path, osm_name, 'osm_mask.tif')
    else:
        osm_path = os.path.join(work_path, osm_name, 'osm_mask_no_settlement.tif')
    mask_path = os.path.join(work_path, f'{result_name}_{"_".join(years)}', 'intersect_esa_dw_esri_2years.tif')
    save_path = os.path.join(work_path, f'{result_name}_{"_".join(years)}', 'intersect_osm.tif')
    merge([osm_path, mask_path], nodata=nodata, dst_path=save_path)
    print('ready')


if __name__ == '__main__':
    region_name = 'Краснодарский край'  # 'Адыгея' 'Краснодарский край'
    years = ('2020', '2021')
    use_settlement = False
    work_path = os.path.join('/home/user/usage_project/data', region_name)
    result_name = 'result'
    osm_name = 'osm_new'
    # intersect_esa_dw_esri_2years(work_path, region_name, years)
    # osm_to_vector(work_path, region_name, osm_name=osm_name)
    osm_to_raster(work_path, region_name, use_settlement=use_settlement, osm_name=osm_name)
    # intersect_masks(work_path, years, use_settlement=use_settlement, osm_name=osm_name)
