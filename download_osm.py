import os
import pandas as pd
import geopandas as gpd
import osmnx as ox
import matplotlib.pyplot as plt
import contextily as cx
import numpy as np


osm_path = 'data/osm'
if not os.path.exists(osm_path):
    os.makedirs(osm_path)
# https://wiki.openstreetmap.org/wiki/Map_features
map_features = {
    'usage': {
        'amenity': True,
        'building': True,
        'aeroway': True,
        'geological': True,
        'highway': True,
        'landuse': [
            'commercial', 'construction', 'education', 'industrial', 'residential', 'retail', 'institutional',
            'aquaculture', 'allotments', 'farmland', 'farmyard', 'flowerbed', 'forest', 'greenhouse_horticulture',
            'meadow', 'orchard', 'plant_nursery', 'vineyard',
            'basin', 'reservoir', 'salt_pond',
            'brownfield', 'cemetery', 'depot', 'grass', 'greenfield', 'landfill', 'port', 'quarry', 'railway',
            'recreation_ground', 'religious', 'village_green'
        ],
        'leisure': ['beach_resort', 'garden', 'park', 'pitch', 'swimming_pool'],
        'man_made': True,
        'military': ['academy', 'airfield', 'barracks', 'office', 'trench'],
        'natural': True,
        'office': True,
        'public_transport': True,
        'railway': True,
        'waterway': True
    },
    'admin': {
        'place': [
            'state', 'region', 'province',
            'district', 'county',
            'municipality',
            'city',
            'town', 'village', 'hamlet',
            'borough', 'suburb',
            'quarter', 'neighbourhood', 'block'
            'isolated_dwelling', 'farm', 'allotments'
        ]
    }
}

need_columns = {
    'usage': [
        'geometry', 'element_type', 'osmid', 'name', 'name:ru', 'name:uk',
        *map_features['usage'].keys()
    ],
    'admin': [
        'geometry', 'element_type', 'osmid', 'name', 'name:ru', 'name:uk',
        'is_in', 'population', 'admin_level'
    ]
}

region_names = ['Шебекино']
# region_names = ['Донецкая область', 'Луганская область', 'Запорожская область', 'Херсонская область']
for region_name in region_names:
    region_gdf = ox.geocode_to_gdf(region_name)
    utm_crs = region_gdf.estimate_utm_crs()
    aoi = region_gdf.loc[0].geometry
    for layer, features in map_features.items():
        layer_path = os.path.join(osm_path, f'osm_{layer}_{region_name}.gpkg')
        layer_gdf = None
        for key, values in features.items():
            print(key, end=' ')
            objects = ox.geometries_from_polygon(aoi, tags={key: values})
            objects = objects.loc[(objects.geometry.type != 'Point') & (objects.geometry.type != 'MultiPoint')]
            objects_lines = objects.loc[(objects.geometry.type == 'LineString') |
                                        (objects.geometry.type == 'MultiLineString')]
            objects = objects.loc[~((objects.geometry.type == 'LineString') |
                                    (objects.geometry.type == 'MultiLineString'))]
            if len(objects_lines) > 0:
                objects_lines.to_crs(utm_crs, inplace=True)
                objects_lines['geometry'] = objects_lines['geometry'].buffer(3, cap_style=2, join_style=3)
                objects_lines.to_crs(objects.crs, inplace=True)
                objects = gpd.GeoDataFrame(
                    pd.concat((objects, objects_lines), ignore_index=False),  crs=objects.crs)
            objects.reset_index(inplace=True)
            if len(objects) > 0:
                if layer_gdf is None:
                    layer_gdf = objects.copy()
                else:
                    layer_gdf = gpd.GeoDataFrame(
                        pd.concat((layer_gdf, objects), ignore_index=False).drop_duplicates(
                            subset=['element_type', 'osmid']), crs=layer_gdf.crs)
            print(len(objects))
        print('готово', len(layer_gdf))
        # delete_columns = np.array(layer_gdf.columns)[~np.isin(np.array(layer_gdf.columns), need_columns[layer])]
        # layer_gdf.drop(columns=delete_columns, inplace=True)
        layer_gdf = layer_gdf[need_columns[layer]]
        layer_gdf.to_file(layer_path)
