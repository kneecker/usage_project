import geopandas as gpd
import numpy as np
from shapely.geometry import mapping, box


def make_grid(geometry, deltas=(25000, 25000), only_intersecting=True,
              buffer=None, input_crs='EPSG:3857', out_crs='EPSG:3857', cells_crs=None,
              crop_by_geometry=True, return_list=False):
    geometry_series = gpd.GeoSeries([geometry], crs=input_crs)
    if cells_crs is None:
        cells_crs = geometry_series.estimate_utm_crs()
    geometry_series = geometry_series.to_crs(cells_crs)
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
    if return_list:
        grids = grid_gdf.geometry.tolist()
        return grids
    else:
        return grid_gdf


if __name__ == '__main__':
    REGIONS_PATH = 'work_rus_regions_2023_with_buffer_1km.gpkg'
    REGIONS_NAME_PROPERTY = 'name'
    region_name = 'Краснодарский край'
    need_crs = 'EPSG:3857'
    grids_file = f'data/{region_name}/grid_{need_crs}.gpkg'
    regions_gdf = gpd.read_file(REGIONS_PATH).to_crs(need_crs)
    region_geometry = regions_gdf.loc[regions_gdf[REGIONS_NAME_PROPERTY] == region_name].iloc[0].geometry
    grids = make_grid(region_geometry, cells_crs=need_crs)
    grids.to_file(grids_file)
