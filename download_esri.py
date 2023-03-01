import requests
import os
import geopandas as gpd
from tqdm import tqdm
from utils import run_function_in_multiprocess


def download_esri(utm_cell_paths):
    downloaded_cells = {}
    for utm_cell_path in utm_cell_paths:
        utm_cell_file = os.path.split(utm_cell_path)[1]
        utm_cell, year = utm_cell_file.split('_')
        year = year[:4]
        if os.path.exists(utm_cell_path):
            print(utm_cell_file, 'уже загружен')
        else:
            download_url = f'https://lulctimeseries.blob.core.windows.net/lulctimeseriespublic/' \
                           f'lc{year}/{utm_cell_file}'
            r = requests.get(download_url, headers={'User-agent': 'Mozilla/5.0'}, stream=True)  # allow_redirects=True,
            total_length = int(r.headers.get('content-length', 0))
            with open(utm_cell_path, 'wb') as wbf, tqdm(
                    desc=utm_cell_file,
                    total=total_length,
                    unit='iB',
                    unit_scale=True,
                    unit_divisor=1024) as bar:
                try:
                    for chunk in r.iter_content(chunk_size=1024):
                        if chunk:
                            size = wbf.write(chunk)
                            bar.update(size)
                            wbf.flush()
                except Exception as exc:
                    print(utm_cell_file, 'ошибка', exc)
                    if os.path.exists(utm_cell_path):
                        os.remove(utm_cell_path)
                    continue
                if r.status_code == 200:
                    print(utm_cell_file, 'готово')
                    if year in downloaded_cells.keys():
                        downloaded_cells[year].append(utm_cell)
                    else:
                        downloaded_cells[year] = [utm_cell]
                else:
                    print(utm_cell_file, 'ошибка', r.content)
                    if os.path.exists(utm_cell_path):
                        os.remove(utm_cell_path)
    return downloaded_cells


if __name__ == '__main__':
    regions_path = 'work_se_regions_2023.gpkg'
    regions_name_property = 'Name'
    # regions_path = 'work_rus_regions_2023_with_buffer_1km.gpkg'
    # regions_name_property = 'name'
    regions_gdf = gpd.read_file(regions_path).to_crs('EPSG:4326')
    utm_path = 'vectors/ESRI_cells.gpkg'
    utm_name_property = 'name'
    utm_gdf = gpd.read_file(utm_path).to_crs('EPSG:4326')
    masks_path = '/home/user/usage_project/data/esri'

    # region_name = 'Приморский край'
    region_name = 'Херсонская область'
    region_geometry = regions_gdf.loc[regions_gdf[regions_name_property] == region_name].iloc[0].geometry
    region_utm_cells = utm_gdf.loc[utm_gdf.intersects(region_geometry)][utm_name_property].tolist()
    # years = list(range(2022, 2016, -1))
    years = list(range(2021, 2019, -1))
    num_threads = 2
    utm_cell_paths = []
    for year in years:
        year_cell_paths = [os.path.join(masks_path, str(year), f'{region_utm_cell}_{year}0101-{year + 1}0101.tif')
                           for region_utm_cell in region_utm_cells
                           if not os.path.exists(os.path.join(masks_path, str(year),
                                                              f'{region_utm_cell}_{year}0101-{year + 1}0101.tif'))]
        print(year, ': загружаются ячейки', ', '.join([os.path.split(u)[1].split('_')[0] for u in year_cell_paths]))
        utm_cell_paths += year_cell_paths
    downloaded_cells = run_function_in_multiprocess(download_esri, utm_cell_paths,
                                                    None, num_threads, 'dict')
