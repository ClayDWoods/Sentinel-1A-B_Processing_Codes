import sys, os

sys.path.append("..")
os.environ['PROJ_LIB'] = '/projects/clwo4142/.conda_pkgss/proj4-4.9.3-hc8507d1_7/share/proj'

import sys, csv
import os, os.path
import rasterio
import xml.dom.minidom
import pyspark
import tempfile, shutil
import time
import fiona
import pandas as pd
import glob
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from rasterio.plot import show
import seaborn as sns

plt.style.use('seaborn-whitegrid')
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys, csv
import os, os.path
import rasterio
import pyspark
import tempfile, shutil
import time
import geopandas as gpd
import pandas as pd
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from rasterio.plot import show
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors
from rasterio.warp import calculate_default_transform, reproject, Resampling
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from cartopy import config
import cartopy.crs as ccrs
import math
from matplotlib.patches import PathPatch
from matplotlib.path import Path as Path_m
import sys, csv
import os, os.path
from os import mkdir
import rasterio
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pyspark
import subprocess
import urllib.request
import re
from descartes import PolygonPatch
import shapefile as shp
import tempfile, shutil
import base64
import time
import getpass
import ssl
import signal
import xml.etree.ElementTree as ET
import pandas as pd
import glob
import zipfile
import scipy
from tqdm import tqdm_notebook as tqdm
import math
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio import Affine as Affine
import itertools
from multiprocessing import Pool
# from rasterio.windows import Window
import statsmodels.api as sm
import statsmodels.formula.api as smf
from sys import getsizeof
from time import sleep
from IPython.display import clear_output
import datetime
from time import mktime
import time
import multiprocessing as mp
import timeit
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon as Polygon_mpl
from shapely import geometry
import shapely
from matplotlib.collections import LineCollection

tqdm().pandas()

import shutil
from pathlib import Path
from general_set_up_functions import rasterio_basic_functions, general_functions
import math
from statistics import mean
import numpy as np
from datetime import datetime as dt


def read_geologic_map_shp(geologic_map_path):
    map_df = gpd.read_file(geologic_map_path).sort_values(by=['GLG']).reset_index(drop=True)
    not_use_units = ['H2O', 'lake', 'U', 'QTdd', 'QTdi', 'QTdt', 'QTg', 'Qsm', 'Tb', 'Tba', 'Tbb', 'Tj', 'Tt', 'dm',
                     'csd', 'ava', 'afo']
    # , 'dsd', 'afy', 'asd
    geol_units = list(map_df.GLG.unique())
    geol_units_new = geol_units.copy()
    for geol_unit in geol_units:
        for n_un in not_use_units:
            if geol_unit == n_un:
                geol_units_new.remove(geol_unit)
    geol_units_ints = list(range(0, len(geol_units_new)))
    dfs = []
    for geol_unit, geol_unit_int in zip(geol_units_new, geol_units_ints):
        partit = map_df[map_df.GLG == str(geol_unit)]
        partit['GLG_int'] = geol_unit_int
        dfs.append(partit)
    frame = pd.concat(dfs).sort_values(by=['GLG'], ascending=True).reset_index(drop=True)
    cmap = plt.get_cmap('terrain', len(geol_units_ints))
    vmin, vmax = 0, 26
    return frame, cmap, vmin, vmax


def normalize(array, percentile):
    lower_percentile, upper_percentile = percentile, 100 - percentile
    array_min, array_max = np.nanpercentile(array, lower_percentile), np.nanpercentile(array, upper_percentile)
    array_norm = (array - array_min) / (array_max - array_min)
    array_norm[array_norm > 1] = 1
    array_norm[array_norm < 0] = 0
    return array_norm


def sinusoid(x, A, offset, omega, phase):
    return A * np.sin(omega * x + phase) + offset


def get_p0(x, y, T):
    A0 = (max(y[0:T]) - min(y[0:T])) / 2
    offset0 = y[0]
    phase0 = 0
    omega0 = 2. * np.pi / T
    return [A0, offset0, omega0, phase0]


def return_frame_coords(list_of_frame_xmls):
    polygons_coordinates = []
    for frame_xml in list_of_frame_xmls:
        tree = ET.parse(frame_xml)
        root = tree.getroot()
        for coordinate in root.iter('coordinates'):
            coordinates_text = coordinate.text
            parse_1_coordinates = coordinates_text.split(' ')
            parse_2_coordinates = [x.split(',') for x in parse_1_coordinates]
            coordinates_correct_order = [parse_2_coordinates[2], parse_2_coordinates[3], parse_2_coordinates[0],
                                         parse_2_coordinates[1]]
            coord_float_pairs = []
            for coordinate_pair in coordinates_correct_order:
                lon = float(coordinate_pair[0])
                lat = float(coordinate_pair[1])
                coord_pair_float = [lon, lat]
                coord_float_pairs.append(coord_pair_float)
            polygon_of_coordinates = np.array(coord_float_pairs)
        polygons_coordinates.append(polygon_of_coordinates)
    return polygons_coordinates


def make_cross_section_line(path_2_defo, deformation_array, topo_array, point_lons, point_lats, width_):
    lon1, lat1 = point_lons[0], point_lats[0]
    lon2, lat2 = point_lons[1], point_lats[1]
    start = -1
    end = 1
    width = end - start
    base = 0.1
    if lon1 == lon2:
        if lat1 > lat2:
            lat_1_temp = lat2
            lat_2_temp = lat1
            lat1 = lat_1_temp
            lat2 = lat_2_temp
        col1, row1 = rasterio_basic_functions.lon_lat_2_col_row(path_2_defo, lon1, lat1)
        col2, row2 = rasterio_basic_functions.lon_lat_2_col_row(path_2_defo, lon2, lat2)
        col2 = col1
        if row1 > row2:
            temp_row_1 = row2
            temp_row_2 = row1
            row1 = temp_row_1
            row2 = temp_row_2
        index_ = range(row1, row2)
        if width_ > 0:
            deformation_array_extract = np.mean(deformation_array[row1:row2, (col2 - width_):(col2 + width_ + 1)],
                                                axis=1)
            topo_array_extract = np.mean(topo_array[row1:row2, (col2 - width_):(col2 + width_ + 1)], axis=1)
        else:
            deformation_array_extract = deformation_array[row1:row2, col2]
            topo_array_extract = topo_array[row1:row2, col2]
        end_coord, start_coord = (general_functions.nearest_bound_coord(lat2, base, 'up'),
                                  general_functions.nearest_bound_coord(lat1, base, 'down'))
        deformation_array_extract[deformation_array_extract == 0] = np.nan
        mean_def = np.nanmean(deformation_array_extract)
    else:
        if lon1 > lon2:
            lon_1_temp = lon2
            lon_2_temp = lon1
            lon1 = lon_1_temp
            lon2 = lon_2_temp
        col1, row1 = rasterio_basic_functions.lon_lat_2_col_row(path_2_defo, lon1, lat1)
        col2, row2 = rasterio_basic_functions.lon_lat_2_col_row(path_2_defo, lon2, lat2)
        row2 = row1
        if col1 > col2:
            temp_col_1 = col2
            temp_col_2 = col1
            col1 = temp_col_1
            col2 = temp_col_2
        index_ = range(col1, col2)
        if width_ > 0:
            deformation_array_extract = np.mean(deformation_array[(row2 - width_):(row2 + width_ + 1), col1:col2],
                                                axis=0)
            topo_array_extract = np.mean(topo_array[(row2 - width_):(row2 + width_ + 1), col1:col2], axis=0)
        else:
            deformation_array_extract = deformation_array[row2, col1:col2]
            topo_array_extract = topo_array[row2, col1:col2]

        end_coord, start_coord = (general_functions.nearest_bound_coord(lon2, base, 'up'),
                                  general_functions.nearest_bound_coord(lon1, base, 'down'))
        deformation_array_extract[deformation_array_extract == 0] = np.nan
        mean_def = np.nanmean(deformation_array_extract)
    return ((deformation_array_extract + 2 * mean_def) * 1000,
            topo_array_extract, index_, round(start_coord, 1), round(end_coord, 1))


# def make_cross_section_line(path_2_defo, deformation_array, topo_array, point_lons, point_lats, width_):
#     lon1, lat1 = point_lons[0], point_lats[0]
#     lon2, lat2 = point_lons[1], point_lats[1]
#     start = -1
#     end = 1
#     width = end - start
#     base = 0.1
#     if lon1 == lon2:
#         if lat1 > lat2:
#             lat_1_temp = lat2
#             lat_2_temp = lat1
#             lat1 = lat_1_temp
#             lat2 = lat_2_temp
#         col1, row1 = rasterio_basic_functions.lon_lat_2_col_row(path_2_defo, lon1, lat1)
#         col2, row2 = rasterio_basic_functions.lon_lat_2_col_row(path_2_defo, lon2, lat2)
#         col2 = col1
#         if row1 > row2:
#             temp_row_1 = row2
#             temp_row_2 = row1
#             row1 = temp_row_1
#             row2 = temp_row_2
#         index_ = range(row1, row2)
#         if width_ > 0:
#             deformation_array_extract = np.mean(deformation_array[row1:row2, (col2-width_):(col2+width_)], axis=1)
#             topo_array_extract = np.mean(topo_array[row1:row2, (col2-width_):(col2+width_)], axis=1)
#         else:
#             deformation_array_extract = deformation_array[row1:row2, col2]
#             topo_array_extract = topo_array[row1:row2, col2]
#         def_min = np.min(deformation_array_extract)
#         def_max = np.max(deformation_array_extract)
#         topo_min = np.min(topo_array_extract)
#         topo_max = np.max(topo_array_extract)
#         deformation_array_extract_norm = (((deformation_array_extract - def_min)
#                                            / (def_max - def_min)) * width + start)
#         topo_array_extract_norm = (((topo_array_extract - topo_min)
#                                     / (topo_max - topo_min)) * width + start)
#         end_coord, start_coord = (general_functions.nearest_bound_coord(lat2, base, 'up'),
#                                   general_functions.nearest_bound_coord(lat1, base, 'down'))
#     else:
#         if lon1 > lon2:
#             lon_1_temp = lon2
#             lon_2_temp = lon1
#             lon1 = lon_1_temp
#             lon2 = lon_2_temp
#         col1, row1 = rasterio_basic_functions.lon_lat_2_col_row(path_2_defo, lon1, lat1)
#         col2, row2 = rasterio_basic_functions.lon_lat_2_col_row(path_2_defo, lon2, lat2)
#         row2 = row1
#         if col1 > col2:
#             temp_col_1 = col2
#             temp_col_2 = col1
#             col1 = temp_col_1
#             col2 = temp_col_2
#         index_ = range(col1, col2)
#         if width_ > 0:
#             deformation_array_extract =  np.mean(deformation_array[(row2-width_):(row2+width_), col1:col2], axis=0)
#             topo_array_extract =  np.mean(topo_array[(row2-width_):(row2+width_), col1:col2], axis=0)
#         else:
#             deformation_array_extract = deformation_array[row2, col1:col2]
#             topo_array_extract = topo_array[row2, col1:col2]
#         def_min = np.min(deformation_array_extract)
#         def_max = np.max(deformation_array_extract)
#         topo_min = np.min(topo_array_extract)
#         topo_max = np.max(topo_array_extract)
#         deformation_array_extract_norm = (((deformation_array_extract - def_min)
#                                            / (def_max - def_min)) * width + start)
#         topo_array_extract_norm = (((topo_array_extract - topo_min)
#                                     / (topo_max - topo_min)) * width + start)
#         end_coord, start_coord = (general_functions.nearest_bound_coord(lon2, base, 'up'),
#                                   general_functions.nearest_bound_coord(lon1, base, 'down'))
#     return deformation_array_extract_norm, topo_array_extract_norm, index_, round(start_coord, 1), round(end_coord, 1)

# def make_cross_section_line(path_2_defo, deformation_array, topo_array, point_lons, point_lats):
#     lon1, lat1 = point_lons[0], point_lats[0]
#     lon2, lat2 = point_lons[1], point_lats[1]
#     start = -1
#     end = 1
#     width = end - start
#     base = 0.1
#     if lon1 == lon2:
#         if lat1 > lat2:
#             lat_1_temp = lat2
#             lat_2_temp = lat1
#             lat1 = lat_1_temp
#             lat2 = lat_2_temp
#         col1, row1 = rasterio_basic_functions.lon_lat_2_col_row(path_2_defo, lon1, lat1)
#         col2, row2 = rasterio_basic_functions.lon_lat_2_col_row(path_2_defo, lon2, lat2)
#         col2 = col1
#         if row1 > row2:
#             temp_row_1 = row2
#             temp_row_2 = row1
#             row1 = temp_row_1
#             row2 = temp_row_2
#         index_ = range(row1, row2)
#         deformation_array_extract = deformation_array[row1:row2, col2]
#         topo_array_extract = topo_array[row1:row2, col2]
#         def_min = np.min(deformation_array_extract)
#         def_max = np.max(deformation_array_extract)
#         topo_min = np.min(topo_array_extract)
#         topo_max = np.max(topo_array_extract)
#         deformation_array_extract_norm = (((deformation_array_extract - def_min)
#                                            / (def_max - def_min)) * width + start)
#         topo_array_extract_norm = (((topo_array_extract - topo_min)
#                                     / (topo_max - topo_min)) * width + start)
#         start_coord, end_coord = (general_functions.nearest_bound_coord(lat2, base, 'up'),
#                                   general_functions.nearest_bound_coord(lat1, base, 'down'))
#     else:
#         if lon1 > lon2:
#             lon_1_temp = lon2
#             lon_2_temp = lon1
#             lon1 = lon_1_temp
#             lon2 = lon_2_temp
#         col1, row1 = rasterio_basic_functions.lon_lat_2_col_row(path_2_defo, lon1, lat1)
#         col2, row2 = rasterio_basic_functions.lon_lat_2_col_row(path_2_defo, lon2, lat2)
#         row2 = row1
#         if col1 > col2:
#             temp_col_1 = col2
#             temp_col_2 = col1
#             col1 = temp_col_1
#             col2 = temp_col_2
#         index_ = range(col1, col2)
#         deformation_array_extract = deformation_array[row2, col1:col2]
#         topo_array_extract = topo_array[row2, col1:col2]
#         def_min = np.min(deformation_array_extract)
#         def_max = np.max(deformation_array_extract)
#         topo_min = np.min(topo_array_extract)
#         topo_max = np.max(topo_array_extract)
#         deformation_array_extract_norm = (((deformation_array_extract - def_min)
#                                            / (def_max - def_min)) * width + start)
#         topo_array_extract_norm = (((topo_array_extract - topo_min)
#                                     / (topo_max - topo_min)) * width + start)
#         start_coord, end_coord = (general_functions.nearest_bound_coord(lon2, base, 'up'),
#                                   general_functions.nearest_bound_coord(lon1, base, 'down'))
#     return deformation_array_extract_norm, topo_array_extract_norm, index_, round(start_coord, 1), round(end_coord, 1)


def mean_sliced_list(numpy_file_list, half_size, lon_cols, lat_rows, time_vector_month_shift, m_2_mm):
    mean_sliced_all, coeff_all = [], []
    for lon_col, lat_row in zip(lon_cols, lat_rows):
        mean_sliced = []
        for numpy_file in numpy_file_list:
            arr_open = general_functions.open_numpy_file(numpy_file)
            row_start = lat_row - half_size
            row_end = lat_row + half_size + 1
            col_start = lon_col - half_size
            col_end = lon_col + half_size + 1
            arr_slice = arr_open[row_start:row_end, col_start:col_end]
            arr_mean = np.mean(arr_slice) * m_2_mm
            mean_sliced.append(arr_mean)
        coeff = np.polyfit(time_vector_month_shift, mean_sliced, 1)
        mean_sliced_all.append(mean_sliced)
        coeff_all.append(coeff)
    return mean_sliced_all, coeff_all


def mean_sliced_list_no_poly(numpy_file_list, half_size, lon_cols, lat_rows, m_2_mm):
    mean_sliced_all = []
    for lon_col, lat_row in zip(lon_cols, lat_rows):
        mean_sliced = []
        for numpy_file in numpy_file_list:
            arr_open = general_functions.open_numpy_file(numpy_file)
            row_start = lat_row - half_size
            row_end = lat_row + half_size + 1
            col_start = lon_col - half_size
            col_end = lon_col + half_size + 1
            arr_slice = arr_open[row_start:row_end, col_start:col_end]
            arr_mean = np.mean(arr_slice) * m_2_mm
            mean_sliced.append(arr_mean)
        mean_sliced_all.append(mean_sliced)
    return mean_sliced_all


def find_average_value_at_point(base_array, lon_, lat_, half_size_, base_tif_path):
    col, row = rasterio_basic_functions.lon_lat_2_col_row(base_tif_path, lon_, lat_)
    row_start = row - half_size_
    row_end = row + half_size_ + 1
    col_start = col - half_size_
    col_end = col + half_size_ + 1
    arr_slice = base_array[row_start:row_end, col_start:col_end]
    arr_mean = float(np.nanmean(arr_slice))
    return round(arr_mean, 5)


def read_df_gps(txt_file_name):
    df_open = pd.read_csv(txt_file_name, sep='\t', skiprows=36, header=None)
    return df_open


def test_numeric(num_):
    try:
        float(num_)
        isnum = True
    except ValueError:
        isnum = False
    return isnum


def toYearFraction(date):
    def sinceEpoch(date):  # returns seconds since epoch
        return time.mktime(date.timetuple())

    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year + 1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed / yearDuration

    return date.year + fraction


def convert_YYYY_MM_DD_HH_MM_SS_2_YYYY_YYY(YYYY_MM_DD, HH_MM_SS):
    YYYY, MM_, DD = int(str(int(YYYY_MM_DD))[0:4]), int(str(int(YYYY_MM_DD))[4:6]), int(str(int(YYYY_MM_DD))[6:])
    HH, MM, SS = int(str(int(HH_MM_SS))[0:2]), int(str(int(HH_MM_SS))[2:4]), int(str(int(HH_MM_SS))[4:])
    date_ = dt(YYYY, MM_, DD, HH, MM, SS)
    return date_


def parse_gps_txt_file(gps_path_):
    test_df = read_df_gps(gps_path_[0])
    df_parsed_store = []
    for index, row in test_df.iterrows():
        parsed_rows = [x for x in row[0].split(' ') if x != '']
        df_parsed_store.append(parsed_rows)
    column_names = df_parsed_store[0]
    data_ = df_parsed_store[1:]
    data_float = []
    for x_ in data_:
        float_x_ = [float(x) if test_numeric(x) else str(x) for x in x_]
        data_float.append(float_x_)
    df_final = pd.DataFrame(data_float, columns=column_names)
    date_s = list(map(toYearFraction, list(
        map(convert_YYYY_MM_DD_HH_MM_SS_2_YYYY_YYY, list(df_final.iloc[:, 0]), list(df_final.iloc[:, 1])))))
    df_final['YYYY_YYY'] = date_s
    return df_final


def read_df_camp_gps(txt_file_name):
    df_open = pd.read_excel(txt_file_name)
    return df_open


class combine_all_tifs():
    def __init__(self, linear_or_coh, tif_name, geol_overlay, topo_box_coords, ndvi_overlay, composite, plot_points,
                 make_hist, saturate_linear, topo_correlation, find_statistics, base_map_make, make_time_series,
                 make_examples):
        self.base_path_time_series = None
        self.kwargs_transform_path = None
        self.binary_spatial_array_path = None
        self.tsout_path_months = None
        self.tsout_path_all = None
        self.time_series_path = None
        self.time_series_save_name = None
        self.to_mm = None
        self.half_size = 5
        kwargs_open = None
        self.transform_list = None
        tsout_open_month = None
        tsout_open_all = None
        self.times_month = None
        self.times_month_shifted = None
        self.times_all = None
        self.times_all_shifted = None
        msbas_used_np_month = None
        msbas_used_np_all = None
        test_nps_months = None
        test_nps_all = None
        self.binary_spatial_arrays_months = None
        self.binary_spatial_arrays_all = None
        self.def_months = None
        self.coeff_months = None
        self.def_all = None
        self.coeff_all = None
        self.ylim_min = None
        self.ylim_max = None
        self.base_path_time_series_1 = None
        self.base_path_time_series_2 = None

        self.kwargs_transform_path_1 = None
        self.binary_spatial_array_path_1 = None
        self.tsout_path_all_1 = None
        kwargs_open_1 = None
        self.transform_list_1 = None
        tsout_open_all_1 = None
        self.times_all_1 = None
        self.times_all_shifted_1 = None
        msbas_used_np_all_1 = None
        test_nps_all_1 = None
        self.binary_spatial_arrays_all_1 = None

        self.kwargs_transform_path_2 = None
        self.binary_spatial_array_path_2 = None
        self.tsout_path_all_2 = None
        kwargs_open_2 = None
        self.transform_list_2 = None
        tsout_open_all_2 = None
        self.times_all_2 = None
        self.times_all_shifted_2 = None
        msbas_used_np_all_2 = None
        test_nps_all_2 = None
        self.binary_spatial_arrays_all_2 = None

        self.times_all_use = None
        self.times_all_shifted_use = None

        self.base_path = Path.cwd().absolute().parent
        self.plot_points = plot_points
        self.make_time_series = make_time_series
        self.dem_path = self.base_path / 'DEM_NDVI_data'
        self.base_kwargs_paths = [Path(str(x)) for x in
                                  self.base_path.glob(
                                      'Pair_P*_F*_P*_F???/Plots_Regression/tif_meta_data/kwargs_tif.json')]
        self.base_kwargs_paths.sort()
        self.ndvi_overlay = ndvi_overlay
        self.time_series_points_path = self.base_path / 'Combined_Pairs/Time_Series_Points/Plot_TIme_Series_points.txt'
        self.make_hist = make_hist
        self.composite = composite
        self.return_statistics = find_statistics
        self.base_map_make = base_map_make
        self.saturate_linear = saturate_linear
        self.make_topo_correlation = topo_correlation
        self.gps_time_series_path = self.base_path / 'Combined_Pairs/GPS_Time_Series'
        self.gps_campaign_path = self.base_path / 'Combined_Pairs/Campaign_GPS/SoB_GAMITsubsidence.xls'
        self.gps_named_all = self.base_path / 'Combined_Pairs/Campaign_GPS/gps_Rates_edited.xls'
        camp_open = read_df_camp_gps(self.gps_campaign_path)
        gps_named_all_open = read_df_camp_gps(self.gps_named_all)
        self.gps_time_series_paths = ([Path(str(x)) for x in self.gps_time_series_path.glob('*.eos.final.itr14.pos')])
        self.combine_pairs_path = self.base_path / 'Combined_Pairs/Path_114_12_referenced'
        self.referenced_pairs_save_path = self.combine_pairs_path / 'referenced_path_frame'
        self.combined_tifs_path = self.base_path / 'Combined_Pairs/Path_114_12_combined'
        self.path_2_sentinel_composite = self.dem_path / 'sentinel_composite_geometry_cropped.tif'
        self.ndvi_path = self.dem_path / 'ndvi_crop_resampled.tif'
        self.path_2_kmls = self.base_path / 'make_base_map'
        self.path_2_topo_ndvi = self.base_path / 'Topo_NDVI_plots_and_code'
        self.path_2_dem = self.dem_path / 'combined_dem.tif'
        self.sundarbans_shp = self.base_path / 'Combined_Pairs/Sundarbans/Sundarbans_2015poly.shp'
        self.json_time_series_paths = self.base_path / 'Combined_Pairs/Path_114_12_combined/Time_Series_Plots'
        self.json_time_series_files = None
        self.combined_linear_path_dem_crop = self.dem_path / 'linear_rates_combined_all_dem_crop.tif'
        self.kml_paths = [Path(str(x)) for x in self.path_2_kmls.glob('*.kml')]
        self.corners = return_frame_coords(self.kml_paths)
        self.make_examples = make_examples
        self.linear_rate_tif_paths = ([Path(str(x)) for x in
                                       self.referenced_pairs_save_path.glob('linear*.tif')])
        self.coh_paths = ([Path(str(x)) for x in self.combine_pairs_path.glob('coherence_mask_*.tif')])
        self.linear_rates_tif_combined_path_save = self.combined_tifs_path / 'linear_rates_combined_all.tif'
        self.coh_tif_combined_path_save = self.combined_tifs_path / 'coherence_combined_all.tif'
        self.geologic_map_path = str(self.combine_pairs_path / 'geo8bg/geo8bg.shp')
        self.district_boundaries_path = (self.base_path /
                                         'Combined_Pairs/District_Boundaries/stanford-ps807dh8348-geojson.json')
        self.combined_list = self.linear_rate_tif_paths + self.coh_paths
        self.example_ifg_path = (self.base_path /
                                 'Pair_P114_F76_P150_F510/Path_114_Frame_76/Main/Processed_Data/20210109_20210121')
        self.example_ifg_wrapped_path = self.example_ifg_path / '20210109_20210121_w_phase_VV.tif'
        self.example_ifg_unwrapped_path = self.example_ifg_path / '20210109_20210121_unw_phase_VV.tif'
        self.example_coh_path = self.example_ifg_path / '20210109_20210121_coherence_VV.tif'
        self.example_non_gacos_path = self.example_ifg_path / '20210109_20210121_displacement_VV_rcut_mcoh.tif'
        self.example_gacos_path = self.example_ifg_path / '20210109_20210121_displacement_VV_rcut_mcoh_gacos.tif'
        self.example_ifg_masked_save_path = self.combined_tifs_path / 'ifg_example.png'
        self.example_non_gacos_gacos_save_path = self.combined_tifs_path / 'gacos_non_gacos_example.png'
        self.gps_sinusoid_name = self.combined_tifs_path / 'gps_defo_sinusoid.png'
        self.points_of_interest_sinusoid_name = self.combined_tifs_path / 'pt_interest_sinusoid.png'
        self.topo_box_coords = topo_box_coords
        if self.make_topo_correlation == 'true':
            self.topo_correlation_save_name = self.combined_tifs_path / 'topo_correlation_plot_all.png'
            self.topo_correlation_save_name_all = self.combined_tifs_path / 'topo_cross_section_all.png'
            self.topo_correlation_save_name_1 = self.combined_tifs_path / 'topo_cross_section_1.png'
            self.topo_correlation_save_name_2 = self.combined_tifs_path / 'topo_cross_section_2.png'
            self.topo_correlation_save_name_3 = self.combined_tifs_path / 'topo_cross_section_3.png'
            self.topo_correlation_save_name_4 = self.combined_tifs_path / 'topo_cross_section_4.png'
            self.topo_correlation_save_name_5 = self.combined_tifs_path / 'topo_cross_section_5.png'
            self.topo_correlation_save_name_6 = self.combined_tifs_path / 'topo_cross_section_6.png'
            self.topo_correlation_save_name_7 = self.combined_tifs_path / 'topo_cross_section_7.png'
            self.topo_correlation_save_name_8 = self.combined_tifs_path / 'topo_cross_section_8.png'
            self.topo_correlation_save_name_9 = self.combined_tifs_path / 'topo_cross_section_9.png'
        if self.return_statistics == 'true':
            self.statistics_path = self.combined_tifs_path / "relevant_statistics.json"
            self.gps_def_comparison = self.combined_tifs_path / "gps_def_comparison.csv"
        if not self.topo_box_coords == 'box_empty':
            self.topo_box_coords_vals = np.array([float(x) for x in self.topo_box_coords.split(',')])
            total_num_coords = len(self.topo_box_coords_vals)
            left_lons_ints, top_lats_ints, right_lons_ints, bot_lats_ints = (np.arange(0, total_num_coords, 4),
                                                                             np.arange(1, total_num_coords, 4),
                                                                             np.arange(2, total_num_coords, 4),
                                                                             np.arange(3, total_num_coords, 4))
            self.left_lons, self.top_lats, self.right_lons, self.bot_lats = (list(self.topo_box_coords_vals
                                                                                  [left_lons_ints]),
                                                                             list(self.topo_box_coords_vals
                                                                                  [top_lats_ints]),
                                                                             list(self.topo_box_coords_vals
                                                                                  [right_lons_ints]),
                                                                             list(self.topo_box_coords_vals
                                                                                  [bot_lats_ints]))

        self.geol_overlay = geol_overlay
        self.linear_rate_arrays_use = None
        self.coh_arrays_use = None
        self.master_width = None
        self.master_height = None
        self.master_affine = None
        self.new_kwargs = None
        self.points_lon_campaign = [float(x) for x in list(camp_open['Long'])]
        self.points_lat_campaign = [float(x) for x in list(camp_open['Lat'])]
        self.rates_campaign = [float(x) for x in list(camp_open['Vertical+S  (mm/y)'])]
        self.names_campaign = [str(x) for x in list(camp_open['Sta ID'])]
        self.points_lon_named_all = [float(x) for x in list(gps_named_all_open['Lon'])]
        self.points_lat_named_all = [float(x) for x in list(gps_named_all_open['Lat'])]
        self.rates_named_all = [float(x) for x in list(gps_named_all_open['Vu[mm/yr]'])]
        self.names_named_all = [str(x) for x in list(gps_named_all_open['Name'])]
        self.points_gps_lon = [90.401264,
                               89.45940429,
                               89.43590233,
                               89.53130045]
        self.points_gps_lat = [23.726679,
                               21.81640899,
                               22.50457227,
                               22.85263818]
        # self.points_ts_lon = [89.45940429, 89.43590233, 89.53130045, 87.48962748, 88.51597, 90.355961,
        #                       90.585, 89.424447, 90.17150447]
        # self.points_ts_lat = [21.81640899, 22.50457227, 22.85263818, 22.96315400, 23.01233, 22.2598984,
        #                       24.363917, 24.526216, 23.784344]
        self.points_ts_lon = [89.45940429,
                              89.43590233,
                              89.53130045
                              ]
        self.points_ts_lat = [21.81640899,
                              22.50457227,
                              22.85263818
                              ]
        points_list_read_file_df = general_functions.read_df_tab(self.time_series_points_path)
        points_ts_lon_file_list = list(points_list_read_file_df['Point Longitude'])
        points_ts_lat_file_list = list(points_list_read_file_df['Point_lattitude'])
        self.points_ts_lon_all = self.points_ts_lon + points_ts_lon_file_list
        self.points_ts_lat_all = self.points_ts_lat + points_ts_lat_file_list
        self.gps_labels = ['DHAK', 'HRNP', 'PD32', 'KLNA']
        self.time_series_labels = ['HRNP', 'PD32', 'KLNA']
        self.gps_rates_3 = [-7.42, -4.87, -4.32]
        self.gps_rates_3_names = ['HRNP_gps', 'PD32_gps', 'KLNA_gps']
        time_series_names = list(points_list_read_file_df['Point Name'])
        self.time_series_labels_all = self.time_series_labels + time_series_names
        self.polygon_point_list_lons = [87.13, 89.1, 88.87, 89, 91.47, 91.29, 90.8, 91.02, 90.74, 88.33, 88.23,
                                        87.45, 87.13]
        self.polygon_point_list_lats = [23.39, 23.74, 24.82, 25.38, 25.06, 24.13, 24.02, 22.89, 21.9, 21.48,
                                        22.04, 21.89, 23.39]
        self.cross_section_labels_start = [('C' + str(x)) for x in range(1, 10)]
        self.cross_section_labels_end = [('C' + str(x) + '`') for x in range(1, 10)]
        self.cross_section_lon_list = [[87.47272941, 87.47272941],
                                       [87.27835342, 91],
                                       [88.16137596, 88.16137596],
                                       [89.87743883, 89.87743883],
                                       [89.59975880, 90.03849321],
                                       [90.59385351, 90.59385351],
                                       [89, 89.75],
                                       [88.2, 89.8],
                                       [89.4, 89.4]]
        self.cross_section_lat_list = [[21.94213776, 23.3],
                                       [23.02, 23.02],
                                       [22.63078425, 23.08062608],
                                       [24.44125841, 24.96329690],
                                       [24.68561684, 24.68561684],
                                       [24.16357826, 24.56343760],
                                       [21.95, 21.95],
                                       [22.44, 22.44],
                                       [22.3, 23.3]]
        self.fault_lat_lons = [91.541820181714, 25.13096291590124,
                               91.14621369211464, 24.87031271464589,
                               90.74536635179479, 24.51171298730961,
                               90.45242609456923, 23.86854051228257,
                               90.4618093587257, 22.92105347595646,
                               90.78328050366275, 21.89515333662549,
                               91.38613397334895, 20.60155021309592,
                               92.10287121238154, 19.87936379969205,
                               92.40280296854121, 19.5748775834888,
                               92.84217626887933, 19.26861740191356,
                               93.05224293761468, 18.91060951407074,
                               93.31974321514502, 18.66671357380149,
                               93.31286340111504, 18.3275256879741,
                               93.4384938313968, 18.06597525285159,
                               93.69831068623874, 17.96349910647801,
                               93.9484833580179, 17.86042512651765,
                               94.00792372874173, 17.58898020169614,
                               94.02818037813435, 17.23537176013848,
                               93.90036752037889, 16.94657800982252]
        i = 0
        j = 1
        lons = []
        lats = []
        len_ = len(self.fault_lat_lons)
        while j <= len_:
            lon_, lat_ = self.fault_lat_lons[i], self.fault_lat_lons[j]
            i = i + 2
            j = j + 2
            lons.append(lon_)
            lats.append(lat_)
        self.fault_lons = lons
        self.fault_lats = lats
        self.linear_or_coh = linear_or_coh
        self.num_of_std = 2
        self.base_linear = 5
        self.base_coh = 0.1
        if self.linear_or_coh == 'linear':
            self.paths_to_use = self.combined_tifs_path / tif_name
            tif_name_split = tif_name.split('.')[0]
            if self.geol_overlay == 'true':
                self.path_save_name = self.combined_tifs_path / (tif_name_split + '_geo_overlay' + '.png')
            elif not self.topo_box_coords == 'box_empty':
                self.path_save_name = self.combined_tifs_path / (tif_name_split + '_topo_overlay' + '.png')
            elif self.ndvi_overlay == 'true':
                self.path_save_name = self.combined_tifs_path / (tif_name_split + '_ndvi_overlay' + '.png')
            elif self.composite == 'true':
                self.path_save_name = self.combined_tifs_path / (tif_name_split + '_composite' + '.png')
            elif self.base_map_make == 'true':
                self.path_save_name = self.combined_tifs_path / (tif_name_split + '_base_map' + '.png')
            elif self.saturate_linear == 'true':
                self.path_save_name = self.combined_tifs_path / (tif_name_split + '_saturate' + '.png')
            else:
                self.path_save_name = self.combined_tifs_path / (tif_name_split + '.png')
            self.path_save_name_hist = self.combined_tifs_path / (tif_name_split + '_hist.png')
            self.path_save_name_hist_m_std = self.combined_tifs_path / (tif_name_split + '_hist_m_std.png')
        else:
            self.path_to_coh = self.combined_tifs_path
            self.paths_to_use = self.path_to_coh / tif_name
            tif_name_split = tif_name.split('.')[0]
            if self.geol_overlay == 'true':
                self.path_save_name = self.combined_tifs_path / (tif_name_split + '_geo_overlay' + '.png')
            elif not self.topo_box_coords == 'box_empty':
                self.path_save_name = self.combined_tifs_path / (tif_name_split + '_topo_overlay' + '.png')
            elif self.ndvi_overlay == 'true':
                self.path_save_name = self.combined_tifs_path / (tif_name_split + '_ndvi_overlay' + '.png')
            else:
                self.path_save_name = self.combined_tifs_path / (tif_name_split + '.png')
            self.path_save_name_hist = self.combined_tifs_path / (tif_name_split + '_hist.png')
            self.path_save_name_hist_m_std = self.combined_tifs_path / (tif_name_split + '_hist_m_std.png')
        self.save_name_title = tif_name_split
        print('linear_or_coh=', self.linear_or_coh,
              'geol_overlay=', self.geol_overlay,
              'topo_box_plot=', self.topo_box_coords,
              'ndvi_overlay=', self.ndvi_overlay,
              'composite_make=', self.composite,
              'self_plot_points=', self.plot_points,
              'make_hists=', self.make_hist,
              'saturate_linear=', self.saturate_linear,
              'make_topo_correlations=', self.make_topo_correlation,
              'make_stats=', self.return_statistics,
              'make_base_map=', self.base_map_make,
              'make_time_series=', self.make_time_series)

    def define_master_resample_parameters(self):
        longitude_spacing_list, upper_left_corner_longitudes_list = [], []
        lattitude_spacing_list, upper_left_corner_lattitude_list = [], []
        lower_right_corner_longitude_list = []
        lower_right_corner_lattitude_list = []

        for asc_dsc in self.combined_list:
            affine_transform = rasterio_basic_functions.get_tif_transform(asc_dsc)
            tif_height, tif_width = (rasterio_basic_functions.get_tif_height(asc_dsc),
                                     rasterio_basic_functions.get_tif_width(asc_dsc))
            lower_right_corner_longitude = (affine_transform * (tif_width, tif_height))[0]
            lower_right_corner_lattitude = (affine_transform * (tif_width, tif_height))[1]
            longitude_spacing, upper_left_corner_longitude = affine_transform.a, affine_transform.c
            lattitude_spacing, upper_left_corner_lattitude = affine_transform.e, affine_transform.f
            upper_left_corner_longitudes_list.append(upper_left_corner_longitude)
            upper_left_corner_lattitude_list.append(upper_left_corner_lattitude)
            lower_right_corner_longitude_list.append(lower_right_corner_longitude)
            lower_right_corner_lattitude_list.append(lower_right_corner_lattitude)
            longitude_spacing_list.append(longitude_spacing)
            lattitude_spacing_list.append(lattitude_spacing)

        max_longitude = max(lower_right_corner_longitude_list) + 0.05
        min_longitude = min(upper_left_corner_longitudes_list) - 0.05
        max_lattitude = max(upper_left_corner_lattitude_list) + 0.05
        min_lattitude = min(lower_right_corner_lattitude_list) - 0.05

        mean_longitude_spacing = mean(longitude_spacing_list)
        mean_lattitude_spacing = mean(lattitude_spacing_list)

        self.master_width = math.ceil((abs(max_longitude - min_longitude))
                                      / abs(mean_longitude_spacing))
        self.master_height = math.ceil((abs(max_lattitude - min_lattitude))
                                       / abs(mean_lattitude_spacing))
        self.master_affine = rasterio_basic_functions.make_affine(mean_longitude_spacing,
                                                                  0.0,
                                                                  min_longitude,
                                                                  0.0,
                                                                  mean_lattitude_spacing,
                                                                  max_lattitude)
        ii = 0
        linear_rates_save = []
        coh_save = []
        for linear_rate_tif, coh_tif in zip(self.linear_rate_tif_paths, self.coh_paths):
            linear_rate_array_old = rasterio_basic_functions.tif_2_array(linear_rate_tif)
            coh_array_old = rasterio_basic_functions.tif_2_array(coh_tif)
            linear_transform_old = rasterio_basic_functions.get_tif_transform(linear_rate_tif)
            coh_transform_old = rasterio_basic_functions.get_tif_transform(coh_tif)
            base_crs = rasterio_basic_functions.get_tif_crs(linear_rate_tif)
            base_kwargs = rasterio_basic_functions.get_kwargs(linear_rate_tif)
            base_kwargs.update({'nodata': 0.,
                                'width': self.master_width,
                                'height': self.master_height,
                                'transform': self.master_affine
                                })
            if ii == 0:
                self.new_kwargs = base_kwargs
            ii = ii + 1
            linear_rate_reprojected_array = rasterio_basic_functions.reproject_tif_array(linear_rate_tif,
                                                                                         self.master_height,
                                                                                         self.master_width,
                                                                                         linear_rate_array_old,
                                                                                         linear_transform_old,
                                                                                         base_crs,
                                                                                         self.master_affine)
            coherence_reprojected_array = rasterio_basic_functions.reproject_tif_array(coh_tif,
                                                                                       self.master_height,
                                                                                       self.master_width,
                                                                                       coh_array_old,
                                                                                       coh_transform_old,
                                                                                       base_crs,
                                                                                       self.master_affine)
            linear_rates_save.append(linear_rate_reprojected_array)
            coh_save.append(coherence_reprojected_array)
        self.linear_rate_arrays_use = linear_rates_save
        self.coh_arrays_use = coh_save

    def combine_arrays(self):
        base_linear_array = np.zeros((self.master_height, self.master_width))
        base_linear_divide_array = np.zeros((self.master_height, self.master_width))
        for linear_rate_array in self.linear_rate_arrays_use:
            base_linear_array = base_linear_array + linear_rate_array
            index_array = np.where(linear_rate_array != 0,
                                   1,
                                   0)
            base_linear_divide_array = base_linear_divide_array + index_array
        final_linear_rate_array = np.nan_to_num(base_linear_array / base_linear_divide_array).astype('float32')

        base_coh_array = np.zeros((self.master_height, self.master_width))
        base_index_coh_array = np.zeros((self.master_height, self.master_width))
        for coh_tif in self.coh_arrays_use:
            coh_mask = np.where(coh_tif < 0.15,
                                0,
                                1)
            coh_array_use = np.multiply(coh_mask, coh_tif)
            base_coh_array = base_coh_array + coh_array_use
            index_array = np.where(coh_array_use != 0,
                                   1,
                                   0)
            base_index_coh_array = base_index_coh_array + index_array
        final_coh_array = np.nan_to_num(base_coh_array / base_index_coh_array).astype('float32')

        rasterio_basic_functions.write_reprojected_array_2_tif(final_linear_rate_array,
                                                               self.linear_rates_tif_combined_path_save,
                                                               self.new_kwargs)

        rasterio_basic_functions.write_reprojected_array_2_tif(final_coh_array,
                                                               self.coh_tif_combined_path_save,
                                                               self.new_kwargs)

    def make_plots(self):
        base_array = rasterio_basic_functions.tif_2_array(self.paths_to_use)
        base_array[base_array == 0.] = np.nan
        if self.ndvi_overlay == 'true':
            base_ndvi_array = rasterio_basic_functions.tif_2_array(self.ndvi_path)
            base_ndvi_array[base_ndvi_array == 0.] = np.nan
            mean_tif_ndvi = np.nanmean(base_ndvi_array)
            std_tif_ndvi = np.nanstd(base_ndvi_array)
            lower_bound_u_z_ndvi, upper_bound_u_z_ndvi = general_functions.set_bounds(mean_tif_ndvi,
                                                                                      std_tif_ndvi,
                                                                                      self.num_of_std,
                                                                                      0.2)
            cmap_ndvi = mpl.cm.RdYlGn
            norm_ndvi = mpl.colors.Normalize(vmin=lower_bound_u_z_ndvi, vmax=upper_bound_u_z_ndvi)

        red_band = rasterio_basic_functions.tif_2_array_multiband(self.path_2_sentinel_composite, 1)
        green_band = rasterio_basic_functions.tif_2_array_multiband(self.path_2_sentinel_composite, 2)
        blue_band = rasterio_basic_functions.tif_2_array_multiband(self.path_2_sentinel_composite, 3)
        left_, bottom_, right_, top_ = rasterio_basic_functions.return_bounds(self.path_2_sentinel_composite)
        red_band_norm = normalize(red_band, 2)
        green_band_norm = normalize(green_band, 2)
        blue_band_norm = normalize(blue_band, 2)
        rgb_stack = np.dstack((red_band_norm, green_band_norm, blue_band_norm))

        if self.linear_or_coh == 'linear':
            base_array = base_array * 1000
            mean_tif = np.nanmean(base_array)
            std_tif = np.nanstd(base_array)
            cmap = mpl.cm.rainbow
            if self.saturate_linear == 'false':
                lower_bound_u_z, upper_bound_u_z = general_functions.set_bounds(mean_tif,
                                                                                std_tif,
                                                                                self.num_of_std,
                                                                                self.base_linear)
            else:
                base_ = 1
                lower_b = mean_tif - std_tif
                upper_b = mean_tif + std_tif
                lower_b_u = general_functions.nearest_bound_coord(lower_b, base_, 'down')
                upper_b_u = general_functions.nearest_bound_coord(upper_b, base_, 'up')
                lower_bound_u_z, upper_bound_u_z = lower_b_u, upper_b_u
        else:
            cmap = mpl.cm.Blues.reversed()
            min_coh = np.nanmin(base_array)
            max_coh = np.nanmax(base_array)
            mean_tif = np.nanmean(base_array)
            std_tif = np.nanstd(base_array)
            # lower_bound_u_z, upper_bound_u_z = (general_functions.nearest_bound_coord(min_coh, self.base_coh, 'down'),
            #                                     general_functions.nearest_bound_coord(max_coh, self.base_coh, 'up'))
            #
            lower_bound_u_z, upper_bound_u_z = general_functions.set_bounds_coh(mean_tif,
                                                                                std_tif,
                                                                                self.num_of_std,
                                                                                self.base_coh)
            if lower_bound_u_z < 0:
                lower_bound_u_z = 0
            if upper_bound_u_z > 1:
                upper_bound_u_z = 1

        norm = mpl.colors.Normalize(vmin=lower_bound_u_z, vmax=upper_bound_u_z)
        left, bottom, right, top = rasterio_basic_functions.return_bounds(self.paths_to_use)
        lon_spacing = abs(right - left) / 5
        lat_spacing = abs(top - bottom) / 5
        if 0.25 < lon_spacing < 0.5:
            lon_spacing = 0.25
        elif 0.5 < lon_spacing < 1:
            lon_spacing = 0.5
        else:
            lon_spacing = 0.1
        if 0.25 < lat_spacing < 0.5:
            lat_spacing = 0.25
        elif 0.5 < lat_spacing < 1:
            lat_spacing = 0.5
        else:
            lat_spacing = 0.1
        base_degree_lon = lon_spacing
        base_degree_lat = lat_spacing
        left_use, bottom_use, right_use, top_use = (general_functions.nearest_bound_coord(left,
                                                                                          base_degree_lat, 'down'),
                                                    general_functions.nearest_bound_coord(bottom,
                                                                                          base_degree_lon, 'down'),
                                                    general_functions.nearest_bound_coord(right,
                                                                                          base_degree_lat, 'up'),
                                                    general_functions.nearest_bound_coord(top,
                                                                                          base_degree_lon, 'up'))
        lon_mid = (left + right) / 2
        patterns = ["/", "....", ".", "**", "+", "x", "+/", "\\\\", "|", "////", '+\\', '||', '--']
        frame, cmap_, vmin, vmax = read_geologic_map_shp(self.geologic_map_path)
        # fig_size_col, fig_size_row = 20, 22
        fig_size_col, fig_size_row = 20, 22
        if self.ndvi_overlay == 'true':
            fig, (axs1, axs2) = plt.subplots(1, 2, figsize=(fig_size_col, fig_size_row), facecolor='w', edgecolor='k')
            m1 = Basemap(epsg=4326, llcrnrlat=bottom, urcrnrlat=top,
                         llcrnrlon=left, urcrnrlon=right, resolution='f', lon_0=lon_mid, ax=axs1)
            m2 = Basemap(epsg=4326, llcrnrlat=bottom, urcrnrlat=top,
                         llcrnrlon=left, urcrnrlon=right, resolution='f', lon_0=lon_mid, ax=axs2)
            # m.arcgisimage(service='ESRI_Imagery_World_2D', xpixels=2500, verbose=True, dpi=300)
            m1.imshow(rgb_stack, extent=[left_, right_, bottom_, top_], origin="upper", interpolation=None)
            m2.imshow(rgb_stack, extent=[left_, right_, bottom_, top_], origin="upper", interpolation=None)
            im1 = m1.imshow(base_ndvi_array, interpolation=None,
                            vmin=lower_bound_u_z_ndvi, vmax=upper_bound_u_z_ndvi, cmap=cmap_ndvi,
                            norm=norm_ndvi,
                            origin="upper", extent=[left_, right_, bottom_, top_])
            im2 = m2.imshow(base_array, extent=[left_, right_, bottom_, top_], interpolation=None,
                            vmin=lower_bound_u_z, vmax=upper_bound_u_z, cmap=cmap, norm=norm, origin="upper")

        elif self.composite == 'true':
            fig, axs1 = plt.subplots(1, 1, figsize=(fig_size_col, fig_size_row), facecolor='w', edgecolor='k')
            m = Basemap(epsg=4326, llcrnrlat=bottom, urcrnrlat=top,
                        llcrnrlon=left, urcrnrlon=right, resolution='f', lon_0=lon_mid, ax=axs1)
            im1 = m.imshow(rgb_stack, origin="upper", interpolation=None,
                           extent=[left_, right_, bottom_, top_])
        elif self.base_map_make == 'true':
            fig, axs1 = plt.subplots(1, 1, figsize=(fig_size_col, fig_size_row), facecolor='w', edgecolor='k')
            m = Basemap(epsg=4326, llcrnrlat=bottom, urcrnrlat=top,
                        llcrnrlon=left, urcrnrlon=right, resolution='f', lon_0=lon_mid, ax=axs1)
            im1 = m.imshow(rgb_stack, origin="upper", interpolation=None,
                           extent=[left_, right_, bottom_, top_])
            df_places = gpd.read_file(str(self.district_boundaries_path))
            patches = []
            for index, row in df_places.iterrows():
                district_name, poly = row['name_1'], row['geometry']
                if poly.geom_type == 'Polygon':
                    mpoly = shapely.ops.transform(m, poly)
                    patches.append(PolygonPatch(mpoly, edgecolor='k', facecolor='gray', fill=True, lw=2.0))
                elif poly.geom_type == 'MultiPolygon':
                    for subpoly in poly:
                        mpoly = shapely.ops.transform(m, subpoly)
                        patches.append(PolygonPatch(mpoly, edgecolor='k', facecolor='gray', fill=True, lw=2.0))
            patches_s = []
            sundarbans_shp = gpd.read_file(str(self.sundarbans_shp))
            for index, row in sundarbans_shp.iterrows():
                poly = row['geometry']
                if poly.geom_type == 'Polygon':
                    mpoly = shapely.ops.transform(m, poly)
                    patches_s.append(PolygonPatch(mpoly, edgecolor='r', fill=False, lw=2.0))
                elif poly.geom_type == 'MultiPolygon':
                    for subpoly in poly:
                        mpoly = shapely.ops.transform(m, subpoly)
                        patches_s.append(PolygonPatch(mpoly, edgecolor='r', fill=False, lw=2.0))
            axs1.add_collection(PatchCollection(patches, match_original=True, alpha=0.6))
            axs1.add_collection(PatchCollection(patches_s, match_original=True, facecolor='none'))
            plt.rcParams.update({'text.color': "black",
                                 'font.size': 30})
            df_places['coords'] = df_places['geometry'].apply(lambda x: x.representative_point().coords[:])
            df_places['coords'] = [coords[0] for coords in df_places['coords']]
            bad_names = ['Chittagong', 'Rangpur', 'Sylhet']
            for idx, row in df_places.iterrows():
                name_1 = row['name_1']
                if not any(bad_name == name_1 for bad_name in bad_names):
                    if name_1 == 'Barisal':
                        coord_x, coord_y = row['coords']
                        coord_x_u, coord_y_u = coord_x + 0.15, coord_y + 0.2
                    elif name_1 == 'Dhaka':
                        coord_x, coord_y = row['coords']
                        coord_x_u, coord_y_u = coord_x - 0.1, coord_y
                    elif name_1 == 'Khulna':
                        coord_x, coord_y = row['coords']
                        coord_x_u, coord_y_u = coord_x, coord_y
                    elif name_1 == 'Rajshahi':
                        coord_x, coord_y = row['coords']
                        coord_x_u, coord_y_u = coord_x + 0.3, coord_y
                    else:
                        continue
                    plt.annotate(s=name_1, xy=(coord_x_u, coord_y_u),
                                 horizontalalignment='center')
            mpl.rcParams.update(mpl.rcParamsDefault)
            axs2 = inset_axes(axs1, "30%", "30%", loc="upper left")
            inmap = Basemap(epsg=4326, llcrnrlat=3, urcrnrlat=36,
                            llcrnrlon=66, urcrnrlon=96, resolution='f', lon_0=(96 - 66) / 2, ax=axs2)
            inmap.drawcoastlines()
            inmap.fillcontinents(color='grey', lake_color='white')
            inmap.drawmapboundary(fill_color='lightblue')
            inmap.drawcountries(color='black', linewidth=1)
            bx, by = inmap(m.boundarylons, m.boundarylats)
            xy = list(zip(bx, by))
            patches_frame_overlay = []
            patches_overlay = []
            for lon_, lat_ in zip(self.polygon_point_list_lons, self.polygon_point_list_lats):
                item_ = [lon_, lat_]
                patches_frame_overlay.append(item_)
            polygon_array = np.array(patches_frame_overlay)
            patches_overlay.append(Polygon_mpl(polygon_array))
            axs1.add_collection(PatchCollection(patches_overlay, facecolor='none', edgecolor='k', linewidths=3.0))
            mapboundary = Polygon_mpl(xy, edgecolor='k', linewidth=2, fill=False)
            inmap.ax.add_patch(mapboundary)
            m.plot(self.fault_lons, self.fault_lats, linewidth=2, color='r', linestyle='--')
        else:
            fig, axs1 = plt.subplots(1, 1, figsize=(fig_size_col, fig_size_row), facecolor='w', edgecolor='k')
            m = Basemap(epsg=4326, llcrnrlat=bottom, urcrnrlat=top,
                        llcrnrlon=left, urcrnrlon=right, resolution='f', lon_0=lon_mid, ax=axs1)
            m.imshow(rgb_stack, extent=[left_, right_, bottom_, top_], interpolation=None, origin="upper")
            im1 = m.imshow(base_array, extent=[left_, right_, bottom_, top_], interpolation=None,
                           vmin=lower_bound_u_z, vmax=upper_bound_u_z, cmap=cmap, norm=norm, origin="upper")
            color_ = 'blue'
            color_text = 'blue'
            if ((self.linear_or_coh == 'linear' and not self.geol_overlay == 'true') and
                    (self.linear_or_coh == 'linear' and not self.saturate_linear == 'true')):
                pts_ = []
                for lons_, lats_ in zip(self.cross_section_lon_list, self.cross_section_lat_list):
                    lon1, lon2 = lons_[0], lons_[1]
                    lat1, lat2 = lats_[0], lats_[1]
                    pt = [(lon1, lat1), (lon2, lat2)]
                    pts_.append(pt)
                axs1.add_collection(LineCollection(pts_, edgecolor=color_, linewidth=3))
                plt.rcParams.update({'text.color': color_text,
                                     'font.size': 20})
                for lons_, lats_, start_label, end_label in zip(self.cross_section_lon_list,
                                                                self.cross_section_lat_list,
                                                                self.cross_section_labels_start,
                                                                self.cross_section_labels_end):
                    lon1, lon2 = lons_[0], lons_[1]
                    lat1, lat2 = lats_[0], lats_[1]
                    if lon1 == lon2:
                        axs1.annotate(s=start_label, xy=(lon1 + 0.015, lat1 - 0.04))
                        axs1.annotate(s=end_label, xy=(lon2 + 0.015, lat2 - 0.03))
                    else:
                        axs1.annotate(s=start_label, xy=(lon1 - 0.04, lat1 + 0.015))
                        axs1.annotate(s=end_label, xy=(lon2 - 0.03, lat2 + 0.015))
                mpl.rcParams.update(mpl.rcParamsDefault)

        if self.geol_overlay == 'true':
            for index, row in frame.iterrows():
                polygon_geom, geol_int = row[6], row[7]
                pattern = patterns[int(geol_int)]
                exterior_ring_coords = [polygon_geom.exterior.coords[:]]
                interior_rings_coords = [inner.coords[:] for inner in polygon_geom.interiors]
                interior_rings_coords_flat_list = []
                for sublist in interior_rings_coords:
                    for item_ in sublist:
                        interior_rings_coords_flat_list.append(item_)

                # number_of_exterior_ring = len(exterior_ring_coords)
                # number_of_inner_rings = len(interior_rings_coords)
                length_list_of_inner_rings = []
                length_list_of_exterion_rings = []

                for item in exterior_ring_coords:
                    len_ = len(item)
                    length_list_of_exterion_rings.append(len_)
                for item in interior_rings_coords:
                    len_ = len(item)
                    length_list_of_inner_rings.append(len_)

                length_list_of_all = length_list_of_exterion_rings + length_list_of_inner_rings

                len_intervals = [0]
                sum_ = 0.
                for len_number in length_list_of_all:
                    sum_ = sum_ + len_number
                    len_intervals.append(sum_)

                code_move = Path_m.MOVETO
                code_line_2 = Path_m.LINETO
                code_close_poly = Path_m.CLOSEPOLY
                all_rings = exterior_ring_coords[0] + interior_rings_coords_flat_list
                all_codes = []

                for len_of_poly in length_list_of_all:
                    idx = 0
                    code_list = []
                    while idx < len_of_poly:
                        if idx == 0:
                            code_list.append(code_move)
                        elif idx < (len_of_poly - 1):
                            code_list.append(code_line_2)
                        else:
                            code_list.append(code_close_poly)
                        idx = idx + 1
                    all_codes.append(code_list)

                code_flat_list = []
                for sublist in all_codes:
                    for item_ in sublist:
                        code_flat_list.append(item_)
                my_poly_path = Path_m(all_rings, code_flat_list)
                patch = PathPatch(my_poly_path, edgecolor='k', fill=False, hatch=pattern, lw=2.0)
                axs1.add_patch(patch)
            geol_units_new_ = list(frame.GLG.unique())
            legend_com = []
            for geol_, pat in zip(geol_units_new_, patterns):
                legen_c = mpatches.Patch(hatch=pat, label=geol_, edgecolor='k', facecolor='w', linewidth=3, fill=False)
                legend_com.append(legen_c)

        if not self.topo_box_coords == 'box_empty':
            patches = []
            for b_lat, t_lat, l_lon, r_lon in zip(self.bot_lats, self.top_lats, self.left_lons, self.right_lons):
                box_list_of_coords = np.array([[l_lon, t_lat], [r_lon, t_lat], [r_lon, b_lat], [l_lon, b_lat]])
                patches.append(Polygon_mpl(box_list_of_coords))
                del box_list_of_coords
            axs1.add_collection(PatchCollection(patches, facecolor='none', edgecolor='r', linewidths=3.0))

        points_lon_use, points_lat_use, labels_use = (self.points_ts_lon_all[3:], self.points_ts_lat_all[3:],
                                                      self.time_series_labels_all[3:])
        color_ = 'blue'
        color_text = 'blue'
        if self.plot_points == 'true':
            if self.ndvi_overlay == 'true':
                s1, s2 = 60, 30
                plt.rcParams.update({'text.color': color_text,
                                     'font.size': 14})
                for p_lon, p_lat, gps_label in zip(self.points_gps_lon, self.points_gps_lat, self.gps_labels):
                    m1.scatter(p_lon, p_lat, facecolor='none', edgecolor=color_, s=s1, marker='*', linewidths=3.0,
                               zorder=100)
                    axs1.annotate(s=gps_label, xy=(p_lon + 0.03, p_lat), horizontalalignment='left')
                for p_lon_ts, p_lat_ts, ts_label in zip(points_lon_use, points_lat_use,
                                                        labels_use):
                    m1.scatter(p_lon_ts, p_lat_ts, facecolor='none', edgecolor=color_, s=s2, marker='o', linewidths=3.0,
                               zorder=101)
                    axs1.annotate(s=ts_label, xy=(p_lon_ts + 0.03, p_lat_ts), horizontalalignment='left')
                for p_lon, p_lat, gps_label in zip(self.points_gps_lon, self.points_gps_lat, self.gps_labels):
                    m2.scatter(p_lon, p_lat, facecolor='none', edgecolor=color_, s=s1, marker='*', linewidths=3.0,
                               zorder=100)
                    axs2.annotate(s=gps_label, xy=(p_lon + 0.03, p_lat), horizontalalignment='left')
                for p_lon_ts, p_lat_ts, ts_label in zip(points_lon_use, points_lat_use,
                                                        labels_use):
                    m2.scatter(p_lon_ts, p_lat_ts, facecolor='none', edgecolor=color_, s=s2, marker='o', linewidths=3.0,
                               zorder=101)
                    axs2.annotate(s=ts_label, xy=(p_lon_ts + 0.03, p_lat_ts), horizontalalignment='left')
            else:
                s1, s2 = 100, 60
                if self.linear_or_coh == 'coh':
                    plt.rcParams.update({'text.color': color_text,
                                         'font.size': 18})
                    for p_lon, p_lat, gps_label in zip(self.points_gps_lon, self.points_gps_lat, self.gps_labels):
                        m.scatter(p_lon, p_lat, facecolor='none', edgecolor=color_, s=s1, marker='*', linewidths=3.0,
                                  zorder=100)
                        axs1.annotate(s=gps_label, xy=(p_lon + 0.015, p_lat), horizontalalignment='left')
                    for p_lon_ts, p_lat_ts, ts_label in zip(points_lon_use, points_lat_use,
                                                            labels_use):
                        m.scatter(p_lon_ts, p_lat_ts, facecolor='none', edgecolor=color_, s=s2, marker='o',
                                  linewidths=3.0,
                                  zorder=101)
                        axs1.annotate(s=ts_label, xy=(p_lon_ts + 0.015, p_lat_ts), horizontalalignment='left')
                else:
                    plt.rcParams.update({'text.color': color_text,
                                         'font.size': 24})
                    for p_lon, p_lat, gps_label in zip(self.points_gps_lon, self.points_gps_lat, self.gps_labels):
                        m.scatter(p_lon, p_lat, facecolor='none', edgecolor=color_, s=s1, marker='*', linewidths=3.0,
                                  zorder=100)
                        axs1.annotate(s=gps_label, xy=(p_lon + 0.015, p_lat), horizontalalignment='left')
                    for p_lon_ts, p_lat_ts, ts_label in zip(points_lon_use, points_lat_use,
                                                            labels_use):
                        m.scatter(p_lon_ts, p_lat_ts, facecolor='none', edgecolor=color_, s=s2, marker='o',
                                  linewidths=3.0,
                                  zorder=101)
                        axs1.annotate(s=ts_label, xy=(p_lon_ts + 0.015, p_lat_ts), horizontalalignment='left')
            mpl.rcParams.update(mpl.rcParamsDefault)
        if self.composite == 'true':
            if not self.corners == None:
                patches = []
                for Frame in self.corners[0:2]:
                    patches.append(Polygon_mpl(Frame))
                axs1.add_collection(PatchCollection(patches, facecolor='none', edgecolor='r', linewidths=3.0))
                patches3 = []
                for Frame in self.corners[2:]:
                    patches3.append(Polygon_mpl(Frame))
                axs1.add_collection(PatchCollection(patches3, facecolor='none', edgecolor='r', linewidths=3.0))
        parallels = np.arange(bottom_use, top_use + 1, base_degree_lon)
        meridians = np.arange(left_use, right_use + 1, base_degree_lat)
        if self.ndvi_overlay == 'true':
            plt.rcParams.update({'text.color': "black",
                                 'font.size': 26})
            m1.drawparallels(parallels, labels=[True, False, True, False], linewidth=3.0, fontsize=16)
            m1.drawmeridians(meridians, labels=[True, False, True, False], linewidth=3.0, fontsize=16)
            m2.drawparallels(parallels, labels=[True, False, True, False], linewidth=3.0, fontsize=16)
            m2.drawmeridians(meridians, labels=[True, False, True, False], linewidth=3.0, fontsize=16)
            mpl.rcParams.update(mpl.rcParamsDefault)
        else:
            plt.rcParams.update({'text.color': "black",
                                 'font.size': 38})
            m.drawparallels(parallels, labels=[True, False, True, False], linewidth=3.0, fontsize=26)
            m.drawmeridians(meridians, labels=[True, False, True, False], linewidth=3.0, fontsize=26)
            mpl.rcParams.update(mpl.rcParamsDefault)
        if self.ndvi_overlay == 'true':
            plt.rcParams.update({'text.color': "black",
                                 'font.size': 18})
            divider1 = make_axes_locatable(axs1)
            cax1 = divider1.append_axes("right", size="5%", pad=0.05)
            cb_ticks_ndvi = list(np.arange(lower_bound_u_z_ndvi, upper_bound_u_z_ndvi + 0.2, 0.2))
            cb_ndvi = fig.colorbar(im1, ax=axs1, cax=cax1, extend='both', ticks=cb_ticks_ndvi)
            cb_ndvi.set_label('ndvi')

            divider2 = make_axes_locatable(axs2)
            cax2 = divider2.append_axes("right", size="5%", pad=0.05)
            cb_ticks = list(np.arange(lower_bound_u_z, upper_bound_u_z + 1, self.base_linear))
            cb = fig.colorbar(im2, ax=axs2, cax=cax2, extend='both', ticks=cb_ticks)
            cb.set_label('mm/yr')
            mpl.rcParams.update(mpl.rcParamsDefault)
        elif self.composite == 'true':
            pass
        elif self.base_map_make == 'true':
            pass
        else:
            plt.rcParams.update({'text.color': "black",
                                 'font.size': 24})
            divider1 = make_axes_locatable(axs1)
            cax1 = divider1.append_axes("right", size="5%", pad=0.05)
            if self.linear_or_coh == 'linear':
                cb_ticks = list(np.arange(lower_bound_u_z, upper_bound_u_z + 1, self.base_linear))
                cb = fig.colorbar(im1, ax=axs1, cax=cax1, extend='both', ticks=cb_ticks)
                cb.set_label('mm/yr')
            else:
                cb_ticks = list(np.arange(lower_bound_u_z, upper_bound_u_z + 0.1, self.base_coh))
                cb = fig.colorbar(im1, ax=axs1, cax=cax1, extend='both', ticks=cb_ticks)
                cb.set_label('coherence')
            mpl.rcParams.update(mpl.rcParamsDefault)

        # axs1.set_title(self.save_name_title)
        if not self.ndvi_overlay == 'true':
            box = axs1.get_position()
            axs1.set_position([box.x0, box.y0 + box.height * 0.1,
                               box.width, box.height * 0.9])

        axs1.patch.set_edgecolor('black')
        axs1.patch.set_linewidth('3')
        if self.geol_overlay == 'true':
            plt.rcParams.update({'text.color': "black",
                                 'font.size': 18})
            axs1.legend(handles=legend_com, loc='center right', bbox_to_anchor=(-0.08, 0.5), fancybox=True, shadow=True,
                        ncol=1, handlelength=5, handleheight=5, facecolor="w",
                        edgecolor='k')
            mpl.rcParams.update(mpl.rcParamsDefault)
        if self.path_save_name.exists():
            os.remove(self.path_save_name)
        plt.tight_layout()
        fig.savefig(self.path_save_name, dpi=600, facecolor='w', edgecolor='k',
                    orientation='portrait', bbox_inches='tight',
                    pad_inches=0.3)
        plt.close(fig)
        if self.make_examples == 'true':
            text_label = ['A', 'B', 'C', 'D', 'E', 'F']
            text_label_1 = ['A', 'B']
            fig1, axs1 = plt.subplots(2, 3, figsize=(18, 9.1), sharex=True, sharey=True)
            fig1.subplots_adjust(hspace=0, wspace=0)
            axs1 = axs1.ravel()
            fig2, axs2 = plt.subplots(1, 2, figsize=(14, 10), sharex=True, sharey=True)
            fig2.subplots_adjust(hspace=0, wspace=0)
            axs2 = axs2.ravel()
            example_ifg_wrapped_array = rasterio_basic_functions.tif_2_array(self.example_ifg_wrapped_path)
            example_ifg_wrapped_array[example_ifg_wrapped_array == 0.] = np.nan
            example_ifg_unwrapped_path_array = rasterio_basic_functions.tif_2_array(self.example_ifg_unwrapped_path)
            example_ifg_unwrapped_path_array[example_ifg_unwrapped_path_array == 0.] = np.nan
            example_coh_path_array = rasterio_basic_functions.tif_2_array(self.example_coh_path)
            example_non_gacos_path_array = rasterio_basic_functions.tif_2_array(self.example_non_gacos_path)
            example_non_gacos_path_array[example_non_gacos_path_array == 0.] = np.nan
            example_gacos_path_array = rasterio_basic_functions.tif_2_array(self.example_gacos_path)
            example_gacos_path_array[example_gacos_path_array == 0.] = np.nan
            coh_mask_make = np.where(example_coh_path_array < 0.15,
                                     0,
                                     1)
            example_coh_path_array[example_coh_path_array == 0.] = np.nan
            ifg_example_wrapped_masked = np.multiply(coh_mask_make, example_ifg_wrapped_array)
            ifg_example_wrapped_masked[ifg_example_wrapped_masked == 0.] = np.nan
            ifg_example_unwrapped_masked = np.multiply(coh_mask_make, example_ifg_unwrapped_path_array)
            ifg_example_unwrapped_masked[ifg_example_unwrapped_masked == 0.] = np.nan

            cmap_wrapped = mpl.cm.hsv
            cmap_unwrapped = mpl.cm.rainbow
            cmap_coh = mpl.cm.Blues.reversed()
            cmap_non_gacos_gacos = mpl.cm.rainbow
            cmap_masked_coh = mpl.cm.Greys
            axs1[0].imshow(example_coh_path_array, cmap=cmap_coh, interpolation=None, origin="upper")
            axs1[1].imshow(example_ifg_wrapped_array, cmap=cmap_wrapped, interpolation=None, origin="upper")
            axs1[2].imshow(example_ifg_unwrapped_path_array, cmap=cmap_unwrapped, interpolation=None, origin="upper")
            axs1[3].imshow(coh_mask_make, cmap=cmap_masked_coh, interpolation=None, origin="upper")
            axs1[4].imshow(ifg_example_wrapped_masked, cmap=cmap_wrapped, interpolation=None, origin="upper")
            axs1[5].imshow(ifg_example_unwrapped_masked, cmap=cmap_unwrapped, interpolation=None, origin="upper")

            max_gacos, min_gacos = np.nanmax(example_gacos_path_array), np.nanmin(example_gacos_path_array)
            max_non_gacos, min_non_gacos = (np.nanmax(example_non_gacos_path_array),
                                            np.nanmin(example_non_gacos_path_array))
            if max_gacos < max_non_gacos:
                max_ = max_non_gacos
            else:
                max_ = max_gacos
            if min_gacos < min_non_gacos:
                min_ = min_gacos
            else:
                min_ = min_non_gacos
            norm_ = mpl.colors.Normalize(vmin=min_, vmax=max_)
            axs2[0].imshow(example_non_gacos_path_array, cmap=cmap_non_gacos_gacos, interpolation=None,
                           vmin=min_, vmax=max_, norm=norm_,
                           origin="upper")
            axs2[1].imshow(example_gacos_path_array, cmap=cmap_non_gacos_gacos, interpolation=None,
                           vmin=min_, vmax=max_, norm=norm_,
                           origin="upper")
            for ax1, lab_ in zip(axs1, text_label):
                ax1.text(0.02, 0.93, lab_, fontsize=18, transform=ax1.transAxes)
                ax1.set_xticks([])
                ax1.set_yticks([])
            for ax2, lab_ in zip(axs2, text_label_1):
                ax2.text(0.02, 0.93, lab_, fontsize=18, transform=ax2.transAxes)
                ax2.set_xticks([])
                ax2.set_yticks([])
            fig1.savefig(self.example_ifg_masked_save_path, dpi=600, facecolor='w', edgecolor='k',
                         orientation='portrait', bbox_inches='tight',
                         pad_inches=0.3)
            fig2.savefig(self.example_non_gacos_gacos_save_path, dpi=600, facecolor='w', edgecolor='k',
                         orientation='portrait', bbox_inches='tight',
                         pad_inches=0.3)
            plt.close(fig1)
            plt.close(fig2)

    def make_histograms(self):
        base_array = rasterio_basic_functions.tif_2_array(self.paths_to_use)
        base_array[base_array == 0.] = np.nan
        if self.linear_or_coh == 'linear':
            base_array = base_array * 1000
        h_a = base_array.flatten()[~np.isnan(base_array.flatten())]
        std_ = np.std(h_a)
        mean_ = np.mean(h_a)
        mean_std_1_pos = mean_ + 1 * std_
        mean_std_2_pos = mean_ + 2 * std_
        mean_std_3_pos = mean_ + 3 * std_
        mean_std_1_neg = mean_ - 1 * std_
        mean_std_2_neg = mean_ - 2 * std_
        mean_std_3_neg = mean_ - 3 * std_
        two_std_pos = mean_ + 3 * std_
        two_std_neg = mean_ - 3 * std_
        h, b = np.histogram(h_a, bins=150, density=True)
        fig, axs = plt.subplots(1, 1, figsize=(15, 10), facecolor='w', edgecolor='k')
        sns.histplot(h_a,
                     stat='density',
                     bins=350,
                     kde=True,
                     color='b',
                     ax=axs)
        # axs.hist(h_a, bins=b)
        if self.linear_or_coh == 'linear':
            axs.set_xlabel('mm/yr')
        else:
            axs.set_xlabel('coherence')
        if self.path_save_name_hist.exists():
            os.remove(self.path_save_name_hist)
        axs.set_xlim([two_std_neg, two_std_pos])
        # axs.set_title(self.save_name_title)
        fig.savefig(self.path_save_name_hist, dpi=200, facecolor='w', edgecolor='k',
                    orientation='portrait', bbox_inches='tight', pad_inches=0.3)
        if self.path_save_name_hist_m_std.exists():
            os.remove(self.path_save_name_hist_m_std)
        axs.axvline(x=mean_std_1_pos, color='k')
        axs.axvline(x=mean_std_1_neg, color='k')
        axs.axvline(x=mean_std_2_pos, color='k')
        axs.axvline(x=mean_std_2_neg, color='k')
        axs.axvline(x=mean_std_3_pos, color='k')
        axs.axvline(x=mean_std_3_neg, color='k')
        axs.axvline(x=mean_, color='r')
        fig.savefig(self.path_save_name_hist_m_std, dpi=200, facecolor='w', edgecolor='k',
                    orientation='portrait', bbox_inches='tight', pad_inches=0.3)

    def make_topo_overlay(self):
        deformation_array = rasterio_basic_functions.tif_2_array(self.combined_linear_path_dem_crop)
        topo_array = rasterio_basic_functions.tif_2_array(self.path_2_dem)
        deformation_array_nz = deformation_array[deformation_array != 0].flatten() * 1000
        topo_array_nz = topo_array[topo_array != 0].flatten()
        correlation_coeff = np.round(np.corrcoef(topo_array_nz, deformation_array_nz)[0, 1], 2)
        fig, axs1 = plt.subplots(figsize=(12, 10), facecolor='w', edgecolor='k')
        base = 3
        n_rows = 3
        col_size = 20 / (base - n_rows + 1)
        fig3, axs10 = plt.subplots(figsize=(24, col_size), nrows=n_rows, ncols=3, facecolor='w',
                                   edgecolor='k')
        axs10 = axs10.ravel()
        width_ = 3
        save_names = [self.topo_correlation_save_name_1,
                      self.topo_correlation_save_name_2,
                      self.topo_correlation_save_name_3,
                      self.topo_correlation_save_name_4,
                      self.topo_correlation_save_name_5,
                      self.topo_correlation_save_name_6,
                      self.topo_correlation_save_name_7,
                      self.topo_correlation_save_name_8,
                      self.topo_correlation_save_name_9]
        plt.rcParams.update({'text.color': "black",
                             'font.size': 18})
        kk = 0
        for lons, lats, names_, start_name, end_name in zip(self.cross_section_lon_list,
                                                            self.cross_section_lat_list,
                                                            save_names,
                                                            self.cross_section_labels_start,
                                                            self.cross_section_labels_end):
            fig2, axs2 = plt.subplots(figsize=(16, 16), facecolor='w', edgecolor='k')
            axs4 = axs2.twinx()
            deformation_array_c = rasterio_basic_functions.tif_2_array(self.combined_linear_path_dem_crop)
            topo_array_c = rasterio_basic_functions.tif_2_array(self.path_2_dem)
            defo_norm_cross, topo_norm_cross, index_, start_coord, end_coord = (
                make_cross_section_line(self.combined_linear_path_dem_crop, deformation_array_c,
                                        topo_array_c, lons, lats, width_))
            base = 20
            topo_mean = np.nanmean(topo_norm_cross)
            defo_mean = np.nanmean(defo_norm_cross)
            defo_mean_add = 3.2 * abs(topo_mean - defo_mean)
            y_lim_topo_top = 2 * general_functions.nearest_bound_coord(np.nanmax(topo_norm_cross), base, 'up')
            y_lim_defo_top = (general_functions.nearest_bound_coord(
                np.nanmax(defo_norm_cross) + defo_mean_add, base, 'up'))
            defo_norm_cross[defo_norm_cross == 0] = np.nan
            topo_norm_cross[topo_norm_cross == 0] = np.nan
            axs2.scatter(index_, defo_norm_cross + defo_mean_add, color='r', s=10, label='Deformation')
            axs2.scatter(index_, topo_norm_cross, color='k', s=10, label='DEM')
            axs3 = axs2.twiny()
            x_ = [index_[0], index_[-1]]
            axs3.set_xbound(axs2.get_xbound())
            axs2.set_xticks(x_)
            axs2.set_xticklabels([str(start_coord), str(end_coord)])
            axs3.set_xticklabels([start_name, end_name])
            axs3.set_xticks(x_)
            axs2.set_yticks(np.arange(0, y_lim_defo_top, base))
            axs2.set_yticklabels([str(x) for x in np.arange(0, y_lim_defo_top, base)], fontsize=18)
            axs2.tick_params(labelsize=18)
            axs4.set_yticks(np.arange(0, y_lim_topo_top, base))
            axs2.set_ylabel('Deformation (mm/yr)')
            axs4.set_ylabel('DEM (m)')
            axs2.legend()
            fig2.savefig(names_, dpi=100, facecolor='w', edgecolor='k',
                         orientation='portrait', bbox_inches='tight', pad_inches=0.3)
            plt.close(fig2)
            axs10[kk].scatter(index_, defo_norm_cross + defo_mean_add, color='r', s=10, label='Deformation')
            axs10[kk].scatter(index_, topo_norm_cross, color='k', s=10, label='DEM')
            axs14 = axs10[kk].twinx()
            axs13 = axs10[kk].twiny()
            x_ = [index_[0], index_[-1]]
            axs13.set_xbound(axs10[kk].get_xbound())
            axs10[kk].set_xticks(x_)
            axs10[kk].set_xticklabels([str(start_coord), str(end_coord)], fontsize=18)
            axs13.set_xticklabels([start_name, end_name])
            axs13.set_xticks(x_)
            axs10[kk].set_yticks(np.arange(0, y_lim_defo_top, base))
            axs14.set_yticks(np.arange(0, y_lim_topo_top, base))
            axs10[kk].set_ylabel('Deformation (mm/yr)', fontsize=18)
            axs14.set_ylabel('DEM (m)')
            axs10[kk].legend()
            kk = kk + 1
        plt.tight_layout()
        fig3.savefig(self.topo_correlation_save_name_all, dpi=100, facecolor='w', edgecolor='k',
                     orientation='portrait', bbox_inches='tight', pad_inches=0.3)
        axs1.scatter(topo_array_nz, deformation_array_nz, color='k', s=10)
        axs1.set_ylabel('deformation (mm/yr)', fontsize=16)
        axs1.set_xlabel('dem height (m)', fontsize=16)
        axs1.text(0.01, 0.97, ('pearson correlation coefficient=' + str(correlation_coeff)),
                  fontsize=12, transform=axs1.transAxes)
        fig.savefig(self.topo_correlation_save_name, dpi=100, facecolor='w', edgecolor='k',
                    orientation='portrait', bbox_inches='tight', pad_inches=0.3)
        mpl.rcParams.update(mpl.rcParamsDefault)
        plt.close(fig)
        plt.close(fig3)

    def return_relevant_statistics(self):
        print('finding_stats')
        base_array = rasterio_basic_functions.tif_2_array(self.paths_to_use)
        base_array[base_array == 0.] = np.nan
        base_array = base_array * 1000
        half_size_ = 5
        names_of_average = self.time_series_labels_all[0:3]
        mean_tif = np.nanmean(base_array)
        std_tif = np.nanstd(base_array)
        mean_std_1 = mean_tif + std_tif
        mean_std_2 = mean_tif + 2 * std_tif
        mean_std_3 = mean_tif + 3 * std_tif
        mean_std_4 = mean_tif - std_tif
        mean_std_5 = mean_tif - 2 * std_tif
        mean_std_6 = mean_tif - 3 * std_tif
        param_keys_all = ['mean_linear', 'std_linear', 'mean_std_1', 'mean_std_2', 'mean_std_3', 'mean_std_4',
                          'mean_std_5', 'mean_std_6']
        param_vals_all = [str(mean_tif), str(std_tif), str(mean_std_1), str(mean_std_2), str(mean_std_3),
                          str(mean_std_4), str(mean_std_5), str(mean_std_6)]
        gps_average_rates = []
        deformation_names = []
        deformation_rates = []
        gps_names = []
        lon_average = []
        lat_average = []
        for gps_base_name, gps_base_lon, gps_base_lat, gps_base_rate in zip(names_of_average,
                                                                            self.points_gps_lon[1:],
                                                                            self.points_gps_lat[1:],
                                                                            self.gps_rates_3):
            gps_names.append(gps_base_name)
            lon_average.append(gps_base_lon)
            lat_average.append(gps_base_lat)
            gps_average_rates.append(gps_base_rate)
        for gps_base_name, gps_base_lon, gps_base_lat, gps_base_rate in zip(self.names_named_all,
                                                                            self.points_lon_named_all,
                                                                            self.points_lat_named_all,
                                                                            self.rates_named_all):
            gps_names.append(gps_base_name)
            lon_average.append(gps_base_lon)
            lat_average.append(gps_base_lat)
            gps_average_rates.append(gps_base_rate)
        for gps_base_name, gps_base_lon, gps_base_lat, gps_base_rate in zip(self.names_campaign,
                                                                            self.points_lon_campaign,
                                                                            self.points_lat_campaign,
                                                                            self.rates_campaign):
            gps_names.append(gps_base_name)
            lon_average.append(gps_base_lon)
            lat_average.append(gps_base_lat)
            gps_average_rates.append(gps_base_rate)
        jj = 0
        for lon_, lat_ in zip(lon_average, lat_average):
            average_ = find_average_value_at_point(base_array, lon_, lat_, half_size_, self.paths_to_use)
            deformation_rates.append(str(average_))
            def_name = 'deformation_point_' + str(jj)
            deformation_names.append(def_name)
            jj = jj + 1
        list_of_all = [deformation_names, deformation_rates,
                       gps_names, gps_average_rates,
                       lon_average, lat_average]
        gps_comparison_def_names = ['deformation point name', 'deformation rates (mm/yr)',
                                    'gps name', 'gps deformation rates (mm/yr)',
                                    'longitude', 'latitude']
        gps_comparison_def_names_dict = {}
        for gps_names_def, list_ in zip(gps_comparison_def_names, list_of_all):
            gps_comparison_def_names_dict[gps_names_def] = list_
        df_comparison = pd.DataFrame(gps_comparison_def_names_dict)
        df_comparison = df_comparison[df_comparison['deformation rates (mm/yr)'] != 0]
        df_comparison = df_comparison[df_comparison['deformation rates (mm/yr)'] != 'nan']
        df_comparison = df_comparison[df_comparison['deformation rates (mm/yr)'] != np.nan]
        general_functions.save_df_header(df_comparison, self.gps_def_comparison)
        param_dict = {}
        for param_name, param_val in zip(param_keys_all, param_vals_all):
            param_dict[param_name] = param_val
        general_functions.write_json_from_dict(param_dict, self.statistics_path)

    def make_time_series__(self):
        fig33, axs33 = (plt.subplots(figsize=(14, 10), nrows=1,
                                     ncols=3, facecolor='w', edgecolor='k'))
        axs33 = axs33.ravel()
        points_of_interest = ['P14', 'P19', 'P20', 'P21', 'P22', 'P24']
        fig44, axs44 = (plt.subplots(figsize=(18, 14), nrows=2,
                                     ncols=3, facecolor='w', edgecolor='k'))
        axs44 = axs44.ravel()
        ll = 0
        hh = 0
        for lon_, lat_, t_s_label in zip(self.points_ts_lon_all, self.points_ts_lat_all, self.time_series_labels_all):
            check_list = []
            for base_kwargs_s in self.base_kwargs_paths:
                base_kwargs_test_dict = general_functions.open_json_file(base_kwargs_s)
                base_kwargs_test_dict_transform = base_kwargs_test_dict['transform']
                base_kwargs_test_dict_width = base_kwargs_test_dict['width']
                base_kwargs_test_dict_height = base_kwargs_test_dict['height']
                affine_from_transform = rasterio_basic_functions.make_affine_from_list(base_kwargs_test_dict_transform)
                left, bot, right, top = (rasterio_basic_functions.
                                         find_bounds_from_affine_w_h(affine_from_transform,
                                                                     base_kwargs_test_dict_width,
                                                                     base_kwargs_test_dict_height))
                if (left < lon_ < right) and (bot < lat_ < top):
                    base_value = 1
                    check_list.append(base_value)
                else:
                    base_value = 0
                    check_list.append(base_value)
            check_list_np = np.array(check_list)
            index = np.where(check_list_np == 1)[0]
            files_found = []
            for index_ in index:
                base_file = self.base_kwargs_paths[index_]
                files_found.append(base_file.parent.parent)
            self.to_mm = 1000
            self.lon_list = [lon_]
            self.lat_list = [lat_]
            self.half_size = 5
            if len(files_found) == 1:
                time_series_save_name = t_s_label
                self.base_path_time_series = files_found[0]
                self.kwargs_transform_path = self.base_path_time_series / 'tif_meta_data/kwargs_tif.json'
                self.binary_spatial_array_path = self.base_path_time_series / 'Binary_Spatial_Arrays'
                self.tsout_path_all = self.base_path_time_series / 'Data_used_regression/MSBAS_TSOUT_all.json'
                self.time_series_path = self.combined_tifs_path / 'Time_Series_Plots'
                if not self.time_series_path.exists():
                    os.mkdir(self.time_series_path)
                gps_path_ = [Path(x) for x in self.gps_time_series_paths if time_series_save_name in str(x)]
                if gps_path_:
                    gps_df = parse_gps_txt_file(gps_path_)
                    gps_save_name = (self.time_series_path /
                                     ('Time_Series_' + time_series_save_name + '_gps_overlay.png'))
                self.time_series_save_name = self.time_series_path / ('Time_Series_' + time_series_save_name + '.png')
                time_series_save_name_json = (self.time_series_path / ('Time_Series_' +
                                                                       time_series_save_name + '.json'))
                kwargs_open = general_functions.open_json_file(self.kwargs_transform_path)
                self.transform_list = kwargs_open['transform']
                tsout_open_all = general_functions.open_json_file(self.tsout_path_all)
                self.times_all = np.array([float(x) for x in tsout_open_all['YYYY.YYY']])
                self.times_all_shifted = self.times_all - np.min(self.times_all)
                msbas_used_np_all = tsout_open_all['msbas_file_names']
                test_nps_all = [Path(str(self.binary_spatial_array_path) + '/' + x) for x in msbas_used_np_all]
                test_nps_all.sort()
                self.binary_spatial_arrays_all = test_nps_all
                lon_cols, lat_rows = general_functions.lon_lat_2_col_row_using_transform(self.transform_list,
                                                                                         self.lon_list,
                                                                                         self.lat_list)
                def_all, coeff_all = mean_sliced_list(self.binary_spatial_arrays_all,
                                                      self.half_size,
                                                      lon_cols,
                                                      lat_rows,
                                                      self.times_all_shifted,
                                                      self.to_mm)
                base_a = np.array([])
                for def_r in def_all:
                    np_as = np.array(def_r)
                    base_a = np.append(base_a, np_as)
                base_t = 20
                ylim_min = general_functions.nearest_bound_coord(np.min(base_a), base_t, 'down')
                ylim_max = general_functions.nearest_bound_coord(np.max(base_a), base_t, 'up')
                ylim_min = ylim_min
                ylim_max = ylim_max
                base = 3
                if len(self.lon_list) == 1:
                    n_rows = 1
                    n_cols = 1
                    fig, axs1 = plt.subplots(figsize=(12, 10), nrows=n_rows, ncols=n_cols, facecolor='w', edgecolor='k')
                    if gps_path_:
                        fig2, axs2 = (plt.subplots(figsize=(12, 10), nrows=n_rows,
                                                   ncols=n_cols, facecolor='w', edgecolor='k'))
                else:
                    num_plots = general_functions.nearest_bound_coord(len(self.lon_list), base, 'up')
                    n_rows = int(num_plots / base)
                    col_size = 20 / (base - n_rows + 1)
                    fig, axs1 = plt.subplots(figsize=(24, col_size), nrows=n_rows, ncols=3, facecolor='w',
                                             edgecolor='k')
                    axs1 = axs1.ravel()
                for i in range(len(self.lon_list)):
                    if len(self.lon_list) == 1:
                        axs1.set_ylim([ylim_min, ylim_max])
                        axs1.set_yticks(np.arange(ylim_min, ylim_max, base_t))
                        axs1.scatter(self.times_all, def_all, color='k')
                        poly = np.poly1d(coeff_all[0])
                        rate = round(coeff_all[0][0], 2)
                        new_y1 = poly(self.times_all_shifted)
                        axs1.plot(self.times_all, new_y1, color='b')
                        axs1.set_ylabel('mm', fontsize=16)
                        x_ticks = [2017, 2017.5, 2018, 2018.5, 2019, 2019.5, 2020, 2020.5, 2021, 2021.5]
                        axs1.set_xlim([2017, 2021.5])
                        axs1.set_xticks(x_ticks)
                        axs1.text(0.4, 0.97, (t_s_label + ' linear_rate=' + str(rate) + ' mm/yr'), fontsize=12,
                                  transform=axs1.transAxes)
                        axs1.text(0.4, 0.93, ('lon=' + str(self.lon_list[0]) + ' ' + 'lat=' + str(self.lat_list[0])),
                                  fontsize=12,
                                  transform=axs1.transAxes)
                        if t_s_label in points_of_interest:
                            T = 1
                            param, covariance = curve_fit(sinusoid,
                                                          self.times_all,
                                                          np.array(def_all).flatten(),
                                                          p0=get_p0(self.times_all, np.array(def_all).flatten(), T)
                                                          )
                            sinusoid_values = sinusoid(self.times_all, *param)
                            axs44[hh].set_ylim([ylim_min, ylim_max])
                            axs44[hh].set_yticks(np.arange(ylim_min, ylim_max, base_t))
                            axs44[hh].scatter(self.times_all, def_all, color='k', label='Deformation')
                            axs44[hh].plot(self.times_all, sinusoid_values, color='b', linewidth=1)
                            poly = np.poly1d(coeff_all[0])
                            rate = round(coeff_all[0][0], 2)
                            new_y1 = poly(self.times_all_shifted)
                            axs44[hh].plot(self.times_all, new_y1, color='b')
                            axs44[hh].scatter(self.times_all,
                                              (np.array(def_all).flatten() - sinusoid_values), color='m', s=10,
                                              label='Def - Sinusoid')
                            axs44[hh].set_ylabel('mm', fontsize=16)
                            x_ticks = [2017, 2017.5, 2018, 2018.5, 2019, 2019.5, 2020, 2020.5, 2021, 2021.5]
                            axs44[hh].set_xlim([2017, 2021.5])
                            axs44[hh].set_xticks(x_ticks)
                            axs44[hh].text(0.4, 0.97, (t_s_label + ' linear_rate=' + str(rate) + ' mm/yr'), fontsize=12,
                                           transform=axs44[hh].transAxes)
                            axs44[hh].text(0.4, 0.93,
                                           ('lon=' + str(self.lon_list[0]) + ' ' + 'lat=' + str(self.lat_list[0])),
                                           fontsize=12,
                                           transform=axs44[hh].transAxes)
                            axs44[hh].legend()
                            hh = hh + 1

                        time_series_dict = {'x_Time_Years': self.times_all,
                                            'y_deformation_mm': def_all,
                                            'y_lim_start_end': [ylim_min, ylim_max],
                                            'y_regression': new_y1,
                                            'x_ticks': x_ticks,
                                            'lon_lat': [self.lon_list[0], self.lat_list[0]]
                                            }
                        if gps_path_:
                            T = 1
                            axs2.set_ylim([ylim_min, ylim_max])
                            axs2.set_yticks(np.arange(ylim_min, ylim_max, base_t))
                            axs2.scatter(self.times_all, def_all, color='k')
                            axs2.scatter(gps_df['YYYY_YYY'], gps_df['dU'] * 1000, color='r', s=10)
                            axs2.plot(self.times_all, new_y1, color='b')

                            param, covariance = curve_fit(sinusoid,
                                                          self.times_all,
                                                          np.array(def_all).flatten(),
                                                          p0=get_p0(self.times_all, np.array(def_all).flatten(), T)
                                                          )
                            sinusoid_values = sinusoid(self.times_all, *param)
                            axs33[ll].set_ylim([ylim_min, ylim_max])
                            axs33[ll].set_yticks(np.arange(ylim_min, ylim_max, base_t))
                            axs33[ll].scatter(self.times_all, def_all, color='k', label='Deformation')
                            axs33[ll].scatter(gps_df['YYYY_YYY'], gps_df['dU'] * 1000, color='r', s=10,
                                              label='GPS')
                            axs33[ll].plot(self.times_all, new_y1, color='b')
                            axs33[ll].plot(self.times_all, sinusoid_values, color='b', linewidth=1)
                            axs33[ll].scatter(self.times_all,
                                              (np.array(def_all).flatten() - sinusoid_values), color='m', s=10,
                                              label='Def - Sinusoid')
                            axs33[ll].set_ylabel('mm', fontsize=16)
                            axs33[ll].set_xlim([2017, 2021.5])
                            axs33[ll].set_xticks([2017, 2017.5, 2018, 2018.5, 2019, 2019.5, 2020, 2020.5, 2021, 2021.5])
                            axs33[ll].text(0.4, 0.97, (t_s_label + ' linear_rate=' + str(rate) + ' mm/yr'), fontsize=12,
                                           transform=axs33[ll].transAxes)
                            axs33[ll].text(0.4, 0.93,
                                           ('lon=' + str(self.lon_list[0]) + ' ' + 'lat=' + str(self.lat_list[0])),
                                           fontsize=12,
                                           transform=axs33[ll].transAxes)
                            axs33[ll].legend()
                            axs2.set_ylabel('mm', fontsize=16)
                            axs2.set_xlim([2017, 2021.5])
                            axs2.set_xticks([2017, 2017.5, 2018, 2018.5, 2019, 2019.5, 2020, 2020.5, 2021, 2021.5])
                            axs2.text(0.4, 0.97, (t_s_label + ' linear_rate=' + str(rate) + ' mm/yr'), fontsize=12,
                                      transform=axs2.transAxes)
                            axs2.text(0.4, 0.93,
                                      ('lon=' + str(self.lon_list[0]) + ' ' + 'lat=' + str(self.lat_list[0])),
                                      fontsize=12,
                                      transform=axs2.transAxes)
                            ll = ll + 1
                        general_functions.write_json_from_dict_ts(time_series_dict, time_series_save_name_json)
                    else:
                        axs1[i].set_ylim([ylim_min, ylim_max])
                        axs1[i].scatter(self.times_all, def_all[i], color='k')
                        poly = np.poly1d(coeff_all[i])
                        rate = round(coeff_all[i][0], 2)
                        new_y1 = poly(self.times_all_shifted)
                        axs1[i].plot(self.times_all, new_y1, color='b')
                        axs1[i].set_ylabel('mm', fontsize=16)
                        axs1[i].set_xticks([2017, 2017.5, 2018, 2018.5, 2019, 2019.5, 2020, 2020.5, 2021, 2021.5])
                        axs1[i].text(0.4, 0.97, ('linear_rate=' + str(rate) + ' mm/yr'), fontsize=12,
                                     transform=axs1[i].transAxes)
                        axs1[i].text(0.4, 0.93,
                                     ('lon=' + str(self.lon_list[i]) + ' ' + 'lat=' + str(self.lat_list[i])),
                                     fontsize=12,
                                     transform=axs1[i].transAxes)
                if len(self.lon_list) == 1:
                    pass
                else:
                    num_plots = general_functions.nearest_bound_coord(len(self.lon_list), base, 'up')
                    num_plots_unused = num_plots - len(self.lon_list)
                    for j in range(len(self.lon_list), len(self.lon_list) + num_plots_unused):
                        fig.delaxes(axs1[j])
                fig.tight_layout()
                if self.time_series_save_name.exists():
                    os.remove(self.time_series_save_name)
                print('output time series png to:')
                print(self.time_series_save_name)
                fig.savefig(self.time_series_save_name, dpi=100, facecolor='w', edgecolor='k',
                            orientation='portrait', bbox_inches='tight', pad_inches=0.3)
                if gps_path_:
                    fig2.savefig(gps_save_name, dpi=100, facecolor='w', edgecolor='k',
                                 orientation='portrait', bbox_inches='tight', pad_inches=0.3)
                plt.close(fig)
                plt.close(fig2)
            else:
                base_1 = files_found[0]
                base_2 = files_found[1]
                time_series_save_name = t_s_label
                self.base_path_time_series_1 = base_1
                self.base_path_time_series_2 = base_2
                self.time_series_path = self.combined_tifs_path / 'Time_Series_Plots'
                if not self.time_series_path.exists():
                    os.mkdir(self.time_series_path)
                self.time_series_save_name = (self.time_series_path /
                                              ('Time_Series_' + time_series_save_name + '.png'))
                time_series_save_name_json = (self.time_series_path / ('Time_Series_' +
                                                                       time_series_save_name + '.json'))
                self.kwargs_transform_path_1 = self.base_path_time_series_1 / 'tif_meta_data/kwargs_tif.json'
                self.binary_spatial_array_path_1 = self.base_path_time_series_1 / 'Binary_Spatial_Arrays'
                self.tsout_path_all_1 = self.base_path_time_series_1 / 'Data_used_regression/MSBAS_TSOUT_all.json'
                kwargs_open_1 = general_functions.open_json_file(self.kwargs_transform_path_1)
                self.transform_list_1 = kwargs_open_1['transform']
                tsout_open_all_1 = general_functions.open_json_file(self.tsout_path_all_1)
                self.times_all_1 = np.array([float(x) for x in tsout_open_all_1['YYYY.YYY']])
                self.times_all_shifted_1 = self.times_all_1 - np.min(self.times_all_1)
                msbas_used_np_all_1 = tsout_open_all_1['msbas_file_names']
                test_nps_all_1 = [Path(str(self.binary_spatial_array_path_1) + '/' + x) for x in msbas_used_np_all_1]
                test_nps_all_1.sort()
                self.binary_spatial_arrays_all_1 = test_nps_all_1
                gps_path_ = [Path(x) for x in self.gps_time_series_paths if time_series_save_name in str(x)]
                if gps_path_:
                    gps_df = parse_gps_txt_file(gps_path_)
                    gps_save_name = (self.time_series_path /
                                     ('Time_Series_' + time_series_save_name + '_gps_overlay.png'))
                self.kwargs_transform_path_2 = self.base_path_time_series_2 / 'tif_meta_data/kwargs_tif.json'
                self.binary_spatial_array_path_2 = self.base_path_time_series_2 / 'Binary_Spatial_Arrays'
                self.tsout_path_all_2 = self.base_path_time_series_2 / 'Data_used_regression/MSBAS_TSOUT_all.json'
                kwargs_open_2 = general_functions.open_json_file(self.kwargs_transform_path_2)
                self.transform_list_2 = kwargs_open_2['transform']
                tsout_open_all_2 = general_functions.open_json_file(self.tsout_path_all_2)
                self.times_all_2 = np.array([float(x) for x in tsout_open_all_2['YYYY.YYY']])
                self.times_all_shifted_2 = self.times_all_2 - np.min(self.times_all_2)
                msbas_used_np_all_2 = tsout_open_all_2['msbas_file_names']
                test_nps_all_2 = [Path(str(self.binary_spatial_array_path_2) + '/' + x) for x in msbas_used_np_all_2]
                test_nps_all_2.sort()
                self.binary_spatial_arrays_all_2 = test_nps_all_2

                self.def_all = None
                self.coeff_all = None
                self.ylim_min = None
                self.ylim_max = None
                lon_cols_1, lat_rows_1 = general_functions.lon_lat_2_col_row_using_transform(self.transform_list_1,
                                                                                             self.lon_list,
                                                                                             self.lat_list)
                lon_cols_2, lat_rows_2 = general_functions.lon_lat_2_col_row_using_transform(self.transform_list_2,
                                                                                             self.lon_list,
                                                                                             self.lat_list)
                def_all_first = mean_sliced_list_no_poly(self.binary_spatial_arrays_all_1,
                                                         self.half_size,
                                                         lon_cols_1,
                                                         lat_rows_1,
                                                         self.to_mm)
                def_all_second = mean_sliced_list_no_poly(self.binary_spatial_arrays_all_2,
                                                          self.half_size,
                                                          lon_cols_2,
                                                          lat_rows_2,
                                                          self.to_mm)
                is_all_zero_1 = np.all((np.array(def_all_first) == 0))
                is_all_zero_2 = np.all((np.array(def_all_second) == 0))
                if not is_all_zero_1 and not is_all_zero_2:
                    time_1 = self.times_all_1[0]
                    time_2 = self.times_all_1[1]
                    if time_1 < time_2:
                        time_interp = self.times_all_1
                        def_interp = def_all_first
                        def_n_interp = def_all_second
                        time_2_use_interp = self.times_all_2
                        time_2_use_shifted = self.times_all_shifted_2
                    elif time_1 == time_2:
                        time_interp = self.times_all_2
                        def_interp = def_all_second
                        def_n_interp = def_all_first
                        time_2_use_interp = self.times_all_1
                        time_2_use_shifted = self.times_all_shifted_1
                    else:
                        time_interp = self.times_all_2
                        def_interp = def_all_second
                        def_n_interp = def_all_first
                        time_2_use_interp = self.times_all_1
                        time_2_use_shifted = self.times_all_shifted_1
                    f_2 = interp1d(time_interp, def_interp, kind='nearest', fill_value="extrapolate")
                    def_2_interp = f_2(time_2_use_interp)
                    def_averaged_1_2 = (def_2_interp + def_n_interp) / 2
                    if len(def_averaged_1_2) == 1:
                        def_averaged_1_2 = def_averaged_1_2[0]
                        if len(def_averaged_1_2) == 1:
                            def_averaged_1_2 = def_averaged_1_2[0]
                    coeff_averaged = np.polyfit(time_2_use_shifted, def_averaged_1_2, 1)
                    times_use_shifted = time_2_use_shifted
                    times_use = time_2_use_interp
                elif is_all_zero_1 and not is_all_zero_2:
                    def_averaged_1_2 = def_all_second
                    if len(def_averaged_1_2) == 1:
                        def_averaged_1_2 = def_averaged_1_2[0]
                        if len(def_averaged_1_2) == 1:
                            def_averaged_1_2 = def_averaged_1_2[0]
                    coeff_averaged = np.polyfit(self.times_all_shifted_2, def_averaged_1_2, 1)
                    times_use_shifted = self.times_all_shifted_2
                    times_use = self.times_all_2
                elif not is_all_zero_1 and is_all_zero_2:
                    def_averaged_1_2 = def_all_first
                    if len(def_averaged_1_2) == 1:
                        def_averaged_1_2 = def_averaged_1_2[0]
                        if len(def_averaged_1_2) == 1:
                            def_averaged_1_2 = def_averaged_1_2[0]
                    coeff_averaged = np.polyfit(self.times_all_shifted_1, def_averaged_1_2, 1)
                    times_use_shifted = self.times_all_shifted_1
                    times_use = self.times_all_1
                else:
                    print('both are zero arrays, reset point', str(self.lon_list[0]) + str(self.lat_list[0]))
                self.def_all = def_averaged_1_2
                self.coeff_all = [coeff_averaged]
                self.times_all_use = times_use
                self.times_all_shifted_use = times_use_shifted
                base_a = np.array([])
                for def_r in self.def_all:
                    np_as = np.array(def_r)
                    base_a = np.append(base_a, np_as)
                base_t = 20
                ylim_min = general_functions.nearest_bound_coord(np.min(base_a), base_t, 'down')
                ylim_max = general_functions.nearest_bound_coord(np.max(base_a), base_t, 'up')
                self.ylim_min = ylim_min
                self.ylim_max = ylim_max
                base = 3
                if len(self.lon_list) == 1:
                    n_rows = 1
                    n_cols = 1
                    fig, axs1 = plt.subplots(figsize=(12, 10), nrows=n_rows, ncols=n_cols, facecolor='w',
                                             edgecolor='k')
                    if gps_path_:
                        fig2, axs2 = (plt.subplots(figsize=(12, 10), nrows=n_rows,
                                                   ncols=n_cols, facecolor='w', edgecolor='k'))
                else:
                    num_plots = general_functions.nearest_bound_coord(len(self.lon_list), base, 'up')
                    n_rows = int(num_plots / base)
                    col_size = 20 / (base - n_rows + 1)
                    fig, axs1 = plt.subplots(figsize=(24, col_size), nrows=n_rows, ncols=3, facecolor='w',
                                             edgecolor='k')
                    axs1 = axs1.ravel()
                for i in range(len(self.lon_list)):
                    if len(self.lon_list) == 1:
                        axs1.set_ylim([self.ylim_min, self.ylim_max])
                        axs1.set_yticks(np.arange(self.ylim_min, self.ylim_max, base_t))
                        axs1.scatter(self.times_all_use, self.def_all, color='k')
                        poly = np.poly1d(self.coeff_all[0])
                        rate = round(self.coeff_all[0][0], 2)
                        new_y1 = poly(self.times_all_shifted_use)
                        axs1.plot(self.times_all_use, new_y1, color='b')
                        axs1.set_ylabel('mm', fontsize=16)
                        x_ticks = [2017, 2017.5, 2018, 2018.5, 2019, 2019.5, 2020, 2020.5, 2021, 2021.5]
                        axs1.set_xlim([2017, 2021.5])
                        axs1.set_xticks(x_ticks)
                        axs1.text(0.4, 0.97, (t_s_label + ' linear_rate=' + str(rate) + ' mm/yr'), fontsize=12,
                                  transform=axs1.transAxes)
                        axs1.text(0.4, 0.93,
                                  ('lon=' + str(self.lon_list[0]) + ' ' + 'lat=' + str(self.lat_list[0])),
                                  fontsize=12,
                                  transform=axs1.transAxes)

                        if t_s_label in points_of_interest:
                            T = 1
                            param, covariance = curve_fit(sinusoid,
                                                          self.times_all_use,
                                                          np.array(self.def_all).flatten(),
                                                          p0=get_p0(self.times_all_use,
                                                                    np.array(self.def_all).flatten(), T)
                                                          )
                            sinusoid_values = sinusoid(self.times_all_use, *param)
                            axs44[hh].set_ylim([self.ylim_min, self.ylim_max])
                            axs44[hh].set_yticks(np.arange(self.ylim_min, self.ylim_max, base_t))
                            axs44[hh].scatter(self.times_all_use, self.def_all, color='k', label='Deformation')
                            axs44[hh].plot(self.times_all_use, sinusoid_values, color='b', linewidth=1)
                            poly = np.poly1d(coeff_all[0])
                            rate = round(coeff_all[0][0], 2)
                            new_y1 = poly(self.times_all_shifted_use)
                            axs44[hh].plot(self.times_all_use, new_y1, color='b')
                            axs44[hh].scatter(self.times_all_use,
                                              (np.array(self.def_all).flatten() - sinusoid_values), color='m', s=10,
                                              label='Def - Sinusoid')
                            axs44[hh].set_ylabel('mm', fontsize=16)
                            x_ticks = [2017, 2017.5, 2018, 2018.5, 2019, 2019.5, 2020, 2020.5, 2021, 2021.5]
                            axs44[hh].set_xlim([2017, 2021.5])
                            axs44[hh].set_xticks(x_ticks)
                            axs44[hh].text(0.4, 0.97, (t_s_label + ' linear_rate=' + str(rate) + ' mm/yr'), fontsize=12,
                                           transform=axs44[hh].transAxes)
                            axs44[hh].text(0.4, 0.93,
                                           ('lon=' + str(self.lon_list[0]) + ' ' + 'lat=' + str(self.lat_list[0])),
                                           fontsize=12,
                                           transform=axs44[hh].transAxes)
                            axs44[hh].legend()
                            hh = hh + 1

                        time_series_dict = {'x_Time_Years': self.times_all_use,
                                            'y_deformation_mm': self.def_all,
                                            'y_lim_start_end': [self.ylim_min, self.ylim_max],
                                            'y_regression': new_y1,
                                            'x_ticks': x_ticks,
                                            'lon_lat': [self.lon_list[0], self.lat_list[0]]
                                            }
                        general_functions.write_json_from_dict_ts(time_series_dict, time_series_save_name_json)
                        if gps_path_:
                            T = 1
                            axs2.set_ylim([self.ylim_min, self.ylim_max])
                            axs2.scatter(self.times_all_use, self.def_all, color='k')
                            axs2.set_yticks(np.arange(self.ylim_min, self.ylim_max, base_t))
                            axs2.scatter(gps_df['YYYY_YYY'], gps_df['dU'] * 1000, color='r', s=10)
                            axs2.plot(self.times_all_use, new_y1, color='b')
                            axs2.set_ylabel('mm', fontsize=16)
                            axs2.set_xlim([2017, 2021.5])
                            axs2.set_xticks([2017, 2017.5, 2018, 2018.5, 2019, 2019.5, 2020, 2020.5, 2021, 2021.5])
                            axs2.text(0.4, 0.97, ('linear_rate=' + str(rate) + ' mm/yr'), fontsize=12,
                                      transform=axs2.transAxes)
                            axs2.text(0.4, 0.93,
                                      ('lon=' + str(self.lon_list[0]) + ' ' + 'lat=' + str(self.lat_list[0])),
                                      fontsize=12,
                                      transform=axs2.transAxes)
                            param, covariance = curve_fit(sinusoid,
                                                          self.times_all_use,
                                                          np.array(self.def_all).flatten(),
                                                          p0=get_p0(self.times_all, np.array(self.def_all).flatten(), T)
                                                          )
                            sinusoid_values = sinusoid(self.times_all_use, *param)
                            axs33[ll].set_ylim([self.ylim_min, self.ylim_max])
                            axs33[ll].set_yticks(np.arange(self.ylim_min, self.ylim_max, base_t))
                            axs33[ll].scatter(self.times_all_use, self.def_all, color='k', label='Deformation')
                            axs33[ll].scatter(gps_df['YYYY_YYY'], gps_df['dU'] * 1000, color='r', s=10, label='GPS')
                            axs33[ll].plot(self.times_all_use, new_y1, color='b')
                            axs33[ll].plot(self.times_all_use, sinusoid_values, color='b', linewidth=1)
                            axs33[ll].scatter(self.times_all_use,
                                              (np.array(self.def_all).flatten() - sinusoid_values), color='m', s=10,
                                              label='Def - Sinusoid')
                            axs33[ll].set_ylabel('mm', fontsize=16)
                            axs33[ll].set_xlim([2017, 2021.5])
                            axs33[ll].set_xticks([2017, 2017.5, 2018, 2018.5, 2019, 2019.5, 2020, 2020.5, 2021, 2021.5])
                            axs33[ll].text(0.4, 0.97, (t_s_label + ' linear_rate=' + str(rate) + ' mm/yr'), fontsize=12,
                                           transform=axs33[ll].transAxes)
                            axs33[ll].text(0.4, 0.93,
                                           ('lon=' + str(self.lon_list[0]) + ' ' + 'lat=' + str(self.lat_list[0])),
                                           fontsize=12,
                                           transform=axs33[ll].transAxes)
                            axs33[ll].legend()
                            ll = ll + 1
                    else:
                        axs1[i].set_ylim([self.ylim_min, self.ylim_max])
                        axs1[i].scatter(self.times_all_use, self.def_all[i], color='k')
                        poly = np.poly1d(self.coeff_all[i])
                        rate = round(self.coeff_all[i][0], 2)
                        new_y1 = poly(self.times_all_shifted_use)
                        axs1[i].plot(self.times_all_use, new_y1, color='b')
                        axs1[i].set_ylabel('mm', fontsize=16)
                        axs1[i].set_xticks([2017, 2017.5, 2018, 2018.5, 2019, 2019.5, 2020, 2020.5, 2021, 2021.5])
                        axs1[i].text(0.6, 0.97, (t_s_label + ' linear_rate=' + str(rate) + ' mm/yr'), fontsize=12,
                                     transform=axs1[i].transAxes)
                        axs1[i].text(0.6, 0.93,
                                     ('lon=' + str(self.lon_list[i]) + ' ' + 'lat=' + str(self.lat_list[i])),
                                     fontsize=12,
                                     transform=axs1[i].transAxes)
                if len(self.lon_list) == 1:
                    pass
                else:
                    num_plots = general_functions.nearest_bound_coord(len(self.lon_list), base, 'up')
                    num_plots_unused = num_plots - len(self.lon_list)
                    for j in range(len(self.lon_list), len(self.lon_list) + num_plots_unused):
                        fig.delaxes(axs1[j])
                fig.tight_layout()
                if self.time_series_save_name.exists():
                    os.remove(self.time_series_save_name)
                print('output time series png to:')
                print(self.time_series_save_name)
                fig.savefig(self.time_series_save_name, dpi=100, facecolor='w', edgecolor='k',
                            orientation='portrait', bbox_inches='tight', pad_inches=0.3)
                if gps_path_:
                    fig2.savefig(gps_save_name, dpi=100, facecolor='w', edgecolor='k',
                                 orientation='portrait', bbox_inches='tight', pad_inches=0.3)
                plt.close(fig)
        fig33.savefig(self.gps_sinusoid_name, dpi=100, facecolor='w', edgecolor='k',
                      orientation='portrait', bbox_inches='tight', pad_inches=0.3)
        fig44.savefig(self.points_of_interest_sinusoid_name, dpi=100, facecolor='w', edgecolor='k',
                      orientation='portrait', bbox_inches='tight', pad_inches=0.3)
        plt.close(fig33)
        plt.close(fig44)

        self.json_time_series_files = [Path(str(x)) for x in self.json_time_series_paths.glob("Time_Series_P*.json")]
        path_2_time_series_all_files_cleaned = [x for x in self.json_time_series_files if "PD32" not in str(x)]
        path_2_time_series_all_files_names = [int(x.stem.split('_')[2].split('P')[1:][0]) for x in
                                              path_2_time_series_all_files_cleaned]
        path_2_time_series_all_files_cleaned_sorted = [x for y, x in sorted(
            zip(path_2_time_series_all_files_names, path_2_time_series_all_files_cleaned), key=lambda pair: pair[0])]
        path_2_time_series_all_files_cleaned_sorted_names = [x.stem.split('_')[2] for x in
                                                             path_2_time_series_all_files_cleaned_sorted]
        num_of_time_series = len(path_2_time_series_all_files_cleaned_sorted)
        num_plots_per_axs = 9
        num_of_total_plots = int(num_of_time_series / num_plots_per_axs)
        plot_save_names = []
        plt.rcParams.update({'text.color': "black",
                             'font.size': 14})
        for i in range(0, num_of_total_plots):
            save_name_i = self.json_time_series_paths / ('Time_Series_Combined_' + str(i) + '.png')
            plot_save_names.append(save_name_i)
        for j in range(0, num_of_total_plots):
            save_name_ = plot_save_names[j]
            start_index = j * num_plots_per_axs
            end_index = (j + 1) * num_plots_per_axs
            plot_items_list = path_2_time_series_all_files_cleaned_sorted[start_index:end_index]
            plot_items_names_list = path_2_time_series_all_files_cleaned_sorted_names[start_index:end_index]
            base = 3
            n_rows = 3
            col_size = 20 / (base - n_rows + 1)
            fig, axs1 = plt.subplots(figsize=(24, col_size), nrows=n_rows, ncols=3, facecolor='w',
                                     edgecolor='k')
            axs1 = axs1.ravel()
            for k in range(0, len(plot_items_list)):
                item_2_plot = plot_items_list[k]
                name_2_plot = plot_items_names_list[k]
                dict_json = general_functions.open_json_file(item_2_plot)
                times_ = dict_json['x_Time_Years']
                def_ = dict_json['y_deformation_mm']
                y_lim = dict_json['y_lim_start_end']
                reg = dict_json['y_regression']
                x_ticks = dict_json['x_ticks']
                lon_lat__ = dict_json['lon_lat']
                rate_ = round((reg[-1] - reg[0]) / (np.max(times_ - np.min(times_))), 2)
                axs1[k].set_ylim(y_lim)
                axs1[k].scatter(times_, def_, color='k')
                axs1[k].plot(times_, reg, color='b')
                axs1[k].set_ylabel('mm', fontsize=16)
                axs1[k].set_xticks(x_ticks)
                axs1[k].text(0.4, 0.97, (name_2_plot + ' linear_rate=' + str(rate_) + ' mm/yr'), fontsize=12,
                             transform=axs1[k].transAxes)
                axs1[k].text(0.4, 0.93,
                             ('lon=' + str(lon_lat__[0]) + ' ' + 'lat=' + str(lon_lat__[1])),
                             fontsize=12,
                             transform=axs1[k].transAxes)
            fig.tight_layout()
            fig.savefig(save_name_, dpi=100, facecolor='w', edgecolor='k',
                        orientation='portrait', bbox_inches='tight', pad_inches=0.3)
            plt.close(fig)
        mpl.rcParams.update(mpl.rcParamsDefault)

    def run_all(self):
        self.define_master_resample_parameters()
        if not self.linear_rates_tif_combined_path_save.exists():
            self.combine_arrays()
        self.make_plots()
        if self.make_hist == 'true':
            self.make_histograms()
        if self.make_topo_correlation == 'true':
            self.make_topo_overlay()
        if self.return_statistics == 'true':
            self.return_relevant_statistics()
        if self.make_time_series == 'true':
            self.make_time_series__()


if __name__ == '__main__':
    sys_index_var_1 = sys.argv[1]
    sys_index_var_2 = sys.argv[2]
    sys_index_var_3 = sys.argv[3]
    sys_index_var_4 = sys.argv[4]
    sys_index_var_5 = sys.argv[5]
    sys_index_var_6 = sys.argv[6]
    sys_index_var_7 = sys.argv[7]
    sys_index_var_8 = sys.argv[8]
    sys_index_var_9 = sys.argv[9]
    sys_index_var_10 = sys.argv[10]
    sys_index_var_11 = sys.argv[11]
    sys_index_var_12 = sys.argv[12]
    sys_index_var_13 = sys.argv[13]
    sys_index_var_14 = sys.argv[14]
    combine_tifs = combine_all_tifs(sys_index_var_1, sys_index_var_2, sys_index_var_3, sys_index_var_4,
                                    sys_index_var_5, sys_index_var_6, sys_index_var_7, sys_index_var_8,
                                    sys_index_var_9, sys_index_var_10, sys_index_var_11, sys_index_var_12,
                                    sys_index_var_13, sys_index_var_14)
    combine_tifs.run_all()
