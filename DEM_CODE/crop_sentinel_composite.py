import sys, os

sys.path.append("..")
from helper_functions import rasterio_basic_functions
from helper_functions import general_functions as gen_fun
from pathlib import Path
import pathlib
import numpy as np
import math
from statistics import mean
from rasterio.windows import Window


def window_from_extent(tif_file, left_lon, right_lon, bot_lat, top_lat):
    ul_col, ul_row = rasterio_basic_functions.lon_lat_2_col_row(tif_file, left_lon, top_lat)
    lr_col, lr_row = rasterio_basic_functions.lon_lat_2_col_row(tif_file, right_lon, bot_lat)
    window = Window.from_slices((ul_row, lr_row), (ul_col, lr_col))
    return window


class crop_tif_file():
    def __init__(self):
        self.base_path = Path.cwd().absolute().parent
        self.dem_path = self.base_path / 'DEM_NDVI_data'
        self.comp_path = self.dem_path / 'sentinel_composite_geometry.tif'
        self.combined_pairs_base_path = self.base_path / 'Combined_Pairs/Path_114_12_combined'
        self.combined_linear_path = self.combined_pairs_base_path / 'linear_rates_combined_all.tif'
        self.comp_path_cropped = self.dem_path / 'sentinel_composite_geometry_cropped.tif'

    def crop_tif(self):
        left, bottom, right, top = rasterio_basic_functions.return_bounds(self.combined_linear_path)
        window = window_from_extent(self.comp_path, left, right, bottom, top)
        cropped_array_red = rasterio_basic_functions.tif_2_array_window_bands(self.comp_path, window, 1)
        cropped_array_green = rasterio_basic_functions.tif_2_array_window_bands(self.comp_path, window, 2)
        cropped_array_blue = rasterio_basic_functions.tif_2_array_window_bands(self.comp_path, window, 3)
        window_kwargs = rasterio_basic_functions.windowed_kwargs(self.comp_path, window)
        rasterio_basic_functions.write_reprojected_array_2_tif_bands(cropped_array_red,
                                                                     cropped_array_green,
                                                                     cropped_array_blue,
                                                                     self.comp_path_cropped,
                                                                     window_kwargs)

    def run_all(self):
        self.crop_tif()


if __name__ == '__main__':
    crop_tif_file_ = crop_tif_file()
    crop_tif_file_.run_all()
