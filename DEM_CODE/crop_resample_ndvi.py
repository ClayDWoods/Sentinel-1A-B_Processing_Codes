import sys, os

sys.path.append("..")
from helper_functions import rasterio_basic_functions
from helper_functions import general_functions as gen_fun
from pathlib import Path
import pathlib
import numpy as np
import math
from statistics import mean


class combine_all_dem_tifs():
    def __init__(self):
        self.base_path = Path.cwd().absolute().parent
        self.dem_path = self.base_path / 'DEM_NDVI_data'
        self.ndvi_path = self.dem_path / 'ndvi_mean_overall.tif'
        self.combined_pairs_base_path = self.base_path / 'Combined_Pairs/Path_114_12_combined'
        self.combined_linear_path = self.combined_pairs_base_path / 'linear_rates_combined_all.tif'
        self.final_ndvi_path = self.dem_path / 'ndvi_crop_resampled.tif'
        self.master_width = None
        self.master_height = None
        self.master_affine = None
        self.new_kwargs = None
        self.dem_save = None

    def crop_resample_ndvi(self):
        combined_linear_array = rasterio_basic_functions.tif_2_array(self.combined_linear_path)
        base_mask_array = np.where(combined_linear_array != 0,
                                   1,
                                   0)
        ndvi_transform = rasterio_basic_functions.get_tif_transform(self.ndvi_path)
        ndvi_open = rasterio_basic_functions.tif_2_array(self.ndvi_path)
        reproject_final_ndvi_array = (rasterio_basic_functions.
                                      reproject_tif_array_w_array(self.combined_linear_path,
                                                                  ndvi_open,
                                                                  ndvi_transform))
        final_ndvi_array = np.multiply(base_mask_array, reproject_final_ndvi_array)
        final_kwargs = rasterio_basic_functions.get_kwargs(self.combined_linear_path)
        rasterio_basic_functions.write_reprojected_array_2_tif(final_ndvi_array,
                                                               self.final_ndvi_path,
                                                               final_kwargs)

    def run_all(self):
        self.crop_resample_ndvi()


if __name__ == '__main__':
    combined_dems = combine_all_dem_tifs()
    combined_dems.run_all()
