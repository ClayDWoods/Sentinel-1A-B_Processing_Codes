import sys, os

sys.path.append("..")
from helper_functions import rasterio_basic_functions
from helper_functions import general_functions as gen_fun
from pathlib import Path
import pathlib
import numpy as np
import math
from statistics import mean


def listToString(s):
    # initialize an empty string
    str1 = " "

    # return string
    return str1.join(s)


class combine_all_dem_tifs():
    def __init__(self):
        self.base_path = Path.cwd().absolute().parent
        self.dem_path = self.base_path / 'DEM_combined'
        self.combined_pairs_base_path = self.base_path / 'Combined_Pairs/Path_114_12_combined'
        self.combined_linear_path = self.combined_pairs_base_path / 'linear_rates_combined_all.tif'
        self.dem_path_tifs = [str(x) for x in self.dem_path.glob('n*_e*_1arc_v3.tif')]
        self.dem_path_tifs.sort()
        self.final_dem_path = self.dem_path / 'combined_dem.tif'
        self.merged_dem_path = self.dem_path / 'merged_dem.tif'
        self.gdal_merge_exe = "/projects/clwo4142/software/anaconda/envs/mycustomenv/bin/gdal_merge.py"
        self.master_width = None
        self.master_height = None
        self.master_affine = None
        self.new_kwargs = None
        self.dem_save = None

    def combine_arrays(self):
        combined_linear_array = rasterio_basic_functions.tif_2_array(self.combined_linear_path)
        base_mask_array = np.where(combined_linear_array != 0,
                                   1,
                                   0)
        input_files = listToString([x.split('/')[-1] for x in self.dem_path_tifs])
        merge_cmd = (self.gdal_merge_exe + ' -o merged_dem.tif' + ' ' + input_files)
        gen_fun.subprocess(merge_cmd)
        dem_transform = rasterio_basic_functions.get_tif_transform(self.merged_dem_path)
        merged_dem_open = rasterio_basic_functions.tif_2_array(self.merged_dem_path)
        reproject_final_dem_rate_array = (rasterio_basic_functions.
                                          reproject_tif_array_w_array(self.combined_linear_path,
                                                                      merged_dem_open,
                                                                      dem_transform))
        final_dem_array = np.multiply(base_mask_array, reproject_final_dem_rate_array)
        final_kwargs = rasterio_basic_functions.get_kwargs(self.combined_linear_path)
        rasterio_basic_functions.write_reprojected_array_2_tif(final_dem_array,
                                                               self.final_dem_path,
                                                               final_kwargs)

    def run_all(self):
        self.combine_arrays()


if __name__ == '__main__':
    combined_dems = combine_all_dem_tifs()
    combined_dems.run_all()
