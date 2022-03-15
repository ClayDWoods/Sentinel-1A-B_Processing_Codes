import os
import sys

sys.path.append(".")

from pathlib import Path
from Main.general_set_up_functions import rasterio_basic_functions
from Main.make_asc_dsc_header_txt import auto_find_ref_point_coherence
import math
from statistics import mean
import numpy as np


class gps_reference():
    def __init__(self, path_frame_1_reference_coh, path_frame_2_secondary_coh, num_ref):
        self.path_frame_reference_coh_list = [path_frame_1_reference_coh, path_frame_2_secondary_coh]
        self.linear_rate_arrays_use = None
        self.coh_arrays_use = None
        self.master_width = None
        self.master_height = None
        self.master_affine = None
        self.new_kwargs = None
        self.num_of_ref = num_ref
        self.half_size = 2
        self.coh_min = 0.15
        self.num_of_coh_blocks = 35
        self.num_of_max_to_check = 25
        self.max_cohs = None
        self.c_flag_rows = None
        self.c_flag_cols = None

    def define_master_resample_parameters(self):
        longitude_spacing_list, upper_left_corner_longitudes_list = [], []
        lattitude_spacing_list, upper_left_corner_lattitude_list = [], []
        lower_right_corner_longitude_list = []
        lower_right_corner_lattitude_list = []

        for asc_dsc in self.path_frame_reference_coh_list:
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
        min_longitude = min(upper_left_corner_longitudes_list) + 0.05
        max_lattitude = max(upper_left_corner_lattitude_list) + 0.05
        min_lattitude = min(lower_right_corner_lattitude_list) + 0.05

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
        coh_save = []
        for coh_tif in self.path_frame_reference_coh_list:
            coh_array_old = rasterio_basic_functions.tif_2_array(coh_tif)
            coh_transform_old = rasterio_basic_functions.get_tif_transform(coh_tif)
            base_crs = rasterio_basic_functions.get_tif_crs(coh_tif)
            base_kwargs = rasterio_basic_functions.get_kwargs(coh_tif)
            base_kwargs.update({'nodata': 0.,
                                'width': self.master_width,
                                'height': self.master_height,
                                'transform': self.master_affine
                                })
            if ii == 0:
                self.new_kwargs = base_kwargs
            ii = ii + 1
            coherence_reprojected_array = rasterio_basic_functions.reproject_tif_array(coh_tif,
                                                                                       self.master_height,
                                                                                       self.master_width,
                                                                                       coh_array_old,
                                                                                       coh_transform_old,
                                                                                       base_crs,
                                                                                       self.master_affine)
            coh_save.append(coherence_reprojected_array)
        self.coh_arrays_use = coh_save

    def set_up_arrays_and_reference(self):
        base_coh_array_1 = self.coh_arrays_use[0]
        base_coh_array_2 = self.coh_arrays_use[1]
        binary_coh_array_1 = np.where(base_coh_array_1 != 0,
                                      1,
                                      0)
        binary_coh_array_2 = np.where(base_coh_array_2 != 0,
                                      1,
                                      0)
        coh_sum_array = binary_coh_array_1 + binary_coh_array_2
        coh_mask_array = np.where(coh_sum_array == 2,
                                  1,
                                  0)
        base_coh_1_overlap = np.multiply(base_coh_array_1, coh_mask_array)
        base_coh_2_overlap = np.multiply(base_coh_array_2, coh_mask_array)
        base_coh_sum_overlap = (base_coh_1_overlap + base_coh_2_overlap) / 2
        self.max_cohs, self.c_flag_rows, self.c_flag_cols = (auto_find_ref_point_coherence.auto_find_reference_point
                                                             (base_coh_sum_overlap,
                                                              self.half_size,
                                                              self.num_of_ref,
                                                              self.coh_min,
                                                              self.num_of_coh_blocks,
                                                              self.num_of_max_to_check))
        lon, lat = (rasterio_basic_functions.col_row_2_lon_lat_waffine
                    (self.master_affine, self.c_flag_cols[0], self.c_flag_rows[0]))
        return lon, lat

    def run_gps_reference(self):
        self.define_master_resample_parameters()
        lon, lat = self.set_up_arrays_and_reference()
        return lon, lat
