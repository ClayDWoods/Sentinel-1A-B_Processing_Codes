import os
import shutil
import sys
from pathlib import Path
from general_set_up_functions import rasterio_basic_functions, general_functions
import math
from statistics import mean
import numpy as np
import multiprocessing
import matplotlib.pyplot as plt


def find_nearest(array, value):
    absolute_val_array = np.abs(array - value)
    smallest_difference_index = absolute_val_array.argmin()
    return smallest_difference_index


def check_pixel_diff(diff_array_pre, base_linear_rate_2_overlap, base_linear_rate_1_overlap, base_coh_overlap):
    diff_values_check = []
    base_diff_values = []
    for base_diff_v_use in diff_array_pre:
        base_linear_rate_plus_diff = (general_functions.add_constant_2_array_ignoring_zeros
                                      (base_diff_v_use, base_linear_rate_2_overlap))
        diff_base_1_base_2 = base_linear_rate_1_overlap - base_linear_rate_plus_diff
        diff_base_1_base_2_w = np.multiply(diff_base_1_base_2, base_coh_overlap)
        diff_base_1_base_2_flat_nonzero = diff_base_1_base_2_w[diff_base_1_base_2_w != 0]
        sum_sq = np.sum(np.square(diff_base_1_base_2_flat_nonzero))
        diff_values_check.append(sum_sq)
        base_diff_values.append(base_diff_v_use)
    return diff_values_check, base_diff_values


def overlap_hist_plot(save_fig, data_1, data_2, path_1, path_2):
    data_1_no_zero = data_1[data_1 != 0]
    data_2_no_zero = data_2[data_2 != 0]
    plt.figure(figsize=(8, 6))
    plt.hist(data_1_no_zero, bins=100, alpha=0.5, label=path_1)
    plt.hist(data_2_no_zero, bins=100, alpha=0.5, label=path_2)
    plt.xlabel("value", size=14)
    plt.ylabel("Count", size=14)
    plt.title("overlap comparison plot")
    plt.legend(loc='upper right')
    plt.savefig(save_fig)


def single_hist_plot(save_fig, data_1, path_1, path_2, std_, mean_):
    data_1_no_zero = data_1[data_1 != 0]
    path_name = path_1 + '_' + path_2 + '_diff'
    plt.figure(figsize=(8, 6))
    plt.hist(data_1_no_zero, bins=100, alpha=0.5, label=path_name)
    plt.text(0.1, 0.9, ('mean=' + str(mean_)), fontsize=12)
    plt.text(0.1, 0.8, ('std=' + str(std_)), fontsize=12)
    plt.xlabel("value", size=14)
    plt.ylabel("Count", size=14)
    plt.title("diff_plot")
    plt.legend(loc='upper right')
    plt.savefig(save_fig)


class gps_reference():
    def __init__(self, path_frame_1_reference_rate, path_frame_1_reference_coh,
                 path_frame_2_secondary_rate, path_frame_2_secondary_coh, cores):
        self.num_cores = cores
        self.base_path = Path.cwd().absolute().parent
        self.combine_pairs_path = self.base_path / 'Combined_Pairs/Path_114_12_referenced'
        self.referenced_pairs_save_path = self.combine_pairs_path / 'referenced_path_frame'
        self.overlap_paths_save = self.combine_pairs_path / 'overlap_tifs'
        self.path_frame_1_reference_rate = ([Path(str(x)) for x in
                                             self.combine_pairs_path.glob(path_frame_1_reference_rate)])[0]
        self.path_frame_2_reference_rate = ([Path(str(x)) for x in
                                             self.combine_pairs_path.glob(path_frame_2_secondary_rate)])[0]
        self.path_frame_1_coh = ([Path(str(x)) for x in
                                  self.combine_pairs_path.glob(path_frame_1_reference_coh)])[0]
        self.path_frame_2_coh = ([Path(str(x)) for x in
                                  self.combine_pairs_path.glob(path_frame_2_secondary_coh)])[0]
        self.path_pair_referenced_name = (self.path_frame_2_reference_rate.parent /
                                          (self.path_frame_2_reference_rate.stem + '_ref.tif'))
        self.ref_point_location_name = (self.path_frame_2_reference_rate.parent /
                                        (self.path_frame_2_reference_rate.stem + '_ref.txt'))
        self.path_frame_reference_1 = '_'.join(self.path_frame_1_reference_rate.stem.split('_')[2:6])
        self.path_frame_reference_2 = '_'.join(self.path_frame_2_reference_rate.stem.split('_')[2:6])
        self.hist_save_name = self.combine_pairs_path / (self.path_frame_reference_1 + '_' +
                                                         self.path_frame_reference_2 +
                                                         '_hist.png')
        self.hist_save_name_diff = self.combine_pairs_path / (self.path_frame_reference_1 + '_' +
                                                              self.path_frame_reference_2 +
                                                              '_diff_hist.png')
        self.overlap_save_tif_path_1 = self.overlap_paths_save / (self.path_frame_reference_1 + '_' +
                                                                  self.path_frame_reference_2 +
                                                                  '_1.tif')
        self.overlap_save_tif_path_2 = self.overlap_paths_save / (self.path_frame_reference_1 + '_' +
                                                                  self.path_frame_reference_2 +
                                                                  '_2.tif')
        self.overlap_save_tif_diff = self.overlap_paths_save / (self.path_frame_reference_1 + '_' +
                                                                self.path_frame_reference_2 +
                                                                '_diff.tif')
        self.path_frame_reference_rate_list = [self.path_frame_1_reference_rate, self.path_frame_2_reference_rate]
        self.path_frame_reference_coh_list = [self.path_frame_1_coh, self.path_frame_2_coh]
        self.linear_rate_arrays_use = None
        self.coh_arrays_use = None
        self.master_width = None
        self.master_height = None
        self.master_affine = None
        self.new_kwargs = None
        self.num_of_ref = 15
        self.half_size = 2
        self.coh_min = 0.15
        self.num_of_coh_blocks = 35
        self.num_of_max_to_check = 25

    def define_master_resample_parameters(self):
        longitude_spacing_list, upper_left_corner_longitudes_list = [], []
        lattitude_spacing_list, upper_left_corner_lattitude_list = [], []
        lower_right_corner_longitude_list = []
        lower_right_corner_lattitude_list = []

        for asc_dsc in self.path_frame_reference_rate_list:
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
        for linear_rate_tif, coh_tif in zip(self.path_frame_reference_rate_list, self.path_frame_reference_coh_list):
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

    def set_up_arrays_and_reference(self):
        base_array_2_reference = rasterio_basic_functions.tif_2_array(self.path_frame_2_reference_rate)
        base_linear_rate_1 = self.linear_rate_arrays_use[0]
        base_linear_rate_2 = self.linear_rate_arrays_use[1]
        base_coh_array_1 = self.coh_arrays_use[0]
        base_coh_array_2 = self.coh_arrays_use[1]
        base_linear_array_mask_1 = np.where(base_linear_rate_1 != 0,
                                            1,
                                            0)
        base_linear_array_mask_2 = np.where(base_linear_rate_2 != 0,
                                            1,
                                            0)
        base_linear_array_mask = base_linear_array_mask_1 + base_linear_array_mask_2
        base_linear_array_mask_use = np.where(base_linear_array_mask != 2,
                                              0,
                                              1)
        base_coh_1_overlap = np.multiply(base_coh_array_1, base_linear_array_mask_use)
        base_coh_2_overlap = np.multiply(base_coh_array_2, base_linear_array_mask_use)
        base_coh_sum_overlap = (base_coh_1_overlap + base_coh_2_overlap) / 2
        base_linear_rate_1_overlap = np.multiply(base_linear_rate_1, base_linear_array_mask_use)
        base_linear_rate_2_overlap = np.multiply(base_linear_rate_2, base_linear_array_mask_use)
        base_linear_rate_1_overlap_flat_nozero = base_linear_rate_1_overlap[base_linear_rate_1_overlap != 0]
        base_linear_rate_2_overlap_flat_nozero = base_linear_rate_2_overlap[base_linear_rate_2_overlap != 0]
        base_diff_array_use = base_linear_rate_1_overlap_flat_nozero - base_linear_rate_2_overlap_flat_nozero
        min_diff = np.min(base_diff_array_use)
        max_diff = np.max(base_diff_array_use)
        temp_array_use_diff = np.arange(min_diff, max_diff, (max_diff - min_diff) / 5000)
        #
        pool = multiprocessing.Pool()
        num_of_elements = len(base_linear_rate_1_overlap_flat_nozero)
        block_size = int(num_of_elements / int(self.num_cores))
        block_starts = [x for x in range(0, num_of_elements, block_size)]

        numpy_array_list_1 = [temp_array_use_diff[x:(x + block_size)] for x in block_starts]
        diff_values_all = ([pool.apply(check_pixel_diff,
                                       args=(base_array_1, base_linear_rate_2_overlap,
                                             base_linear_rate_1_overlap, base_coh_sum_overlap)) for
                            base_array_1 in numpy_array_list_1])
        pool.close()
        pool.join()

        diff_values_check = np.array([])
        base_diff_values = np.array([])
        for values in diff_values_all:
            values_diff_values_check, base_diff_values_ = values[0], values[1]
            diff_values_check = np.append(diff_values_check, values_diff_values_check)
            base_diff_values = np.append(base_diff_values, base_diff_values_)
        diff_values_check_list = list(diff_values_check)
        diff_values_check_list_min = min(diff_values_check_list)
        min_index_u = diff_values_check_list.index(diff_values_check_list_min)
        base_diff_values_list = list(base_diff_values)
        min_diff = base_diff_values_list[min_index_u]
        min_index = find_nearest(base_diff_array_use, min_diff)
        best_dif_value = base_diff_array_use[min_index]
        base_1_v_use = base_linear_rate_1_overlap_flat_nozero[min_index]
        index_ = np.where(base_linear_rate_1_overlap == base_1_v_use)
        index_col = index_[1][0]
        index_row = index_[0][0]
        base_array_referenced_2_main_array = (general_functions.
                                              add_constant_2_array_ignoring_zeros(best_dif_value,
                                                                                  base_array_2_reference))
        base_array_overlap_referenced_2_main_array = (general_functions.
                                                      add_constant_2_array_ignoring_zeros(best_dif_value,
                                                                                          base_linear_rate_2_overlap))
        overlap_diff_array = base_linear_rate_1_overlap - base_array_overlap_referenced_2_main_array
        mean_overlap_diff = np.mean(overlap_diff_array)
        std_overlap_diff = np.std(overlap_diff_array)
        single_hist_plot(self.hist_save_name_diff, overlap_diff_array, self.path_frame_reference_1,
                         self.path_frame_reference_2, std_overlap_diff, mean_overlap_diff)
        overlap_hist_plot(self.hist_save_name, base_linear_rate_1_overlap,
                          base_array_overlap_referenced_2_main_array,
                          self.path_frame_reference_1, self.path_frame_reference_2)
        lon, lat = rasterio_basic_functions.col_row_2_lon_lat_waffine(self.master_affine, index_col, index_row)
        base_ref_kwargs = rasterio_basic_functions.get_kwargs(self.path_frame_2_reference_rate)
        rasterio_basic_functions.write_reprojected_array_2_tif(base_array_referenced_2_main_array,
                                                               self.path_pair_referenced_name,
                                                               base_ref_kwargs)
        rasterio_basic_functions.write_reprojected_array_2_tif(base_linear_rate_1_overlap,
                                                               self.overlap_save_tif_path_1,
                                                               self.new_kwargs)
        rasterio_basic_functions.write_reprojected_array_2_tif(base_array_overlap_referenced_2_main_array,
                                                               self.overlap_save_tif_path_2,
                                                               self.new_kwargs)
        rasterio_basic_functions.write_reprojected_array_2_tif(overlap_diff_array,
                                                               self.overlap_save_tif_diff,
                                                               self.new_kwargs)
        shutil.copy(self.path_pair_referenced_name, self.referenced_pairs_save_path)
        txt_str = [str(lon), str(lat), str(best_dif_value)]

        general_functions.write_txt_file(self.ref_point_location_name, txt_str)

    def run_gps_reference(self):
        self.define_master_resample_parameters()
        self.set_up_arrays_and_reference()


if __name__ == '__main__':
    sys_index_var_1 = sys.argv[1]
    sys_index_var_2 = sys.argv[2]
    sys_index_var_3 = sys.argv[3]
    sys_index_var_4 = sys.argv[4]
    sys_index_var_5 = sys.argv[5]
    gps_reference = gps_reference(sys_index_var_1, sys_index_var_2, sys_index_var_3, sys_index_var_4, sys_index_var_5)
    gps_reference.run_gps_reference()
