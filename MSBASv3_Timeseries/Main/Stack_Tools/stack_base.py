import sys
sys.path.append(".")

from Main.general_set_up_functions import rasterio_basic_functions as rbf
import math
import numpy as np
from pathlib import Path
from Main.general_set_up_functions import rasterio_basic_functions

def filter_asc_dsc_list(file_list_base, file_list_base_edit):
    test_paths_ = [str(x) for x in file_list_base_edit]
    test_paths_.sort()

    temp_list = []
    for file in file_list_base:
        date_ = file.name.split('_')[0]
        temp_list.append(date_)
    set_list = list(dict.fromkeys(temp_list))

    matching_list = []
    for date__ in set_list:
        matching = [s for s in test_paths_ if ('/' + date__ + '_') in s]
        matching_list.append(Path(matching[0]))
    matching_list.sort()
    return matching_list


def make_stack_base_asc_dsc(file_list, incidence, c_flag_col, c_flag_row, c_size, coh_mask_tif, num_of_ref, coh_thresh):
    coh_array = rasterio_basic_functions.tif_2_array(coh_mask_tif)
    coherence_mask_apply = np.where(coh_array < coh_thresh,
                                    0,
                                    1)
    half_size_c_flag = c_size
    if c_flag_col and c_flag_row is not None:
        if num_of_ref == 1:
            row_start = int(c_flag_row) - half_size_c_flag
            row_end = int(c_flag_row) + half_size_c_flag + 1
            col_start = int(c_flag_col) - half_size_c_flag
            col_end = int(c_flag_col) + half_size_c_flag + 1
        else:
            for col, row in zip(c_flag_col, c_flag_row):
                row_start_l = []
                row_end_l = []
                col_start_l = []
                col_end_l = []
                row_start = int(row) - half_size_c_flag
                row_end = int(row) + half_size_c_flag + 1
                col_start = int(col) - half_size_c_flag
                col_end = int(col) + half_size_c_flag + 1
                row_start_l.append(row_start)
                row_end_l.append(row_end)
                col_start_l.append(col_start)
                col_end_l.append(col_end)
    mul_factor = np.cos(math.radians(incidence))
    tif_base_file = file_list[0]
    tif_width = rbf.get_tif_width(tif_base_file)
    tif_height = rbf.get_tif_height(tif_base_file)
    base_vertical_array = np.zeros((tif_height, tif_width))
    num_of_images = len(file_list)
    ref_loc_mean_list = []
    for asc_dsc in file_list:
        tif_array = rbf.tif_2_array(asc_dsc)
        tif_array_vertical = tif_array/mul_factor
        if c_flag_col and c_flag_row is not None:
            if num_of_ref == 1:
                tif_array_vertical_ref_subtract = np.mean(tif_array_vertical[row_start:row_end, col_start:col_end])
                ref_loc_mean_list.append(tif_array_vertical_ref_subtract)
            else:
                for row_start, row_end, col_start, col_end in zip(row_start_l, row_end_l, col_start_l, col_end_l):
                    tif_array_vertical_ref_subtract = np.mean(tif_array_vertical[row_start:row_end, col_start:col_end])
                    ref_loc_mean_list.append(tif_array_vertical_ref_subtract)
            tif_array_vertical_ref_mean_subtract = np.mean(ref_loc_mean_list)
            tif_array_vertical_referenced = tif_array_vertical - tif_array_vertical_ref_mean_subtract
            base_vertical_array = base_vertical_array + tif_array_vertical_referenced
        else:
            base_vertical_array = base_vertical_array + tif_array_vertical

    tif_array_vertical_stack = base_vertical_array/num_of_images
    tif_array_vertical_stack_use = np.multiply(tif_array_vertical_stack, coherence_mask_apply).astype('float32')
    return tif_array_vertical_stack_use

