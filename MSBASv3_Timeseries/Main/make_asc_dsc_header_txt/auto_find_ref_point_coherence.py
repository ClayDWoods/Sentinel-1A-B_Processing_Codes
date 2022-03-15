import sys

sys.path.append(".")

import numpy as np
import pandas as pd
from Main.general_set_up_functions import rasterio_basic_functions
import pathlib


def auto_find_reference_point(coherence_average_tif_path, half_size, num_of_ref, coh_min,
                              num_of_ref_coh_blocks, num_of_max_to_check):
    try:
        if isinstance(coherence_average_tif_path, pathlib.PurePath):
            coh_array = rasterio_basic_functions.tif_2_array(coherence_average_tif_path)
            array_width = rasterio_basic_functions.get_tif_width(coherence_average_tif_path)
            array_height = rasterio_basic_functions.get_tif_height(coherence_average_tif_path)
        else:
            coh_array = coherence_average_tif_path
            array_width = coherence_average_tif_path.shape[1]
            array_height = coherence_average_tif_path.shape[0]
        block_half_size = half_size
        y = range(block_half_size - 1, array_width, (2 * block_half_size + 1))
        x = range(block_half_size - 1, array_height, (2 * block_half_size + 1))
        xv, yv = np.meshgrid(x, y, sparse=False, indexing='ij')
        means_list = []
        rows_list = []
        cols_list = []
        for i in range(len(x)):
            for j in range(len(y)):
                row = xv[i, j]
                col = yv[i, j]
                row_start = row - block_half_size + 1
                row_end = row + block_half_size + 1
                col_start = col - block_half_size + 1
                col_end = col + block_half_size + 1
                if row_start < 0:
                    row_start = 0
                if col_start < 0:
                    col_start = 0
                coh_slice_array = coh_array[row_start:row_end, col_start:col_end]
                sliced_coh_average = np.mean(coh_slice_array)
                means_list.append(sliced_coh_average)
                rows_list.append(row)
                cols_list.append(col)
        means_list_coh_use = []
        rows_cut_list = []
        cols_cut_list = []
        avg_coh_sort = (pd.DataFrame(list(zip(means_list, rows_list, cols_list)))
                        .sort_values(by=0).tail(int(num_of_ref_coh_blocks)))
        ii = 1
        for index, row in avg_coh_sort.iterrows():
            v_to_list = row.values.tolist()
            row_u = int(v_to_list[1])
            col_u = int(v_to_list[2])
            row_start = int(row_u - block_half_size + 1)
            row_end = int(row_u + block_half_size + 1)
            col_start = int(col_u - block_half_size + 1)
            col_end = int(col_u + block_half_size + 1)
            if row_start < 0:
                row_start = 0
            if col_start < 0:
                col_start = 0
            #
            coh_slice_array = coh_array[row_start:row_end, col_start:col_end]
            ind = np.unravel_index(np.argsort(coh_slice_array, axis=None), coh_slice_array.shape)
            ind_row = ind[0][-ii]
            ind_col = ind[1][-ii]
            ind_row_use = row_start + ind_row
            ind_col_use = col_start + ind_col
            row_start_use = ind_row_use - block_half_size + 1
            row_end_use = ind_row_use + block_half_size + 1
            col_start_use = ind_col_use - block_half_size + 1
            col_end_use = ind_col_use + block_half_size + 1
            if row_start_use < 0:
                row_start_use = 0
            if col_start_use < 0:
                col_start_use = 0
            coh_slice_use_array = coh_array[row_start_use:row_end_use, col_start_use:col_end_use]
            sliced_coh_average_use = np.mean(coh_slice_use_array)
            if sliced_coh_average_use > coh_min:
                means_list_coh_use.append(sliced_coh_average_use)
                rows_cut_list.append(ind_row_use)
                cols_cut_list.append(ind_col_use)
        max_list_w_ind = (pd.DataFrame(list(zip(means_list_coh_use, rows_cut_list, cols_cut_list)))
                          .sort_values(by=0).tail(int(num_of_ref)))
        max_cohs_list = max_list_w_ind[0].values.tolist()
        max_rows_list = max_list_w_ind[1].values.tolist()
        max_cols_list = max_list_w_ind[2].values.tolist()
        return max_cohs_list, max_rows_list, max_cols_list
    except KeyError:
        print('No valid auto_reference regions found --- trying alternate algorithm')
        print('change parameters if want to find the regions')
        coh_array = rasterio_basic_functions.tif_2_array(coherence_average_tif_path)
        max_list_coh = []
        max_col_list = []
        max_row_list = []
        for ii in range(1, num_of_max_to_check + 1):
            ind = np.unravel_index(np.argsort(coh_array, axis=None), coh_array.shape)
            ind_row = ind[0][-ii]
            ind_col = ind[1][-ii]
            row_start = ind_row - half_size
            row_end = ind_row + half_size
            col_start = ind_col - half_size
            col_end = ind_col + half_size
            sliced_coh_average = np.mean(coh_array[row_start:row_end, col_start:col_end])
            max_list_coh.append(sliced_coh_average)
            max_row_list.append(ind_row)
            max_col_list.append(ind_col)
        max_list_w_ind = (pd.DataFrame(list(zip(max_list_coh, max_row_list, max_col_list)))
                          .sort_values(by=0).tail(int(num_of_ref)))
        max_cohs_list = max_list_w_ind[0].values.tolist()
        max_rows_list = max_list_w_ind[1].values.tolist()
        max_cols_list = max_list_w_ind[2].values.tolist()
        return max_cohs_list, max_rows_list, max_cols_list
