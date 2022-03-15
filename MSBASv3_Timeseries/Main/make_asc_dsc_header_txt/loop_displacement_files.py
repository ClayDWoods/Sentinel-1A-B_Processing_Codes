import sys

sys.path.append(".")

import pathlib
from scipy import stats as st
import numpy as np
from pathlib import Path
from Main.general_set_up_functions import general_functions, rasterio_basic_functions
from Main.make_asc_dsc_header_txt import snaphu_log_parser
import multiprocessing
import math


def return_overall_mean_std_from_files(file_list):
    mean_list = []
    std_list = []
    for file in file_list:
        test_array = rasterio_basic_functions.tif_2_array(file).flatten()
        test_array_nz = test_array[test_array != 0]
        test_array_mean = np.mean(test_array_nz)
        test_array_std = np.std(test_array_nz)
        mean_list.append(test_array_mean)
        std_list.append(test_array_std)
    return mean_list, std_list


def return_overall_mean_from_files(file_list):
    mean_list = []
    for file in file_list:
        test_array = rasterio_basic_functions.tif_2_array(file).flatten()
        test_array_nz = test_array[test_array != 0]
        test_array_mean = np.mean(test_array_nz)
        mean_list.append(test_array_mean)
    mean_array = np.array(mean_list)
    return mean_array


def make_function_parallel(func_, file_list, num_cores):
    pool = multiprocessing.Pool()
    num_of_items = len(file_list)
    block_size = math.ceil(num_of_items / num_cores)
    block_starts = [x for x in range(0, num_of_items, block_size)]
    file_lists = [file_list[x:(x + block_size)] for x in block_starts]
    vals = pool.map(func_, file_lists)
    pool.close()
    pool.join()
    return vals


def find_std_return_lowest(path_1, path_2, mcoh_all, mcoh_gacos_all, num_cores):
    # vals_all_1 = make_function_parallel(return_overall_mean_from_files, mcoh_all, num_cores)
    # vals_all_2 = make_function_parallel(return_overall_mean_from_files, mcoh_gacos_all, num_cores)
    # means_array_1 = np.array([])
    # means_array_2 = np.array([])
    # for values_1, values_2 in zip(vals_all_1, vals_all_2):
    #     means_array_1 = np.append(means_array_1, values_1.flatten())
    #     means_array_2 = np.append(means_array_2, values_2.flatten())
    # mean_mcoh_all = np.mean(means_array_1)
    # mean_mcoh_gacos_all = np.mean(means_array_2)
    # central_mean = (mean_mcoh_all + mean_mcoh_gacos_all) / 2
    # central_mean = mean_mcoh_all
    # central_mean = mean_mcoh_gacos_all
    path_1_array = rasterio_basic_functions.tif_2_array(path_1)
    path_2_array = rasterio_basic_functions.tif_2_array(path_2)
    path_1_array_nz = path_1_array[path_1_array != 0]
    path_2_array_nz = path_2_array[path_2_array != 0]
    path_1_std = np.std(path_1_array_nz)
    path_2_std = np.std(path_2_array_nz)
    # path_1_mean = np.mean(path_1_array_nz)
    # path_2_mean = np.mean(path_2_array_nz)
    # path_1_mean_diff = abs(path_1_mean - central_mean)
    # path_2_mean_diff = abs(path_2_mean - central_mean)
    # if path_1_mean_diff < path_2_mean_diff:
    #     return path_1
    # else:
    #     return path_2
    if path_1_std < path_2_std:
        return path_1
    else:
        return path_2


def find_outlier_distributions(df_to_update, file_list, num_of_stds, num_cores, std_E):
    df_ = df_to_update
    rejected_files = []
    mean_rejected_files = []
    vals_all = make_function_parallel(return_overall_mean_std_from_files, file_list, num_cores)
    means_array = np.array([])
    std_array = np.array([])
    for values in vals_all:
        means_, stds_ = values[0], values[1]
        means_array = np.append(means_array, means_)
        std_array = np.append(std_array, stds_)
    ###
    std_mean_list = np.std(means_array)
    mean_list_mean = np.mean(means_array)
    ###
    std_of_std_list = np.std(std_array)
    mean_of_std_list = np.mean(std_array)
    ###
    lower_bound_mean = mean_list_mean - num_of_stds * std_mean_list
    upper_bound_mean = mean_list_mean + num_of_stds * std_mean_list
    lower_bound_std = mean_of_std_list - num_of_stds * std_of_std_list
    upper_bound_std = mean_of_std_list + num_of_stds * std_of_std_list

    if std_E == 'std':
        for file in file_list:
            test_array = rasterio_basic_functions.tif_2_array(file).flatten()
            test_array_nz = test_array[test_array != 0]
            test_array_mean = np.mean(test_array_nz)
            if test_array_mean < lower_bound_mean or test_array_mean > upper_bound_mean:
                df_ = df_[df_[0] != file]
                rejected_files.append(file)
                mean_rejected_files.append(test_array_mean)
        new_list_means_rejected = [x for x in file_list if x not in rejected_files]
        j = 0
        for file_ in new_list_means_rejected:
            test_array = rasterio_basic_functions.tif_2_array(file_).flatten()
            test_array_nz = test_array[test_array != 0]
            test_array_std_ = np.std(test_array_nz)
            if test_array_std_ < lower_bound_std or test_array_std_ > upper_bound_std:
                df_ = df_[df_[0] != file_]
                if j == 0:
                    rejected_files.append('std')
                    mean_rejected_files.append(0)
                rejected_files.append(file_)
                mean_rejected_files.append(test_array_std_)
                j = j + 1
    return df_, rejected_files, mean_rejected_files


def loop_displacement_files_mcoh(displacement_tif_pairs_mcoh_list, base_path, base_path_pair, msbasv3_path,
                                 flag, num_of_stds, asc_dsc_mcoh_list, asc_dsc_mcoh_gacos_list, num_cores,
                                 gacos_or_non_gacos_path, gacos_or_non_full, std_E):
    rel_mcoh_path_list, baseline_length_list, ref_date_list, second_date_list = [], [], [], []
    gacos_or_non_dict = general_functions.open_json_file(gacos_or_non_gacos_path)
    gacos_or_non_gacos_str_ = [str(x) for x in gacos_or_non_dict['g_ng']]
    for pair_rcut, gacos_or_non_gacos_str in zip(displacement_tif_pairs_mcoh_list, gacos_or_non_gacos_str_):
        reference_date = pair_rcut.name.split('_')[0]
        secondary_date = pair_rcut.name.split('_')[1]
        relative_path_mcoh = pair_rcut.relative_to(base_path)
        relative_path_str_mcoh = str(Path('../' + str(relative_path_mcoh)))
        if gacos_or_non_full == 'gacos_mixed':
            if gacos_or_non_gacos_str == 'g':
                relative_path_str_new = relative_path_str_mcoh.split('.tif')[0] + '_mcoh_gacos.tif'
            else:
                relative_path_str_new = relative_path_str_mcoh.split('.tif')[0] + '_mcoh.tif'
            rel_mcoh_path_list.append(str(relative_path_str_new))
        elif gacos_or_non_full == 'gacos_full':
            if gacos_or_non_gacos_str == 'g':
                relative_path_str_new = relative_path_str_mcoh.split('.tif')[0] + '_mcoh_gacos.tif'
                rel_mcoh_path_list.append(str(relative_path_str_new))
        else:
            if gacos_or_non_gacos_str == 'g':
                relative_path_str_new = relative_path_str_mcoh.split('.tif')[0] + '_mcoh.tif'
                rel_mcoh_path_list.append(str(relative_path_str_new))
        snaphu_base = base_path_pair / ('Main/Snaphu_Logs')
        snaphu_pair_log_path = snaphu_base / ('snaphu_' + reference_date + '_' + secondary_date + '.log')
        base_line_length = snaphu_log_parser.snaphu_log_parser_f(snaphu_pair_log_path)
        baseline_length_list.append(str(base_line_length))
        ref_date_list.append(str(reference_date))
        second_date_list.append(str(secondary_date))
    pair_name_mcoh = msbasv3_path / (flag + '.txt')
    list_of_files_used_name = msbasv3_path / (flag + '_file_list_msbasv3.json')
    list_of_rejected_files_name = msbasv3_path / (flag + '_outliers.json')
    df_asc_dsc = general_functions.df_from_lists(rel_mcoh_path_list, baseline_length_list, ref_date_list,
                                                 second_date_list)
    rel_mcoh_path_list.sort()
    df_updated, rejected_files, mean_rejected_files = (find_outlier_distributions
                                                       (df_asc_dsc, rel_mcoh_path_list, num_of_stds, num_cores, std_E))
    mean_rejected_files_str = [str(x) for x in mean_rejected_files]
    rejected_files_dict = {'outlier_int_mean': rejected_files, 'rejected_mean': mean_rejected_files_str}
    general_functions.write_json_from_dict(rejected_files_dict, list_of_rejected_files_name)
    rel_mcoh_path_list_updated = [x for x in rel_mcoh_path_list if x not in rejected_files]
    rel_mcoh_path_list_updated_dict = {'rel_coh_path_list_updated': rel_mcoh_path_list_updated}
    general_functions.write_json_from_dict(rel_mcoh_path_list_updated_dict, list_of_files_used_name)
    general_functions.save_df(df_updated, pair_name_mcoh)
