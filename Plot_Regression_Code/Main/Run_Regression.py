import os
import shutil
import sys

import numpy as np

sys.path.append("..")

from pathlib import Path
from general_functions_all import rasterio_basic_functions, general_functions
from regression_functions import regression_do

import multiprocessing


class regression_and_time_series():
    def __init__(self, path_1_name, path_2_name, header_dir_name, mode_, EW_UD):
        self.EW_UD = EW_UD
        self.path_1_name = path_1_name
        self.path_2_name = path_2_name
        self.regression_mode = mode_
        self.base_path = Path.cwd().absolute().parent
        self.base_plots_path = self.base_path / 'Plots_Regression'
        self.regression_tif_plots_path = self.base_plots_path / 'Linear_Rate_Spatial_Plots'
        self.kml_tif_plots_path = self.base_plots_path / 'Linear_Rate_Spatial_Plots_kml'
        self.numpy_tif_arrays_path = self.base_plots_path / 'Binary_Spatial_Arrays'
        self.tif_meta_data_path = self.base_plots_path / 'tif_meta_data'
        self.time_series_plots_path = self.base_plots_path / 'Time_Series_Plots'
        self.reference_points_file = self.base_plots_path / 'Reference_point/reference_point.json'
        self.coherence_tif_path = self.base_plots_path / 'coherence_mask/coherence_mask.tif'
        self.path_1_kml_path = self.base_plots_path / (path_1_name + '_kml_coords/map-overlay.kml')
        self.path_2_kml_path = self.base_plots_path / (path_2_name + '_kml_coords/map-overlay.kml')
        self.header_path = self.base_plots_path / header_dir_name
        self.path_2_msbas_tsout = self.header_path / 'MSBAS_TSOUT.txt'
        self.regression_data_used_path = self.base_plots_path / 'Data_used_regression'
        if self.EW_UD == 'UD':
            self.msbas_tif_files = [str(x) for x in self.header_path.glob('MSBAS_2*_UD.tif')]
        else:
            self.msbas_tif_files = [str(x) for x in self.header_path.glob('MSBAS_2*_EW.tif')]
        self.msbas_tif_files.sort()
        self.num_of_cores = 20
        self.years_1 = ['2014', '2015', '2016', '2017', '2018', '2019', '2020', '2021', '2022']
        self.months_1 = ['05', '06', '07', '08', '09', '10']
        self.years_2 = ['2021']
        self.months_2 = ['04', '05', '06', '07', '08', '09', '10']
        self.bucket_size = None
        self.numpy_array_paths = None
        self.kwargs_path = None
        self.base_kwargs = None
        self.msbasv3_json_save_name = None
        self.msbasv3_tsout_dict = None
        self.tsout_dates_yyyy_mm = None
        self.tsout_dates_yyyy_yyy = None
        self.msbas_names = None
        self.tif_height = None
        self.tif_width = None
        self.final_array_use = None
        self.dates_yyyy_mm_list_updated = None
        self.dates_yyyy_yyy_list_updated = None
        self.msbas_np_names_updated = None
        self.time_yyyy_yyy_use = None
        self.regression_tif_save_name = None
        self.coh_save_name = None

    def make_dirs_(self):
        if not self.regression_tif_plots_path.exists():
            os.mkdir(self.regression_tif_plots_path)
        if not self.kml_tif_plots_path.exists():
            os.mkdir(self.kml_tif_plots_path)
        if not self.numpy_tif_arrays_path.exists():
            os.mkdir(self.numpy_tif_arrays_path)
        if not self.tif_meta_data_path.exists():
            os.mkdir(self.tif_meta_data_path)
        if not self.time_series_plots_path.exists():
            os.mkdir(self.time_series_plots_path)
        if not self.regression_data_used_path.exists():
            os.mkdir(self.regression_data_used_path)

    def make_numpy_arrays(self):
        tif_kwargs = rasterio_basic_functions.get_kwargs(self.msbas_tif_files[0])
        self.base_kwargs = tif_kwargs
        tif_kwargs_crs_str = str(tif_kwargs['crs'])
        tif_kwargs['crs'] = tif_kwargs_crs_str
        tif_kwargs_affine_list = list(tif_kwargs['transform'])[0:6]
        tif_kwargs['transform'] = tif_kwargs_affine_list
        tif_kwargs_name = self.tif_meta_data_path / 'kwargs_tif.json'
        general_functions.write_json_from_dict(tif_kwargs, tif_kwargs_name)
        for tif in self.msbas_tif_files:
            tif_path = Path(tif)
            tif_name_np = tif_path.name.split('.')[0] + '.npy'
            tif_numpy_binary_name = self.numpy_tif_arrays_path / tif_name_np
            tif_array = rasterio_basic_functions.tif_2_array(tif)
            general_functions.write_array_2_npy(tif_array, tif_numpy_binary_name)
        if self.EW_UD == 'UD':
            self.numpy_array_paths = [str(x) for x in self.numpy_tif_arrays_path.glob('MSBAS_2*_UD.npy')]
        else:
            self.numpy_array_paths = [str(x) for x in self.numpy_tif_arrays_path.glob('MSBAS_2*_EW.npy')]
        self.numpy_array_paths.sort()
        self.kwargs_path = tif_kwargs_name

    def read_tsout_and_parse(self):
        msbasv3_tsout = general_functions.open_txt_file(self.path_2_msbas_tsout)
        lines_append = []
        for line in msbasv3_tsout:
            lines_append.append(line.split(' '))
        date_list = []
        date_yyyy_yyy_list = []
        msbas_file_name_list = []
        for line_ in lines_append[3:]:
            date_ = line_[0]
            date_yyyy_yyy = line_[1]
            if self.EW_UD == 'UD':
                msbas_file_name = line_[3].split('.')[0] + '.npy'
            else:
                msbas_file_name = line_[2].split('.')[0] + '.npy'
            date_list.append(date_)
            date_yyyy_yyy_list.append(date_yyyy_yyy)
            msbas_file_name_list.append(msbas_file_name)
        dict_lists = {'mode': self.regression_mode, 'YYYYMMDDTHHMMSS': date_list, 'YYYY.YYY': date_yyyy_yyy_list,
                      'msbas_file_names': msbas_file_name_list,
                      'numpy_paths': self.numpy_array_paths}
        self.msbasv3_tsout_dict = dict_lists
        self.tsout_dates_yyyy_mm = date_list
        self.tsout_dates_yyyy_yyy = date_yyyy_yyy_list
        self.msbas_names = msbas_file_name_list
        if self.regression_mode == 'all':
            self.msbasv3_json_save_name = (self.regression_data_used_path /
                                           ('MSBAS_TSOUT_' + self.regression_mode + '.json'))
            general_functions.write_json_from_dict(dict_lists, self.msbasv3_json_save_name)

    def make_large_array_for_regression(self):
        if self.regression_mode == 'all':
            pass
        else:
            self.dates_yyyy_mm_list_updated, self.dates_yyyy_yyy_list_updated, self.msbas_np_names_updated = \
                (regression_do.filter_months_years_regression(self.years_1, self.months_1, self.years_2, self.months_2,
                                                              self.tsout_dates_yyyy_mm, self.tsout_dates_yyyy_yyy,
                                                              self.msbas_names,
                                                              self.regression_mode))
            self.numpy_array_paths = [(self.numpy_tif_arrays_path / x) for x in self.msbas_np_names_updated]
            self.numpy_array_paths.sort()
            numpy_paths_save = [str(x) for x in self.numpy_array_paths]
            dict_updated = {'mode': self.regression_mode, 'YYYYMMDDTHHMMSS': self.dates_yyyy_mm_list_updated,
                            'YYYY.YYY': self.dates_yyyy_yyy_list_updated,
                            'msbas_file_names': self.msbas_np_names_updated,
                            'numpy_paths': numpy_paths_save}
            self.msbasv3_json_save_name = (self.regression_data_used_path /
                                           ('MSBAS_TSOUT_' + self.regression_mode + '.json'))
            general_functions.write_json_from_dict(dict_updated, self.msbasv3_json_save_name)
        base_array = general_functions.open_numpy_file(self.numpy_array_paths[0])
        base_array_flat = base_array.flatten()
        self.tif_height = base_array.shape[0]
        self.tif_width = base_array.shape[1]
        for msbas_np_file in self.numpy_array_paths[1:]:
            np_arr = general_functions.open_numpy_file(msbas_np_file)
            np_arr_flat = np_arr.flatten()
            base_array_flat = np.vstack((base_array_flat, np_arr_flat))
        if self.regression_mode == 'all':
            dates_yyyy_yyy_float_array = np.array([float(x) for x in self.tsout_dates_yyyy_yyy])
            self.time_yyyy_yyy_use = dates_yyyy_yyy_float_array
        else:
            dates_yyyy_yyy_float_array = np.array([float(x) for x in self.dates_yyyy_yyy_list_updated])
            self.time_yyyy_yyy_use = dates_yyyy_yyy_float_array
        self.final_array_use = base_array_flat

    def parallel_polyfit(self):
        pool = multiprocessing.Pool()
        num_of_cols = self.final_array_use.shape[1]
        block_size = int(num_of_cols / self.num_of_cores)
        block_starts = [x for x in range(0, num_of_cols, block_size)]
        numpy_array_list = [self.final_array_use[:, x:(x + block_size)] for x in block_starts]
        time_use = self.time_yyyy_yyy_use - np.min(self.time_yyyy_yyy_use)
        model_slopes = ([pool.apply(regression_do.make_regressions_parallel,
                                    args=(np_array, time_use))
                         for np_array in numpy_array_list])
        pool.close()
        pool.join()
        base_np_arr = np.array([])
        for array in model_slopes:
            base_np_arr = np.append(base_np_arr, array.flatten())
        reshape_array = np.reshape(base_np_arr, (self.tif_height, self.tif_width)).astype('float32')
        name_1_s = self.path_1_name.split('_')
        name_1_u = 'P' + name_1_s[1] + '_F' + name_1_s[-1]
        name_2_s = self.path_2_name.split('_')
        name_2_u = 'P' + name_2_s[1] + '_F' + name_2_s[-1]
        name_use = name_1_u + '_' + name_2_u
        self.coh_save_name = (self.regression_tif_plots_path / ('coherence_mask_' + name_use + '.tif'))
        if not self.coh_save_name.exists():
            shutil.copy(self.coherence_tif_path, self.coh_save_name)
        self.regression_tif_save_name = (self.regression_tif_plots_path /
                                         ('linear_rate_' + name_use + '_' + self.regression_mode + '.tif'))
        rasterio_basic_functions.write_reprojected_array_2_tif(reshape_array,
                                                               self.regression_tif_save_name,
                                                               self.base_kwargs)

    def run_all(self):
        self.make_dirs_()
        self.make_numpy_arrays()
        self.read_tsout_and_parse()
        self.make_large_array_for_regression()
        self.parallel_polyfit()


if __name__ == '__main__':
    sys_index_var_1 = sys.argv[1]
    sys_index_var_2 = sys.argv[2]
    sys_index_var_3 = sys.argv[3]
    sys_index_var_4 = sys.argv[4]
    sys_index_var_5 = sys.argv[5]
    run_regressions = regression_and_time_series(sys_index_var_1, sys_index_var_2, sys_index_var_3, sys_index_var_4,
                                                 sys_index_var_5)
    run_regressions.run_all()
