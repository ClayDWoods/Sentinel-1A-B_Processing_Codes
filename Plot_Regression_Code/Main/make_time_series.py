import sys
sys.path.append("..")

from Main.general_functions_all import general_functions

import numpy as np
import matplotlib.pyplot as plt

plt.style.use('seaborn-whitegrid')
from pathlib import Path


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


class make_time_series_regressions():
    def __init__(self, time_series_save_name, half_size, lon_list, lat_list):
        self.base_path = Path.cwd().absolute().parent
        self.kwargs_transform_path = self.base_path / 'tif_meta_data/kwargs_tif.json'
        self.binary_spatial_array_path = self.base_path / 'Binary_Spatial_Arrays'
        self.tsout_path_months = self.base_path / 'Data_used_regression/MSBAS_TSOUT_months.json'
        self.tsout_path_all = self.base_path / 'Data_used_regression/MSBAS_TSOUT_all.json'
        self.time_series_path = self.base_path / 'Time_Series_Plots'
        self.time_series_save_name = self.time_series_path / ('Time_Series_' + time_series_save_name + '.png')
        self.to_mm = 1000
        self.lon_list = lon_list
        self.lat_list = lat_list
        self.half_size = half_size
        kwargs_open = general_functions.open_json_file(self.kwargs_transform_path)
        self.transform_list = kwargs_open['transform']
        tsout_open_month = general_functions.open_json_file(self.tsout_path_months)
        tsout_open_all = general_functions.open_json_file(self.tsout_path_all)
        self.times_month = np.array([float(x) for x in tsout_open_month['YYYY.YYY']])
        self.times_month_shifted = self.times_month - np.min(self.times_month)
        self.times_all = np.array([float(x) for x in tsout_open_all['YYYY.YYY']])
        self.times_all_shifted = self.times_all - np.min(self.times_all)
        msbas_used_np_month = tsout_open_month['msbas_file_names']
        msbas_used_np_all = tsout_open_all['msbas_file_names']
        test_nps_months = [Path(str(self.binary_spatial_array_path) + '/' + x) for x in msbas_used_np_month]
        test_nps_months.sort()
        test_nps_all = [Path(str(self.binary_spatial_array_path) + '/' + x) for x in msbas_used_np_all]
        test_nps_all.sort()
        self.binary_spatial_arrays_months = test_nps_months
        self.binary_spatial_arrays_all = test_nps_all
        self.def_months = None
        self.coeff_months = None
        self.def_all = None
        self.coeff_all = None
        self.ylim_min = None
        self.ylim_max = None

    def run_mean_lists(self):
        lon_cols, lat_rows = general_functions.lon_lat_2_col_row_using_transform(self.transform_list,
                                                                                 self.lon_list,
                                                                                 self.lat_list)
        def_months, coeff_months = mean_sliced_list(self.binary_spatial_arrays_months,
                                                    self.half_size,
                                                    lon_cols,
                                                    lat_rows,
                                                    self.times_month_shifted,
                                                    self.to_mm)
        def_all, coeff_all = mean_sliced_list(self.binary_spatial_arrays_all,
                                              self.half_size,
                                              lon_cols,
                                              lat_rows,
                                              self.times_all_shifted,
                                              self.to_mm)
        self.def_months = def_months
        self.coeff_months = coeff_months
        self.def_all = def_all
        self.coeff_all = coeff_all

    def find_ylims(self):
        base_a = np.array([])
        for def_r in self.def_all:
            np_as = np.array(def_r)
            base_a = np.append(base_a, np_as)
        base = 5
        ylim_min = general_functions.nearest_bound_coord(np.min(base_a), base, 'down')
        ylim_max = general_functions.nearest_bound_coord(np.max(base_a), base, 'up')
        self.ylim_min = ylim_min
        self.ylim_max = ylim_max

    def make_time_series_plots(self):
        base = 3
        if len(self.lon_list) == 1:
            n_rows = 1
            n_cols = 1
            fig, axs1 = plt.subplots(figsize=(12, 10), nrows=n_rows, ncols=n_cols, facecolor='w', edgecolor='k')
        else:
            num_plots = general_functions.nearest_bound_coord(len(self.lon_list), base, 'up')
            n_rows = int(num_plots / base)
            col_size = 20 / (base - n_rows + 1)
            fig, axs1 = plt.subplots(figsize=(24, col_size), nrows=n_rows, ncols=3, facecolor='w', edgecolor='k')
            axs1 = axs1.ravel()
        for i in range(len(self.lon_list)):
            if len(self.lon_list) == 1:
                axs1.set_ylim([self.ylim_min, self.ylim_max])
                axs1.scatter(self.times_all, self.def_all, color='r')
                poly = np.poly1d(self.coeff_months[0])
                rate = round(self.coeff_months[0][0], 2)
                new_y1 = poly(self.times_month_shifted)
                axs1.plot(self.times_month, new_y1, color='b')
                axs1.set_ylabel('mm', fontsize=16)
                axs1.set_xticks([2017, 2017.5, 2018, 2018.5, 2019, 2019.5, 2020, 2020.5, 2021, 2021.5])
                axs1.text(0.01, 0.97, ('linear_rate=' + str(rate) + ' mm/yr'), fontsize=12, transform=axs1.transAxes)
                axs1.text(0.01, 0.93, ('lon=' + str(self.lon_list[0]) + ' ' + 'lat=' + str(self.lat_list[0])),
                          fontsize=12,
                          transform=axs1.transAxes)
            else:
                axs1[i].set_ylim([self.ylim_min, self.ylim_max])
                axs1[i].scatter(self.times_all, self.def_all[i], color='r')
                poly = np.poly1d(self.coeff_months[i])
                rate = round(self.coeff_months[i][0], 2)
                new_y1 = poly(self.times_month_shifted)
                axs1[i].plot(self.times_month, new_y1, color='b')
                axs1[i].set_ylabel('mm', fontsize=16)
                axs1[i].set_xticks([2017, 2017.5, 2018, 2018.5, 2019, 2019.5, 2020, 2020.5, 2021, 2021.5])
                axs1[i].text(0.01, 0.97, ('linear_rate=' + str(rate) + ' mm/yr'), fontsize=12,
                             transform=axs1[i].transAxes)
                axs1[i].text(0.01, 0.93, ('lon=' + str(self.lon_list[i]) + ' ' + 'lat=' + str(self.lat_list[i])),
                             fontsize=12,
                             transform=axs1[i].transAxes)
        if len(self.lon_list) == 1:
            pass
        else:
            num_plots = general_functions.nearest_bound_coord(len(self.lon_list), base, 'up')
            num_plots_unused = num_plots - len(self.lon_list)
            for j in range(len(self.lon_list), len(self.lon_list) + num_plots_unused):
                fig.delaxes(axs1[j])
        fig.suptitle('Time series plots', y=0.99, fontsize=16, fontweight="bold")
        fig.tight_layout()
        if not self.time_series_save_name.exists():
            print('output time series png to:')
            print(self.time_series_save_name)
            fig.savefig(self.time_series_save_name, dpi=100, facecolor='w', edgecolor='k',
                        orientation='portrait', bbox_inches='tight', pad_inches=0.3)
        else:
            print("file exists")
            print("Change name of save file and re-run")
        plt.show()

    def run_all(self):
        self.run_mean_lists()
        self.find_ylims()
        self.make_time_series_plots()
