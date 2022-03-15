import sys, os

sys.path.append("..")

from pathlib import Path
from max_curvature_codes import rename_msbas_log
from max_curvature_codes import make_norms_list_curvature
from max_curvature_codes import make_L_and_curve_plots
from general_set_up_functions import general_functions

import pandas as pd
import shutil


class max_curvature_extract():
    def __init__(self):
        self.base_path = Path.cwd().absolute().parent
        self.time_series_dir = self.base_path / ("MSBASv3_Timeseries")
        self.L_and_curv_plot_directory = self.time_series_dir / "L_curv_plots"
        self.norms_data_curv_data_directory = self.time_series_dir / "norm_curv_data"
        self.max_curv_msbas_path = self.time_series_dir / 'max_curvature_msbasv3'
        self.header_1_dirs = [Path(str(x)) for x in self.time_series_dir.glob('header_1*/MSBAS_LOG.txt')]
        self.header_2_dirs = [Path(str(x)) for x in self.time_series_dir.glob('header_2*/MSBAS_LOG.txt')]
        self.header_3_dirs = [Path(str(x)) for x in self.time_series_dir.glob('header_3*/MSBAS_LOG.txt')]
        self.header_1_dirs_renamed = None
        self.header_1_dirs_renamed_str = None
        self.header_2_dirs_renamed = None
        self.header_2_dirs_renamed_str = None
        self.header_3_dirs_renamed = None
        self.header_3_dirs_renamed_str = None
        self.list_factor = 0
        self.lam_list_1 = None
        self.curv_1 = None
        self.L_value_1 = None
        self.log_x_1 = None
        self.log_axy_1 = None

        self.lam_list_2 = None
        self.curv_2 = None
        self.L_value_2 = None
        self.log_x_2 = None
        self.log_axy_2 = None

        self.lam_list_3 = None
        self.curv_3 = None
        self.L_value_3 = None
        self.log_x_3 = None
        self.log_axy_3 = None

    def make_directories(self):
        if not self.L_and_curv_plot_directory.exists():
            self.L_and_curv_plot_directory.mkdir()
        if not self.norms_data_curv_data_directory.exists():
            self.norms_data_curv_data_directory.mkdir()
        if not self.max_curv_msbas_path.exists():
            self.max_curv_msbas_path.mkdir()

    def change_names(self):
        if self.header_1_dirs:
            rename_msbas_log.rename_msbas_log_f(self.header_1_dirs)
        if self.header_2_dirs:
            rename_msbas_log.rename_msbas_log_f(self.header_2_dirs)
        if self.header_3_dirs:
            rename_msbas_log.rename_msbas_log_f(self.header_3_dirs)

    def add_new_paths_lists(self):
        if self.header_1_dirs:
            self.header_1_dirs_renamed = [Path(str(x)) for x in self.time_series_dir.glob('header_1*/h*_MSBAS_LOG.txt')]
            self.header_1_dirs_renamed_str = ([str(x) for
                                               x in self.time_series_dir.glob('header_1*/')])
        if self.header_2_dirs:
            self.header_2_dirs_renamed = [Path(str(x)) for x in self.time_series_dir.glob('header_2*/h*_MSBAS_LOG.txt')]
            self.header_2_dirs_renamed_str = ([str(x) for
                                               x in self.time_series_dir.glob('header_2*/')])
        if self.header_3_dirs:
            self.header_3_dirs_renamed = [Path(str(x)) for x in self.time_series_dir.glob('header_3*/h*_MSBAS_LOG.txt')]
            self.header_3_dirs_renamed_str = ([str(x) for
                                               x in self.time_series_dir.glob('header_3*/')])

    def curvature_f(self):
        if self.header_1_dirs_renamed:
            l_value_1 = int(str(self.header_1_dirs_renamed[0].name).split('_')[1])
            norms_file_name_1 = self.norms_data_curv_data_directory / ('norms_axy_x_L_' + str(l_value_1) + '.json')
            curv_file_name_1 = self.norms_data_curv_data_directory / ('curv_L_' + str(l_value_1) + '.json')

            norms_file_name_1_txt = self.norms_data_curv_data_directory / ('norms_axy_x_L_' + str(l_value_1) + '.txt')
            curv_file_name_1_txt = self.norms_data_curv_data_directory / ('curv_L_' + str(l_value_1) + '.txt')

            log_array_norms_x_1, log_array_norms_axy_1, lam_lis_1, lam_axy_1 = \
                make_norms_list_curvature.make_norms_list_f(self.header_1_dirs_renamed, self.list_factor)
            curv_1 = (make_norms_list_curvature.find_curvatures
                      (log_array_norms_x_1, log_array_norms_axy_1, lam_lis_1, lam_axy_1))
            max_lam_1 = make_norms_list_curvature.find_max_curv_lambda(curv_1, lam_lis_1)
            max_lam_str_find_1 = str(max_lam_1).split('.')[1]
            matching_1 = [Path(s) for s in self.header_1_dirs_renamed_str if max_lam_str_find_1 in s][0]
            dst_1 = self.max_curv_msbas_path / matching_1.name
            if dst_1.exists():
                shutil.rmtree(dst_1)
                shutil.copytree(matching_1, dst_1)
            else:
                shutil.copytree(matching_1, dst_1)
            dict_base_norms_1 = {}
            dict_base_curv_1 = {}
            dict_base_norms_1['Lambda'] = list(lam_lis_1)
            dict_base_curv_1['Lambda'] = list(lam_lis_1)
            dict_base_norms_1['AXY_Norms'] = list(log_array_norms_axy_1)
            dict_base_norms_1['X_Norms'] = list(log_array_norms_x_1)
            dict_base_curv_1['curvature'] = list(curv_1)
            text_df_norms_1 = pd.DataFrame.from_dict(dict_base_norms_1)
            text_df_curv_1 = pd.DataFrame.from_dict(dict_base_curv_1)
            general_functions.write_json_from_dict(dict_base_norms_1, norms_file_name_1)
            general_functions.write_json_from_dict(dict_base_curv_1, curv_file_name_1)
            text_df_norms_1.to_csv(norms_file_name_1_txt, sep=' ', index=False)
            text_df_curv_1.to_csv(curv_file_name_1_txt, sep=' ', index=False)
            self.lam_list_1 = lam_lis_1
            self.curv_1 = curv_1
            self.L_value_1 = l_value_1
            self.log_x_1 = log_array_norms_x_1
            self.log_axy_1 = log_array_norms_axy_1
        if self.header_2_dirs_renamed:
            l_value_2 = int(str(self.header_2_dirs_renamed[0].name).split('_')[1])
            norms_file_name_2 = self.norms_data_curv_data_directory / ('norms_axy_x_L_' + str(l_value_2) + '.json')
            curv_file_name_2 = self.norms_data_curv_data_directory / ('curv_L_' + str(l_value_2) + '.json')
            norms_file_name_2_txt = self.norms_data_curv_data_directory / ('norms_axy_x_L_' + str(l_value_2) + '.txt')
            curv_file_name_2_txt = self.norms_data_curv_data_directory / ('curv_L_' + str(l_value_2) + '.txt')
            log_array_norms_x_2, log_array_norms_axy_2, lam_lis_2, lam_axy_2 = \
                make_norms_list_curvature.make_norms_list_f(self.header_2_dirs_renamed, self.list_factor)
            curv_2 = (make_norms_list_curvature.find_curvatures
                      (log_array_norms_x_2, log_array_norms_axy_2, lam_lis_2, lam_axy_2))
            max_lam_2 = make_norms_list_curvature.find_max_curv_lambda(curv_2, lam_lis_2)
            max_lam_str_find_2 = str(max_lam_2).split('.')[1]
            matching_2 = [Path(s) for s in self.header_2_dirs_renamed_str if max_lam_str_find_2 in s][0]
            matching_2_u = Path(self.time_series_dir/ matching_2.stem)
            dst_2 = self.max_curv_msbas_path / matching_2.stem
            if dst_2.exists():
                shutil.rmtree(dst_2)
                shutil.copytree(matching_2_u, dst_2)
            else:
                shutil.copytree(matching_2_u, dst_2)
            dict_base_norms_2 = {}
            dict_base_curv_2 = {}
            dict_base_norms_2['Lambda'] = list(lam_lis_2)
            dict_base_curv_2['Lambda'] = list(lam_lis_2)
            dict_base_norms_2['AXY_Norms'] = list(log_array_norms_axy_2)
            dict_base_norms_2['X_Norms'] = list(log_array_norms_x_2)
            dict_base_curv_2['curvature'] = list(curv_2)
            text_df_norms_2 = pd.DataFrame.from_dict(dict_base_norms_2)
            text_df_curv_2 = pd.DataFrame.from_dict(dict_base_curv_2)
            general_functions.write_json_from_dict(dict_base_norms_2, norms_file_name_2)
            general_functions.write_json_from_dict(dict_base_curv_2, curv_file_name_2)
            text_df_norms_2.to_csv(norms_file_name_2_txt, sep=' ', index=False)
            text_df_curv_2.to_csv(curv_file_name_2_txt, sep=' ', index=False)
            self.lam_list_2 = lam_lis_2
            self.curv_2 = curv_2
            self.L_value_2 = l_value_2
            self.log_x_2 = log_array_norms_x_2
            self.log_axy_2 = log_array_norms_axy_2
        if self.header_3_dirs_renamed:
            l_value_3 = int(str(self.header_3_dirs_renamed[0].name).split('_')[1])
            norms_file_name_3 = self.norms_data_curv_data_directory / ('norms_axy_x_L_' + str(l_value_3) + '.json')
            curv_file_name_3 = self.norms_data_curv_data_directory / ('curv_L_' + str(l_value_3) + '.json')

            norms_file_name_3_txt = self.norms_data_curv_data_directory / ('norms_axy_x_L_' + str(l_value_3) + '.txt')
            curv_file_name_3_txt = self.norms_data_curv_data_directory / ('curv_L_' + str(l_value_3) + '.txt')

            log_array_norms_x_3, log_array_norms_axy_3, lam_lis_3, lam_axy_3 = \
                make_norms_list_curvature.make_norms_list_f(self.header_3_dirs_renamed, self.list_factor)
            curv_3 = (make_norms_list_curvature.find_curvatures
                      (log_array_norms_x_3, log_array_norms_axy_3, lam_lis_3, lam_axy_3))
            max_lam_3 = make_norms_list_curvature.find_max_curv_lambda(curv_3, lam_lis_3)
            max_lam_str_find_3 = str(max_lam_3).split('.')[1]
            matching_3 = [Path(s) for s in self.header_3_dirs_renamed_str if max_lam_str_find_3 in s][0]
            dst_3 = self.max_curv_msbas_path / matching_3.name
            if dst_3.exists():
                shutil.rmtree(dst_3)
                shutil.copytree(matching_3, dst_3)
            else:
                shutil.copytree(matching_3, dst_3)
            dict_base_norms_3 = {}
            dict_base_curv_3 = {}
            dict_base_norms_3['Lambda'] = list(lam_lis_3)
            dict_base_curv_3['Lambda'] = list(lam_lis_3)
            dict_base_norms_3['AXY_Norms'] = list(log_array_norms_axy_3)
            dict_base_norms_3['X_Norms'] = list(log_array_norms_x_3)
            dict_base_curv_3['curvature'] = list(curv_3)
            text_df_norms_3 = pd.DataFrame.from_dict(dict_base_norms_3)
            text_df_curv_3 = pd.DataFrame.from_dict(dict_base_curv_3)
            general_functions.write_json_from_dict(dict_base_norms_3, norms_file_name_3)
            general_functions.write_json_from_dict(dict_base_curv_3, curv_file_name_3)
            text_df_norms_3.to_csv(norms_file_name_3_txt, sep=' ', index=False)
            text_df_curv_3.to_csv(curv_file_name_3_txt, sep=' ', index=False)
            self.lam_list_3 = lam_lis_3
            self.curv_3 = curv_3
            self.L_value_3 = l_value_3
            self.log_x_3 = log_array_norms_x_3
            self.log_axy_3 = log_array_norms_axy_3

    def run_all(self):
        self.make_directories()
        self.change_names()
        self.add_new_paths_lists()
        self.curvature_f()
        if (self.lam_list_1 and self.L_value_1) is not None:
            make_L_and_curve_plots.make_max_curv_plot(self.lam_list_1, self.curv_1, self.L_and_curv_plot_directory,
                                                      self.L_value_1)
        if (self.lam_list_2 and self.L_value_2) is not None:
            make_L_and_curve_plots.make_max_curv_plot(self.lam_list_2, self.curv_2, self.L_and_curv_plot_directory,
                                                      self.L_value_2)
        if (self.lam_list_3 and self.L_value_3) is not None:
            make_L_and_curve_plots.make_max_curv_plot(self.lam_list_3, self.curv_3, self.L_and_curv_plot_directory,
                                                      self.L_value_3)

        make_L_and_curve_plots.make_L_curv_plot(self.log_x_1, self.log_x_2, self.log_x_3, self.log_axy_1,
                                                self.log_axy_2, self.log_axy_3, self.L_and_curv_plot_directory)


if __name__ == '__main__':
    run_max_curv_object = max_curvature_extract()
    run_max_curv_object.run_all()
