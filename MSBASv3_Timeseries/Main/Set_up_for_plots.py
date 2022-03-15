import os
import shutil
import sys
from pathlib import Path


class plotting_set_up():
    def __init__(self, header_dir_name_use, path_1_dir, path_2_dir):
        self.base_path = Path.cwd().absolute().parent
        self.base_plots_path = self.base_path / 'Plots_Regression'
        self.msbasv3_base_path = self.base_path / 'MSBASv3_Timeseries'
        self.coherence_mask_base_path = self.base_path / 'coherence_mask'
        self.coherence_mask_new_path = self.base_plots_path / 'coherence_mask'
        self.stack_path = self.msbasv3_base_path / 'Stack'
        self.stack_path_new = self.base_plots_path / 'Stack'
        self.reference_points_path = self.msbasv3_base_path / 'Reference_point'
        self.reference_points_path_new = self.base_plots_path / 'Reference_point'
        self.kml_path_1 = self.base_path / (path_1_dir + '/' + 'Main/Data/kml_frame')
        self.kml_path_2 = self.base_path / (path_2_dir + '/' + 'Main/Data/kml_frame')
        self.kml_path_1_save = self.base_plots_path / (path_1_dir + '_kml_coords')
        self.kml_path_2_save = self.base_plots_path / (path_2_dir + '_kml_coords')
        self.param_file_1 = self.base_path / (path_1_dir + '/' + 'Main/Parameter_File')
        self.param_file_2 = self.base_path / (path_2_dir + '/' + 'Main/Parameter_File')
        self.param_file_1_save = self.base_plots_path / (path_1_dir + '_Parameter_File')
        self.param_file_2_save = self.base_plots_path / (path_2_dir + '_Parameter_File')
        self.snaphu_logs_1 = self.base_path / (path_1_dir + '/' + 'Main/Snaphu_Logs')
        self.snaphu_logs_2 = self.base_path / (path_2_dir + '/' + 'Main/Snaphu_Logs')
        self.snaphu_logs_1_save = self.base_plots_path / (path_1_dir + '_Snaphu_Logs')
        self.snaphu_logs_2_save = self.base_plots_path / (path_2_dir + '_Snaphu_Logs')
        self.plat_params_1 = self.base_path / (path_1_dir + '/' + 'Main/Time_Series')
        self.plat_params_2 = self.base_path / (path_2_dir + '/' + 'Main/Time_Series')
        self.plat_params_1_save = self.base_plots_path / (path_1_dir + '_Time_Series')
        self.plat_params_2_save = self.base_plots_path / (path_2_dir + '_Time_Series')
        self.header_file_path = self.msbasv3_base_path / header_dir_name_use
        self.header_file_path_new = self.base_plots_path / header_dir_name_use

    def make_dirs(self):
        if not self.base_plots_path.exists():
            os.mkdir(self.base_plots_path)

    def copy_files_and_dirs(self):
        if self.header_file_path_new.exists():
            shutil.rmtree(self.header_file_path_new)
        if self.coherence_mask_new_path.exists():
            shutil.rmtree(self.coherence_mask_new_path)
        if self.stack_path_new.exists():
            shutil.rmtree(self.stack_path_new)
        if self.reference_points_path_new.exists():
            shutil.rmtree(self.reference_points_path_new)
        if self.kml_path_1_save.exists():
            shutil.rmtree(self.kml_path_1_save)
        if self.kml_path_2_save.exists():
            shutil.rmtree(self.kml_path_2_save)
        if self.param_file_1_save.exists():
            shutil.rmtree(self.param_file_1_save)
        if self.param_file_2_save.exists():
            shutil.rmtree(self.param_file_2_save)
        if self.snaphu_logs_1_save.exists():
            shutil.rmtree(self.snaphu_logs_1_save)
        if self.snaphu_logs_2_save.exists():
            shutil.rmtree(self.snaphu_logs_2_save)
        if self.plat_params_1_save.exists():
            shutil.rmtree(self.plat_params_1_save)
        if self.plat_params_2_save.exists():
            shutil.rmtree(self.plat_params_2_save)
        shutil.copytree(self.header_file_path, self.header_file_path_new)
        shutil.copytree(self.coherence_mask_base_path, self.coherence_mask_new_path)
        shutil.copytree(self.stack_path, self.stack_path_new)
        shutil.copytree(self.reference_points_path, self.reference_points_path_new)
        shutil.copytree(self.kml_path_1, self.kml_path_1_save)
        shutil.copytree(self.kml_path_2, self.kml_path_2_save)

        shutil.copytree(self.param_file_1, self.param_file_1_save)
        shutil.copytree(self.param_file_2, self.param_file_2_save)

        shutil.copytree(self.snaphu_logs_1, self.snaphu_logs_1_save)
        shutil.copytree(self.snaphu_logs_2, self.snaphu_logs_2_save)

        shutil.copytree(self.plat_params_1, self.plat_params_1_save)
        shutil.copytree(self.plat_params_2, self.plat_params_2_save)

    def run_all_f(self):
        self.make_dirs()
        self.copy_files_and_dirs()

if __name__ == '__main__':
    sys_index_var_1 = sys.argv[1]
    sys_index_var_2 = sys.argv[2]
    sys_index_var_3 = sys.argv[3]
    run_plot_set_up = plotting_set_up(sys_index_var_1, sys_index_var_2, sys_index_var_3)
    run_plot_set_up.run_all_f()
