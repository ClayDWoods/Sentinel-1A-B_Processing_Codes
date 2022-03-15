import os
import shutil
import sys
from pathlib import Path


class gps_reference_set_up():
    def __init__(self, path_1, path_2):
        self.base_path = Path.cwd().absolute().parent
        self.combine_pairs_path = self.base_path / 'Combined_Pairs/Path_114_12_referenced'
        self.path_pairs_list_months_path_1 = ([Path(str(x)) for x in
                                               self.base_path.glob('Pair_P' + path_1
                                                                   + '_F*/Plots_Regression/Linear_Rate_Spatial_Plots/'
                                                                   + '*all.tif')])
        self.path_pairs_list_months_path_2 = ([Path(str(x)) for x in
                                               self.base_path.glob('Pair_P' + path_2
                                                                   + '_F*/Plots_Regression/Linear_Rate_Spatial_Plots/'
                                                                   + '*all.tif')])
        self.path_pairs_list_coh_path_1 = ([Path(str(x)) for x in
                                            self.base_path.glob('Pair_P' + path_1
                                                                + '_F*/Plots_Regression/Linear_Rate_Spatial_Plots/coh'
                                                                  '*.tif')])
        self.path_pairs_list_coh_path_2 = ([Path(str(x)) for x in
                                            self.base_path.glob('Pair_P' + path_2
                                                                + '_F*/Plots_Regression/Linear_Rate_Spatial_Plots/coh'
                                                                  '*.tif')])
        self.path_pairs_list = (self.path_pairs_list_months_path_1 +
                                self.path_pairs_list_months_path_2 +
                                self.path_pairs_list_coh_path_1 +
                                self.path_pairs_list_coh_path_2)

    def copy_files_2_loc(self):
        for file in self.path_pairs_list:
            shutil.copy(file, self.combine_pairs_path)

    def run_set_up(self):
        self.copy_files_2_loc()


if __name__ == '__main__':
    sys_index_var_1 = sys.argv[1]
    sys_index_var_2 = sys.argv[2]
    run_gps_reference_set_up = gps_reference_set_up(sys_index_var_1, sys_index_var_2)
    run_gps_reference_set_up.run_set_up()
