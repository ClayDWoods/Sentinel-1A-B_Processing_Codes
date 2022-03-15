import os
import shutil
import sys
from pathlib import Path
from general_set_up_functions import rasterio_basic_functions


class gps_reference():
    def __init__(self, path_frame, gps_value):
        self.base_gps_value = float(gps_value)*0.001
        self.base_path = Path.cwd().absolute().parent
        self.combine_pairs_path = self.base_path / 'Combined_Pairs/Path_114_12_referenced'
        self.referenced_pairs_save_path = self.combine_pairs_path / 'referenced_path_frame'
        self.overlap_paths_save = self.combine_pairs_path / 'overlap_tifs'
        self.path_pairs_list_months = ([Path(str(x)) for x in
                                        self.combine_pairs_path.glob('lin*' + path_frame + '*all.tif')])[0]
        self.path_pair_referenced_name = (self.path_pairs_list_months.parent /
                                          (self.path_pairs_list_months.stem + '_ref_gps.tif'))

    def set_up_directories(self):
        if not self.referenced_pairs_save_path.exists():
            os.mkdir(self.referenced_pairs_save_path)
        if not self.overlap_paths_save.exists():
            os.mkdir(self.overlap_paths_save)

    def do_gps_reference(self):
        base_array = rasterio_basic_functions.tif_2_array(self.path_pairs_list_months)
        base_array[base_array != 0] += self.base_gps_value
        base_kwargs = rasterio_basic_functions.get_kwargs(self.path_pairs_list_months)
        rasterio_basic_functions.write_reprojected_array_2_tif(base_array, self.path_pair_referenced_name,
                                                               base_kwargs)
        shutil.copy(self.path_pair_referenced_name, self.referenced_pairs_save_path)

    def run_gps_reference(self):
        self.set_up_directories()
        self.do_gps_reference()


if __name__ == '__main__':
    sys_index_var_1 = sys.argv[1]
    sys_index_var_2 = sys.argv[2]
    gps_reference = gps_reference(sys_index_var_1, sys_index_var_2)
    gps_reference.run_gps_reference()
