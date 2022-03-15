import sys, os

sys.path.append("..")

from pathlib import Path
from make_asc_dsc_header_txt import loop_displacement_files_mcoh, header_txt_file_creation, \
    auto_find_reference_point
from resample_and_mask import Resample_2_same_grid
from general_set_up_functions import general_functions, rasterio_basic_functions
from Stack_Tools import stack_processing
from gacos_functions import make_incidence_tif as mk_in
from gacos_functions import apply_ztd as apz
from make_asc_dsc_header_txt import find_max_coh_lon_lat
import numpy as np


class msbasv3_time_series_set_up():
    def __init__(self, asc_dir_name, dsc_dir_name, path_1_full, path_2_full,
                 ref_type='gps', run_resample='true', ref_points=1, run_dtzd='false',
                 gacos_full='gacos_mixed', std_E='n_std'):
        self.gacos_or_non_gacos_full = gacos_full
        self.std_E = std_E
        self.ref_type = ref_type
        self.run_dtzd = run_dtzd
        self.run_resample = run_resample
        self.base_path = Path.cwd().absolute().parent
        self.base_path_all_dirs = Path.cwd().absolute().parent.parent
        self.coh_path_1 = None
        self.coh_path_2 = None
        self.path_1_full_name = path_1_full
        self.path_2_full_name = path_2_full
        self.time_series_dir = self.base_path / ("MSBASv3_Timeseries")
        self.visual_inspection_gacos_asc = self.time_series_dir / 'visual_remove_gacos_asc.json'
        self.visual_inspection_gacos_dsc = self.time_series_dir / 'visual_remove_gacos_dsc.json'
        self.gacos_or_non_gacos_asc = self.time_series_dir / 'gacos_or_non_gacos_use_asc.json'
        self.gacos_or_non_gacos_dsc = self.time_series_dir / 'gacos_or_non_gacos_use_dsc.json'
        self.time_series_parameters_dir = self.base_path / ("MSBASv3_Timeseries/Reference_point")
        self.time_series_parameters_save_file = (self.time_series_parameters_dir / ('reference_point_' +
                                                                                    path_2_full + '.json'))
        self.stack_path = self.base_path / ("MSBASv3_Timeseries/Stack")
        self.coherence_mask_path = self.base_path / ('coherence_mask/coherence_mask.tif')
        self.asc_pair_path = self.base_path / asc_dir_name
        self.dsc_pair_path = self.base_path / dsc_dir_name
        self.asc_snaphu_log_path = self.asc_pair_path / ('Main/Snaphu_Logs')
        self.dsc_snaphu_log_path = self.dsc_pair_path / ('Main/Snaphu_Logs')
        self.meta_path_asc = self.asc_pair_path / ('Main/Data/Meta_Data')
        self.meta_path_dsc = self.dsc_pair_path / ('Main/Data/Meta_Data')
        self.path_2_asc_ztd = self.base_path / 'MSBASv3_Timeseries/asc_ztd'
        self.path_2_dsc_ztd = self.base_path / 'MSBASv3_Timeseries/dsc_ztd'
        self.path_2_inc_asc_ztd = self.path_2_asc_ztd / 'incidence_rcut_mcoh_asc.tif'
        self.path_2_inc_dsc_ztd = self.path_2_dsc_ztd / 'incidence_rcut_mcoh_dsc.tif'
        self.asc_time_series_parameters = self.asc_pair_path / ('Main/Time_Series/platform_parameters.json')
        self.dsc_time_series_parameters = self.dsc_pair_path / ('Main/Time_Series/platform_parameters.json')
        self.asc_displacement_tif_pairs_list = (self.asc_pair_path.
                                                glob('Main/Processed_Data/20*20*/*displacement*V?.tif'))
        self.dsc_displacement_tif_pairs_list = (self.dsc_pair_path.
                                                glob('Main/Processed_Data/20*20*/*displacement*V?.tif'))
        self.asc_displacement_tif_pairs_list_ = (self.asc_pair_path.
                                                 glob('Main/Processed_Data/20*20*/*displacement*V?.tif'))
        self.dsc_displacement_tif_pairs_list_ = (self.dsc_pair_path.
                                                 glob('Main/Processed_Data/20*20*/*displacement*V?.tif'))
        self.asc_coh_tif_pairs_list = self.asc_pair_path.glob('Main/Processed_Data/20*20*/*coherence*V?.tif')
        self.dsc_coh_tif_pairs_list = self.dsc_pair_path.glob('Main/Processed_Data/20*20*/*coherence*V?.tif')
        self.asc_coh_tif_pairs_list_ = self.asc_pair_path.glob('Main/Processed_Data/20*20*/*coherence*V?.tif')
        self.dsc_coh_tif_pairs_list_ = self.dsc_pair_path.glob('Main/Processed_Data/20*20*/*coherence*V?.tif')
        self.asc_dsc_displacement_tif_combined_list = (list(self.asc_displacement_tif_pairs_list_)
                                                       + list(self.dsc_displacement_tif_pairs_list_))
        self.asc_dsc_coherence_tif_combined_list = (list(self.asc_coh_tif_pairs_list_) +
                                                    list(self.dsc_coh_tif_pairs_list_))
        self.asc_flag = 'asc'
        self.dsc_flag = 'dsc'
        self.parsed_list_mcoh_asc_name = self.time_series_dir / (self.asc_flag + '_file_list_msbasv3.json')
        self.parsed_list_mcoh_dsc_name = self.time_series_dir / (self.dsc_flag + '_file_list_msbasv3.json')
        self.c_size = 4
        if self.ref_type == 'gps':
            self.gps_lon = 90.401264
            self.gps_lat = 23.726679
        self.coherence_threshold = 0.15
        # self.lambda__ = np.linspace(0.01, 0.6, 25)
        self.lambda__ = np.array([0.01, 0.10833333, 0.20666667, 0.30500000, 0.40333333, 0.50166667, 0.60000000])
        self.L = [2]
        self.num_of_stds = 2
        self.num_of_ref = int(ref_points)
        self.half_size = 8
        self.coh_min = 0.5
        self.num_of_coh_blocks = 100
        self.num_of_max_to_check = 15
        self.num_cores = 16
        self.max_coh = None
        self.c_flag_col = None
        self.c_flag_row = None
        self.master_height = None
        self.master_width = None
        self.asc_heading = None
        self.asc_incidence = None
        self.asc_time = None
        self.dsc_heading = None
        self.dsc_incidence = None
        self.dsc_time = None
        self.asc_displacement_tif_pairs_list_cut = None
        self.dsc_displacement_tif_pairs_list_cut = None
        self.asc_displacement_tif_pairs_list_cut_ = None
        self.dsc_displacement_tif_pairs_list_cut_ = None
        self.asc_coh_tif_pairs_list_cut = None
        self.dsc_coh_tif_pairs_list_cut = None
        self.asc_coh_tif_pairs_list_cut_ = None
        self.dsc_coh_tif_pairs_list_cut_ = None
        self.asc_dsc_displacement_tif_combined_list_cut = None
        self.asc_dsc_coherence_tif_combined_list_cut = None
        self.asc_coh_tif_pairs_list_mcoh_gacos = None
        self.dsc_coh_tif_pairs_list_mcoh_gacos = None
        self.asc_coh_tif_pairs_list_mcoh_ = None
        self.dsc_coh_tif_pairs_list_mcoh_ = None
        self.asc_mcoh_4_ztd_base = None
        self.dsc_mcoh_4_ztd_base = None
        self.asc_rcut_4_ztd_base = None
        self.dsc_rcut_4_ztd_base = None
        self.asc_rcut_mcoh_4_ztd = None
        self.dsc_rcut_mcoh_4_ztd = None
        self.asc_mcoh_all = None
        self.dsc_mcoh_all = None
        self.asc_mcoh_gacos_all = None
        self.dsc_mcoh_gacos_all = None

    def read_asc_dsc_parameters(self):
        asc_params_dict = general_functions.open_json_file(self.asc_time_series_parameters)
        dsc_params_dict = general_functions.open_json_file(self.dsc_time_series_parameters)
        self.asc_heading = asc_params_dict['platformHeading']
        self.asc_incidence = asc_params_dict['incidenceAngleMidSwath']
        self.asc_time = asc_params_dict['startTime']
        self.dsc_heading = dsc_params_dict['platformHeading']
        self.dsc_incidence = dsc_params_dict['incidenceAngleMidSwath']
        self.dsc_time = dsc_params_dict['startTime']

    def add_new_tif_names(self):
        if self.visual_inspection_gacos_asc.exists() and self.visual_inspection_gacos_dsc.exists():
            visual_inspection_asc_dict = general_functions.open_json_file(self.visual_inspection_gacos_asc)
            visual_inspection_dsc_dict = general_functions.open_json_file(self.visual_inspection_gacos_dsc)
            visual_inspection_asc_list = [Path(x) for x in visual_inspection_asc_dict['trash_files']]
            visual_inspection_dsc_list = [Path(x) for x in visual_inspection_dsc_dict['trash_files']]
            visual_inspection_asc_list_converted = ([(self.base_path / (
                    x.absolute().parent.parent.parent.parent.stem + '/Main/Processed_Data/' +
                    x.absolute().parent.stem + '/' + x.absolute().parent.stem + '_displacement_VV_rcut.tif'))
                                                     for x in visual_inspection_asc_list])
            visual_inspection_dsc_list_converted = ([(self.base_path / (
                    x.absolute().parent.parent.parent.parent.stem + '/Main/Processed_Data/' +
                    x.absolute().parent.stem + '/' + x.absolute().parent.stem + '_displacement_VV_rcut.tif'))
                                                     for x in visual_inspection_dsc_list])
            visual_inspection_asc_list_converted_coh = ([(self.base_path / (
                    x.absolute().parent.parent.parent.parent.stem + '/Main/Processed_Data/' +
                    x.absolute().parent.stem + '/' + x.absolute().parent.stem + '_coherence_VV_rcut.tif'))
                                                         for x in visual_inspection_dsc_list])
            visual_inspection_dsc_list_converted_coh = ([(self.base_path / (
                    x.absolute().parent.parent.parent.parent.stem + '/Main/Processed_Data/' +
                    x.absolute().parent.stem + '/' + x.absolute().parent.stem + '_coherence_VV_rcut.tif'))
                                                         for x in visual_inspection_dsc_list])

            ######################
            asc_displacement_tif_pairs_list_cut_1 = [Path(str(x)) for x in self.asc_pair_path.glob(
                'Main/Processed_Data/20*20*/*displacement*rcut.tif')]
            self.asc_displacement_tif_pairs_list_cut = ([x for x in asc_displacement_tif_pairs_list_cut_1 if x not
                                                         in visual_inspection_asc_list_converted])
            self.asc_displacement_tif_pairs_list_cut.sort()
            ######################
            ######################
            dsc_displacement_tif_pairs_list_cut_1 = [Path(str(x)) for x in self.dsc_pair_path.glob(
                'Main/Processed_Data/20*20*/*displacement*rcut.tif')]
            self.dsc_displacement_tif_pairs_list_cut = ([x for x in dsc_displacement_tif_pairs_list_cut_1 if x not
                                                         in visual_inspection_dsc_list_converted])
            self.dsc_displacement_tif_pairs_list_cut.sort()
            ######################

            ######################
            asc_coh_tif_pairs_list_cut_1 = [Path(str(x)) for x in self.asc_pair_path.glob(
                'Main/Processed_Data/20*20*/*coherence*rcut.tif')]
            self.asc_coh_tif_pairs_list_cut = ([x for x in asc_coh_tif_pairs_list_cut_1 if x not
                                                in visual_inspection_asc_list_converted_coh])
            self.asc_coh_tif_pairs_list_cut.sort()

            dsc_coh_tif_pairs_list_cut_1 = [Path(str(x)) for x in self.dsc_pair_path.glob(
                'Main/Processed_Data/20*20*/*coherence*rcut.tif')]
            self.dsc_coh_tif_pairs_list_cut = ([x for x in dsc_coh_tif_pairs_list_cut_1 if x not
                                                in visual_inspection_dsc_list_converted_coh])
            self.dsc_coh_tif_pairs_list_cut.sort()
            ######################

            self.asc_dsc_displacement_tif_combined_list_cut = (self.asc_displacement_tif_pairs_list_cut +
                                                               self.dsc_displacement_tif_pairs_list_cut)
            self.asc_dsc_coherence_tif_combined_list_cut = (self.asc_coh_tif_pairs_list_cut +
                                                            self.dsc_coh_tif_pairs_list_cut)
        else:
            self.asc_displacement_tif_pairs_list_cut = [Path(str(x)) for x in self.asc_pair_path.glob(
                'Main/Processed_Data/20*20*/*displacement*rcut.tif')]
            self.asc_displacement_tif_pairs_list_cut.sort()
            self.dsc_displacement_tif_pairs_list_cut = [Path(str(x)) for x in self.dsc_pair_path.glob(
                'Main/Processed_Data/20*20*/*displacement*rcut.tif')]
            self.dsc_displacement_tif_pairs_list_cut.sort()

            self.asc_coh_tif_pairs_list_cut = [Path(str(x)) for x in self.asc_pair_path.glob(
                'Main/Processed_Data/20*20*/*coherence*rcut.tif')]
            self.asc_coh_tif_pairs_list_cut.sort()
            self.dsc_coh_tif_pairs_list_cut = [Path(str(x)) for x in self.dsc_pair_path.glob(
                'Main/Processed_Data/20*20*/*coherence*rcut.tif')]
            self.dsc_coh_tif_pairs_list_cut.sort()
            self.asc_dsc_displacement_tif_combined_list_cut = (self.asc_displacement_tif_pairs_list_cut +
                                                               self.dsc_displacement_tif_pairs_list_cut)
            self.asc_dsc_coherence_tif_combined_list_cut = (self.asc_coh_tif_pairs_list_cut +
                                                            self.dsc_coh_tif_pairs_list_cut)

        asc_mcoh_4_ztd_base = [Path(str(x)) for x in self.asc_pair_path.glob(
            'Main/Processed_Data/20*20*/*displacement*rcut_mcoh.tif')]
        asc_mcoh_4_ztd_base.sort()
        self.asc_mcoh_4_ztd_base = asc_mcoh_4_ztd_base[0]

        dsc_mcoh_4_ztd_base = [Path(str(x)) for x in self.dsc_pair_path.glob(
            'Main/Processed_Data/20*20*/*displacement*rcut_mcoh.tif')]
        dsc_mcoh_4_ztd_base.sort()
        self.dsc_mcoh_4_ztd_base = dsc_mcoh_4_ztd_base[0]

        self.asc_rcut_4_ztd_base = self.asc_displacement_tif_pairs_list_cut[0]
        self.dsc_rcut_4_ztd_base = self.dsc_displacement_tif_pairs_list_cut[0]

        self.asc_rcut_mcoh_4_ztd = [Path(str(x)) for x in self.asc_pair_path.glob(
            'Main/Processed_Data/20*20*/*displacement*rcut_mcoh.tif')]
        self.asc_rcut_mcoh_4_ztd.sort()

        self.dsc_rcut_mcoh_4_ztd = [Path(str(x)) for x in self.dsc_pair_path.glob(
            'Main/Processed_Data/20*20*/*displacement*rcut_mcoh.tif')]
        self.dsc_rcut_mcoh_4_ztd.sort()

        self.coh_path_1 = [Path(str(x)) for x in self.base_path_all_dirs.glob('Pair_' + self.path_1_full_name +
                                                                              '/coherence_mask/coherence_mask.tif')][0]
        self.coh_path_2 = [Path(str(x)) for x in self.base_path_all_dirs.glob('Pair_' + self.path_2_full_name +
                                                                              '/coherence_mask/coherence_mask.tif')][0]

    def make_incidence_files(self):
        mk_in.incidence_tif_file_creation(self.meta_path_asc,
                                          self.asc_rcut_4_ztd_base,
                                          self.asc_mcoh_4_ztd_base,
                                          self.path_2_asc_ztd,
                                          'asc')
        mk_in.incidence_tif_file_creation(self.meta_path_dsc,
                                          self.dsc_rcut_4_ztd_base,
                                          self.dsc_mcoh_4_ztd_base,
                                          self.path_2_dsc_ztd,
                                          'dsc')

    def apply_ztd(self):
        for base_rcut_mcoh_displacement_path in self.asc_rcut_mcoh_4_ztd:
            apz.apply_gacos_correction(base_rcut_mcoh_displacement_path,
                                       self.path_2_inc_asc_ztd,
                                       self.path_2_asc_ztd)
        for base_rcut_mcoh_displacement_path in self.dsc_rcut_mcoh_4_ztd:
            apz.apply_gacos_correction(base_rcut_mcoh_displacement_path,
                                       self.path_2_inc_dsc_ztd,
                                       self.path_2_dsc_ztd)

    def make_asc_dsc_txt_files_mcoh(self):
        self.asc_mcoh_all = [Path(str(x)) for x in self.asc_pair_path.glob(
            'Main/Processed_Data/20*20*/20*_20*_displacement_VV_rcut_mcoh.tif')]
        self.dsc_mcoh_all = [Path(str(x)) for x in self.dsc_pair_path.glob(
            'Main/Processed_Data/20*20*/20*_20*_displacement_VV_rcut_mcoh.tif')]
        self.asc_mcoh_gacos_all = [Path(str(x)) for x in self.asc_pair_path.glob(
            'Main/Processed_Data/20*20*/20*_20*_displacement_VV_rcut_mcoh_gacos.tif')]
        self.dsc_mcoh_gacos_all = [Path(str(x)) for x in self.dsc_pair_path.glob(
            'Main/Processed_Data/20*20*/20*_20*_displacement_VV_rcut_mcoh_gacos.tif')]
        loop_displacement_files_mcoh(self.asc_displacement_tif_pairs_list_cut,
                                     self.base_path,
                                     self.asc_pair_path,
                                     self.time_series_dir,
                                     self.asc_flag,
                                     self.num_of_stds,
                                     self.asc_mcoh_all,
                                     self.asc_mcoh_gacos_all,
                                     self.num_cores,
                                     self.gacos_or_non_gacos_asc,
                                     self.gacos_or_non_gacos_full,
                                     self.std_E)
        loop_displacement_files_mcoh(self.dsc_displacement_tif_pairs_list_cut,
                                     self.base_path,
                                     self.dsc_pair_path,
                                     self.time_series_dir,
                                     self.dsc_flag,
                                     self.num_of_stds,
                                     self.dsc_mcoh_all,
                                     self.dsc_mcoh_gacos_all,
                                     self.num_cores,
                                     self.gacos_or_non_gacos_dsc,
                                     self.gacos_or_non_gacos_full,
                                     self.std_E)

    def add_mcoh_tif_names(self):
        asc_mcoh_gacos_list = (general_functions.open_json_file
        (self.parsed_list_mcoh_asc_name)['rel_coh_path_list_updated'])
        dsc_mcoh_gacos_list = (general_functions.open_json_file
        (self.parsed_list_mcoh_dsc_name)['rel_coh_path_list_updated'])
        self.asc_coh_tif_pairs_list_mcoh_gacos = [Path(x) for x in asc_mcoh_gacos_list]
        self.dsc_coh_tif_pairs_list_mcoh_gacos = [Path(x) for x in dsc_mcoh_gacos_list]

    def make_ref_point(self):
        if not self.time_series_parameters_dir.exists():
            os.mkdir(self.time_series_parameters_dir)
        base_tif_use = self.asc_coh_tif_pairs_list_mcoh_gacos[0]
        self.master_width = rasterio_basic_functions.get_tif_width(base_tif_use)
        self.master_height = rasterio_basic_functions.get_tif_height(base_tif_use)
        if self.ref_type == 'gps':
            self.c_flag_col, self.c_flag_row = (rasterio_basic_functions.
                                                lon_lat_2_col_row(base_tif_use, self.gps_lon, self.gps_lat))
            gps_dict = {'lat_lon': [str(self.gps_lon), str(self.gps_lat)]}
            general_functions.write_json_from_dict(gps_dict, self.time_series_parameters_save_file)
        elif self.ref_type == 'auto':
            self.max_coh, self.c_flag_row, self.c_flag_col = (auto_find_reference_point
                                                              (self.coherence_mask_path,
                                                               self.half_size,
                                                               self.num_of_ref,
                                                               self.coh_min,
                                                               self.num_of_coh_blocks,
                                                               self.num_of_max_to_check))
            dict_lon_lat = {}
            hh = 1
            for col, row in zip(self.c_flag_col, self.c_flag_row):
                lon, lat = rasterio_basic_functions.col_row_2_lon_lat(base_tif_use, col, row)
                dict_lon_lat[('lat_lon_' + str(hh))] = [str(lat), str(lon)]
                hh = hh + 1
            general_functions.write_json_from_dict(dict_lon_lat, self.time_series_parameters_save_file)
        elif self.ref_type == 'overlap':
            overlap_gps_points = find_max_coh_lon_lat.gps_reference(self.coh_path_1, self.coh_path_2, self.num_of_ref)
            lon, lat = overlap_gps_points.run_gps_reference()
            self.c_flag_col, self.c_flag_row = rasterio_basic_functions.lon_lat_2_col_row(base_tif_use, lon, lat)
            overlap_dict = {'lat_lon': [str(lon), str(lat)]}
            general_functions.write_json_from_dict(overlap_dict, self.time_series_parameters_save_file)
        else:
            self.c_flag_col = None
            self.c_flag_row = None

    def run_time_series_set_up(self):
        if self.run_resample == 'true':
            resample_object = Resample_2_same_grid.resample_2_common_grid(self.asc_dsc_displacement_tif_combined_list,
                                                                          self.asc_dsc_coherence_tif_combined_list,
                                                                          self.base_path,
                                                                          self.coherence_threshold)
            resample_object.run_resample()
        self.read_asc_dsc_parameters()
        self.add_new_tif_names()
        if self.run_dtzd == 'true':
            self.make_incidence_files()
            self.apply_ztd()
        self.make_asc_dsc_txt_files_mcoh()
        self.add_mcoh_tif_names()
        self.make_ref_point()
        header_txt_file_creation.header_txt_file_creation_f(self.lambda__, self.L, self.master_width,
                                                            self.master_height, self.asc_time, self.asc_heading,
                                                            self.asc_incidence, self.dsc_time, self.dsc_heading,
                                                            self.dsc_incidence, self.time_series_dir,
                                                            self.c_size, self.c_flag_col, self.c_flag_row,
                                                            self.num_of_ref,
                                                            self.ref_type)
        stack_object = stack_processing.make_stack_asc_dsc(self.stack_path, self.asc_coh_tif_pairs_list_mcoh_gacos,
                                                           self.dsc_coh_tif_pairs_list_mcoh_gacos,
                                                           self.asc_incidence,
                                                           self.dsc_incidence,
                                                           self.c_flag_col,
                                                           self.c_flag_row,
                                                           self.c_size,
                                                           self.coherence_mask_path,
                                                           self.num_of_ref,
                                                           self.coherence_threshold)
        stack_object.run_functions()


if __name__ == '__main__':
    sys_index_var_1 = sys.argv[1]
    sys_index_var_2 = sys.argv[2]
    sys_index_var_3 = sys.argv[3]
    sys_index_var_4 = sys.argv[4]
    sys_index_var_5 = sys.argv[5]
    sys_index_var_6 = sys.argv[6]
    sys_index_var_7 = int(sys.argv[7])
    sys_index_var_8 = sys.argv[8]
    sys_index_var_9 = sys.argv[9]
    sys_index_var_10 = sys.argv[10]
    run_time_series_object = msbasv3_time_series_set_up(sys_index_var_1, sys_index_var_2, sys_index_var_3,
                                                        sys_index_var_4, ref_type=sys_index_var_5,
                                                        run_resample=sys_index_var_6, ref_points=sys_index_var_7,
                                                        run_dtzd=sys_index_var_8,
                                                        gacos_full=sys_index_var_9,
                                                        std_E=sys_index_var_10)
    run_time_series_object.run_time_series_set_up()
