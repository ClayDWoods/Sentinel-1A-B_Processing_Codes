import os, sys
from general_set_up_functions import rasterio_basic_functions as rb
from general_set_up_functions import general_functions as gf
from pathlib import Path


class gacos_set_up():
    def __init__(self, asc_name, dsc_name):
        self.base_path = Path.cwd().absolute().parent
        self.base_path_asc = self.base_path / asc_name
        self.base_path_dsc = self.base_path / dsc_name
        self.parameter_file_asc = self.base_path_asc / 'Main/Time_Series/platform_parameters.json'
        self.parameter_file_dsc = self.base_path_dsc / 'Main/Time_Series/platform_parameters.json'
        self.path_2_asc_ztd = self.base_path / 'MSBASv3_Timeseries/asc_ztd'
        self.path_2_dsc_ztd = self.base_path / 'MSBASv3_Timeseries/dsc_ztd'
        self.asc_file_download_name = self.path_2_asc_ztd / 'asc_downloads.txt'
        self.dsc_file_download_name = self.path_2_dsc_ztd / 'dsc_downloads.txt'
        self.processed_dirs_asc = [Path(str(x)) for x in
                                   self.base_path_asc.glob('Main/Processed_Data/20*20*/*displacement_VV_rcut_mcoh.tif')]
        self.processed_dirs_asc.sort()
        self.processed_dirs_dsc = [Path(str(x)) for x in
                                   self.base_path_dsc.glob('Main/Processed_Data/20*20*/*displacement_VV_rcut_mcoh.tif')]
        self.processed_dirs_dsc.sort()
        self.base_file_asc = self.processed_dirs_asc[0]
        self.base_file_dsc = self.processed_dirs_dsc[0]
        self.left_asc, self.bottom_asc, self.right_asc, self.top_asc = None, None, None, None
        self.left_dsc, self.bottom_dsc, self.right_dsc, self.top_dsc = None, None, None, None
        self.time_of_day_asc, self.time_of_day_dsc = None, None
        self.asc_dates_list, self.dsc_dates_list = None, None

    def make_dirs(self):
        if not self.path_2_asc_ztd.exists():
            os.mkdir(self.path_2_asc_ztd)
        if not self.path_2_dsc_ztd.exists():
            os.mkdir(self.path_2_dsc_ztd)

    def bounds_return(self):
        self.left_asc, self.bottom_asc, self.right_asc, self.top_asc = rb.tif_bounds_ztd(self.base_file_asc)
        self.left_dsc, self.bottom_dsc, self.right_dsc, self.top_dsc = rb.tif_bounds_ztd(self.base_file_dsc)

    def time_of_day_return(self):
        time_of_day_asc, time_of_day_dsc = (gf.open_json_file(self.parameter_file_asc)['startTime'],
                                            gf.open_json_file(self.parameter_file_dsc)['startTime'])
        self.time_of_day_asc, self.time_of_day_dsc = (gf.HHMMSS_2_HHMM(time_of_day_asc),
                                                      gf.HHMMSS_2_HHMM(time_of_day_dsc))

    def dates_return(self):
        self.asc_dates_list, self.dsc_dates_list = (gf.unique_dates(self.processed_dirs_asc),
                                                    gf.unique_dates(self.processed_dirs_dsc))

    def write_files(self):
        asc_list = [[self.left_asc, self.bottom_asc, self.right_asc, self.top_asc],
                    [self.time_of_day_asc],
                    self.asc_dates_list
                    ]
        dsc_list = [[self.left_dsc, self.bottom_dsc, self.right_dsc, self.top_dsc],
                    [self.time_of_day_dsc],
                    self.dsc_dates_list]
        gf.write_txt_file_gacos(self.asc_file_download_name, asc_list)
        gf.write_txt_file_gacos(self.dsc_file_download_name, dsc_list)

    def run_all(self):
        self.make_dirs()
        self.bounds_return()
        self.time_of_day_return()
        self.dates_return()
        self.write_files()

if __name__ == '__main__':
    sys_index_var_1 = sys.argv[1]
    sys_index_var_2 = sys.argv[2]
    run_gacos_set_up = gacos_set_up(sys_index_var_1, sys_index_var_2)
    run_gacos_set_up.run_all()



