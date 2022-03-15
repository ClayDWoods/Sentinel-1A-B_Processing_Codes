import sys, os
import os.path

sys.path.append("..")
import setup_code
import helper_functions as gen_fun
import rasterio_composite_functions as rc
from pathlib import Path
import pathlib
from math import floor
import xml.etree.ElementTree as ET
import numpy as np
import shutil
import datetime


class insar_processing_snap():
    # -------------------------------------------------------------------------------------------------------------------------------
    # 1.
    # exe paths
    base_snap_exe = "singularity exec /projects/clwo4142/containers/snap_8_0.sif gpt"
    base_snap_exe_no_aux = "singularity exec /projects/clwo4142/containers/snap_8_0_1.sif gpt"
    snaphu_exe = "/projects/clwo4142/SNAPHU_Install/snaphu-v2.0.4_linux/bin/snaphu"
    gdal_exe = "/projects/clwo4142/software/anaconda/envs/mycustomenv/bin/gdal_translate"
    exe_keys = ['base_snap_exe', 'base_snap_no_aux', 'snaphu_exe', 'gdal_exe']
    exe_vals = [base_snap_exe, base_snap_exe_no_aux, snaphu_exe, gdal_exe]
    #
    # 1.
    # -------------------------------------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------------------------------------
    # 2.
    # max seperation of pairs in days
    max_delta_t_keys = ['max_delta_t']
    max_delta_t = 24
    max_delta_t_vals = [max_delta_t]
    #
    # 2.
    # -------------------------------------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------------------------------------
    # 3.
    # number of threads to use
    processors_keys = ['number_of_processors', 'number_of_processors_snaphu']
    q = '-q 40 '
    qs = 40
    processors_vals = [q, qs]
    #
    # 3.
    # -------------------------------------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------------------------------------
    # 4.
    # orbit polynomial degree to use
    Orbit_Keys = ['orbit_degree_apply_orbit', 'orbit_degree_ifg', 'orbit_degree_tpr']
    Orbit_degree_apply_orbit, orbit_degree_ifg, orbit_degree_tpr = 3, 3, 3  # 1, 2, 3, 4 or 5 apply or ifg
    # 1-10 for tpr
    Orbit_Vals = [Orbit_degree_apply_orbit, orbit_degree_ifg, orbit_degree_tpr]
    #
    # 4.
    # -------------------------------------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------------------------------------
    # 5.
    #
    topsar_split_swaths_keys = ['swaths', 'burststart', 'burstend']
    burststart = [str(6), str(6), str(6)]
    burstend = [str(9), str(9), str(9)]
    topsar_split_swaths_vals = [['IW1', 'IW2', 'IW3'], burststart,
                                burstend]  # list of swaths to use IW123, IW12, IW23, IW1, IW2, IW3
    #
    # 5.
    # -------------------------------------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------------------------------------
    # 6.
    #
    polarisation_keys = ['polarisation']
    polarisation_vals = ['VV']  # polarisation to use
    #
    # 6.
    # -------------------------------------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------------------------------------
    # 7.
    #
    Back_Geocoding_Param_Keys = ['dem_resample_method', 'image_resample_method']
    # Parameter Values
    #                                  ------------------------------------------------------------------------------------------------------------------------
    # DEM_resample_method : <string>   Sets parameter 'demResamplingMethod' to <string>
    #                                  Default value is 'BICUBIC_INTERPOLATION'
    #                                  Value must be one of 'NEAREST_NEIGHBOUR', 'BILINEAR_INTERPOLATION', 'CUBIC_CONVOLUTION', 'BISINC_5_POINT_INTERPOLATION',
    #                                  'BISINC_11_POINT_INTERPOLATION', 'BISINC_21_POINT_INTERPOLATION', 'BICUBIC_INTERPOLATION'
    #                                  ------------------------------------------------------------------------------------------------------------------------
    # Image_Resample_Method : <string> The method to be used when resampling the slave grid onto the master grid
    #                                  Default value is 'BISINC_5_POINT_INTERPOLATION'
    #                                  Value must be one of 'NEAREST_NEIGHBOUR', 'BILINEAR_INTERPOLATION', 'CUBIC_CONVOLUTION', 'BISINC_5_POINT_INTERPOLATION',
    #                                  'BISINC_11_POINT_INTERPOLATION', 'BISINC_21_POINT_INTERPOLATION', 'BICUBIC_INTERPOLATION'
    #                                  ------------------------------------------------------------------------------------------------------------------------
    DEM_resample_method = 'BISINC_5_POINT_INTERPOLATION'
    Image_Resample_Method = 'BISINC_5_POINT_INTERPOLATION'
    # Parameter Values List
    Back_Geocoding_Param_Vals = [DEM_resample_method, Image_Resample_Method]
    #
    # 7.
    # -------------------------------------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------------------------------------
    #
    # 8.
    Enhanced_Spectral_Diversity_Param_Keys = ['fine_win_width', 'fine_win_height', 'fine_win_acc_azimuth',
                                              'fine_win_acc_range', 'fine_win_oversampling', 'xCorr_threshold',
                                              'coh_threshold', 'num_blocks_per_overlap']
    # Parameter Values
    #                                --------------------------------------------------------------------
    # Fine_Win_Width : <int>         Sets parameter 'fineWinWidthStr' to <string>
    #                                Value must be one of '32', '64', '128', '256', '512', '1024', '2048'
    #                                Default value is '512'
    #                                --------------------------------------------------------------------
    # Fine_Win_Height : <int>        Sets parameter 'fineWinHeightStr' to <string>
    #                                Value must be one of '32', '64', '128', '256', '512', '1024', '2048'
    #                                Default value is '512'
    #                                --------------------------------------------------------------------
    # Fine_Win_Acc_Azimuth : <int>   Sets parameter 'fineWinAccAzimuth' to <string>
    #                                Value must be one of '2', '4', '8', '16', '32', '64'
    #                                Default value is '16'
    #                                --------------------------------------------------------------------
    # Fine_Win_Acc_Range : <int>     Sets parameter 'fineWinAccRange' to <string>
    #                                Value must be one of '2', '4', '8', '16', '32', '64'
    #                                Default value is '16'
    #                                --------------------------------------------------------------------
    # Fine_Win_Oversampling : <int>  Sets parameter 'fineWinOversampling' to <string>
    #                                Value must be one of '32', '64', '128', '256'
    #                                Default value is '128'
    #                                --------------------------------------------------------------------
    # XCorr_Threshold : <double>     The peak cross-correlation threshold
    #                                Valid interval is (0, *)
    #                                Default value is '0.1'
    #                                --------------------------------------------------------------------
    # Coh_Threshold : <double>       The coherence threshold for outlier removal
    #                                Valid interval is (0, 1]
    #                                Default value is '0.3'
    #                                --------------------------------------------------------------------
    # Num_Blocks_Per_Overlap : <int> The number of windows per overlap for ESD
    #                                Valid interval is [1, 20]
    #                                Default value is '10'
    #                                --------------------------------------------------------------------
    Fine_Win_Width = 512
    Fine_Win_Height = 512
    Fine_Win_Acc_Azimuth = 16
    Fine_Win_Acc_Range = 16
    Fine_Win_Oversampling = 256
    XCorr_Threshold = 0.1
    Coh_Threshold = 0.3
    Num_Blocks_Per_Overlap = 10
    # Parameter Values List
    Enhanced_Spectral_Diversity_Param_Vals = [Fine_Win_Width, Fine_Win_Height, Fine_Win_Acc_Azimuth,
                                              Fine_Win_Acc_Range, Fine_Win_Oversampling, XCorr_Threshold,
                                              Coh_Threshold, Num_Blocks_Per_Overlap]
    #
    # 8.
    # -------------------------------------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------------------------------------
    #
    # 9.
    Interferogram_Param_Keys = ['srp_polynomial_degree', 'srp_number_points']
    # Parameter Values
    #                               -----------------------------------------------------------------------------
    # Srp_Polynomial_Degree : <int> Order of 'Flat earth phase' polynomial
    #                               Value must be one of '1', '2', '3', '4', '5', '6', '7', '8'
    #                               Default value is '5'
    #                               -----------------------------------------------------------------------------
    # Srp_Number_Points : <int>     Number of points for the 'flat earth phase' polynomial estimation
    #                               Value must be one of '301', '401', '501', '601', '701', '801', '901', '1001'
    #                               Default value is '501'
    #                               -----------------------------------------------------------------------------
    Srp_Polynomial_Degree = 6
    Srp_Number_Points = 601

    # Parameter Values List
    Interferogram_Param_Vals = [Srp_Polynomial_Degree, Srp_Number_Points]
    #
    # 9.
    # -------------------------------------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------------------------------------
    #
    # 10.
    Multilook_Param_Keys = ['num_rg_Looks']
    # Parameter Values
    #                 ---------------------------------------------
    # Num_Rg_Looks : <int> The user defined number of azimuth looks
    #                 Valid interval is [1, *)
    #                 Default value is '1'
    #                 ---------------------------------------------
    # Num_Az_Looks : <int> The user defined number of range looks
    #                 Valid interval is [1, *)
    #                 Default value is '1'
    #                 ---------------------------------------------
    Num_Rg_Looks = 28
    # Num_Az_Looks = 3
    # Parameter Values List
    Multilook_Param_Vals = [Num_Rg_Looks]
    #
    # 10.
    # -------------------------------------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------------------------------------
    #
    # 11.
    Goldstein_Phase_Filter_Param_Keys = ['alpha', 'fft_size', 'window_size', 'coherence_mask', 'coh_threshold_gpf']
    # Parameter Values
    #                                ---------------------------------------------
    # Alpha : <double>               Adaptive filter exponent
    #                                Valid interval is (0, 1]
    #                                Default value is '1.0'
    #                                ---------------------------------------------
    # FFT_Size : <string>            Sets parameter 'FFTSizeString' to <string>
    #                                Value must be one of '32', '64', '128', '256'
    #                                Default value is '64'
    #                                ---------------------------------------------
    # Window_Size : <string>         Sets parameter 'windowSizeString' to <string>
    #                                Value must be one of '3', '5', '7'
    #                                Default value is '3'
    #                                ---------------------------------------------
    # Coherence_Mask : <boolean>     Use coherence mask
    #                                Default value is 'false'
    #                                ---------------------------------------------
    # Coherence_Threshold : <double> The coherence threshold
    #                                Valid interval is [0, 1]
    #                                Default value is '0.2'
    #                                ---------------------------------------------
    Alpha = 1
    FFT_Size = 256
    Window_Size = 7
    Coherence_Mask = 'false'
    Coherence_Threshold_gpf = 0.2
    # Parameter Values List
    Goldstein_Phase_Filter_Param_Vals = [Alpha, FFT_Size, Window_Size, Coherence_Mask, Coherence_Threshold_gpf]
    #
    # 11.
    # -------------------------------------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------------------------------------
    # 12.
    #
    # -PgeoRegion=<geometry>                            The subset region in geographical coordinates using WKT-format,
    #                                                     e.g. POLYGON((<lon1> <lat1>, <lon2> <lat2>, ..., <lon1> <lat1>))
    #                                                     (make sure to quote the option due to spaces in <geometry>).
    #                                                     If not given, the entire scene is used.
    # POLYGON((<lon1> <lat1>, <lon2> <lat2>, ..., <lon1> <lat1>))
    #   -PsubSamplingX=<int>                              The pixel sub-sampling step in X (horizontal image direction)
    #                                                     Default value is '1'.
    #   -PsubSamplingY=<int>                              The pixel sub-sampling step in Y (vertical image direction)
    #                                                     Default value is '1'.

    lon_1, lat_1, lat_2, lon_3 = None, None, None, None
    if lon_1 == None and lat_1 == None and lat_2 == None and lon_3 == None:
        spatial_subset_keys = ['wkt_polygon', 'subSamplingX', 'subSamplingY']
        subSamplingX, subSamplingY = ' ', ' '
        spatial_subset_vals = [' ', subSamplingX, subSamplingY]
    else:
        spatial_subset_keys = ['wkt_polygon', 'subSamplingX', 'subSamplingY']
        subSamplingX, subSamplingY = 1, 1
        spatial_subset_vals = 'POLYGON((' + lon_1 + ' ' + lat_1 + ', ' + lon_1 + ' ' + lat_2 + ', ' + lon_3 + ' ' + \
                              lat_2 + ', ' + lon_3 + ' ' + lat_1 + ', ' + lon_1 + ' ' + lat_1 + '))'

    #
    # 12.
    # -------------------------------------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------------------------------------
    #
    # 13.
    Snaphu_Export_Param_Keys = ['stat_cost_mode', 'init_method', 'number_of_tile_rows', 'number_of_tile_cols',
                                'row_over_lap', 'col_over_lap', 'tile_cost_threshold']
    # Parameter Values
    #                              -----------------------------------------------------------------------------------------
    # Stat_Cost_Mode : <string>    Size of coherence estimation window in Azimuth direction
    #                              Value must be one of 'TOPO', 'DEFO', 'SMOOTH', 'NOSTATCOSTS'
    #                              Default value is 'DEFO'
    #                              -----------------------------------------------------------------------------------------
    # Init_Method : <string>       Algorithm used for initialization of the wrapped phase values
    #                              Value must be one of 'MST', 'MCF'
    #                              Default value is 'MST'
    #                              -----------------------------------------------------------------------------------------
    # Number_Of_Tile_Rows : <int>  Divide the image into tiles and process in parallel. Set to 1 for single tiled
    #                              Default value is '10'
    #                              -----------------------------------------------------------------------------------------
    # Number_Of_Tile_Cols : <int>  Divide the image into tiles and process in parallel. Set to 1 for single tiled
    #                              Default value is '10'
    #                              -----------------------------------------------------------------------------------------
    # Number_Of_Processors : <int> Number of concurrent processing threads. Set to 1 for single threaded
    #                              Default value is '4'
    #                              -----------------------------------------------------------------------------------------
    # Row_Over_Lap : <int>         Overlap, in pixels, between neighboring tiles
    #                              Default value is '200'
    #                              -----------------------------------------------------------------------------------------
    # Col_Over_Lap : <int>         Overlap, in pixels, between neighboring tiles
    #                              Default value is '200'
    #                              -----------------------------------------------------------------------------------------
    # Tile_Cost_Threshold  <int>   Cost threshold to use for determining boundaries of reliable regions
    #                              (long, dimensionless; scaled according to other cost constants)
    #                              Larger cost threshold implies smaller regions---safer, but more expensive computationally
    #                              Default value is '500'
    #                              -----------------------------------------------------------------------------------------
    Stat_Cost_Mode = 'DEFO'
    Init_Method = 'MST'
    Number_Of_Tile_Rows = 5
    Number_Of_Tile_Cols = 6
    Row_Over_Lap = 75
    Col_Over_Lap = 75
    Tile_Cost_Threshold = 800
    # Parameter Values List
    snaphu_export_param_vals = [Stat_Cost_Mode, Init_Method, Number_Of_Tile_Rows, Number_Of_Tile_Cols,
                                Row_Over_Lap, Col_Over_Lap, Tile_Cost_Threshold]
    #
    # 13.
    # -------------------------------------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------------------------------------
    # default cache size
    cache_size_keys = ['tile_cache_size', 'disable_tile_cache', 'use_file_tile_cache', 'tie_point_geocoding',
                       'pixel_geocoding_max_fraction']
    tile_cache_size = 2048
    disable_tile_cache = 'true'
    use_file_tile_cache = 'true'
    tie_point_geocoding = 'true'
    pixel_geocoding_fraction = 'true'
    tile_cache_size_vals = [tile_cache_size, disable_tile_cache, use_file_tile_cache, tie_point_geocoding,
                            pixel_geocoding_fraction]
    #
    # -------------------------------------------------------------------------------------------------------------------------------

    # Parameter Keys and Values list
    param_keys_all = [exe_keys, max_delta_t_keys, processors_keys, Orbit_Keys, topsar_split_swaths_keys,
                      polarisation_keys, Back_Geocoding_Param_Keys,
                      Enhanced_Spectral_Diversity_Param_Keys, Interferogram_Param_Keys,
                      Multilook_Param_Keys, Goldstein_Phase_Filter_Param_Keys, spatial_subset_keys,
                      Snaphu_Export_Param_Keys, cache_size_keys]
    param_vals_all = [exe_vals, max_delta_t_vals, processors_vals, Orbit_Vals, topsar_split_swaths_vals,
                      polarisation_vals, Back_Geocoding_Param_Vals,
                      Enhanced_Spectral_Diversity_Param_Vals, Interferogram_Param_Vals,
                      Multilook_Param_Vals, Goldstein_Phase_Filter_Param_Vals, spatial_subset_vals,
                      snaphu_export_param_vals, tile_cache_size_vals]
    # -------------------------------------------------------------------------------------------------------------------------------
    # To be passed to main function to input parameters into each processing step
    all_Parameters = gen_fun.make_dictionary_from_params(param_keys_all, param_vals_all)

    def __init__(self, sys_argv_int, param_dict=all_Parameters):
        self.sys_argv_int = int(sys_argv_int)
        self.param_dict = param_dict
        self.base_path = Path.cwd().absolute()
        self.base_path_name = self.base_path.name
        self.pair_path = self.base_path.parent
        self.aux_poe = '/gpfs/summit/scratch/clwo4142/SNAP_RUNS/Snap_Runs/orbits/aux_poeorb'
        self.res_orb = '/gpfs/summit/scratch/clwo4142/SNAP_RUNS/Snap_Runs/orbits/aux_resorb'
        self.data_after_processed_path = self.pair_path / (self.base_path_name + '_data_temp')
        self.data_path = self.base_path / 'Main/Data'
        self.pairs_path = self.base_path / 'Main/Pairs'
        try:
            self.DEM_path = [str(x) for x in self.base_path.glob('Main/DEM/*')][0]
        except IndexError:
            self.DEM_path = [str(x) for x in self.base_path.glob('Main/DEM/')][0]
        self.meta_data_path = self.data_path / 'Meta_Data'
        self.kml_frame_coords_kml_path = self.data_path / 'kml_frame'
        self.snaphu_logs_path = self.base_path / 'Main/Snaphu_Logs'
        self.time_series_path = self.base_path / 'Main/Time_Series'
        self.json_pairs_list_path = self.data_path / 'Pairs_List.json'
        self.pair_fp = self.pairs_path / ("Pair_" + str(self.sys_argv_int) + ".json")
        self.processed_data_path = self.base_path / 'Main/Processed_Data'
        self.number_of_processors = self.param_dict['number_of_processors']
        self.number_of_processors_snaphu = self.param_dict['number_of_processors_snaphu']
        self.orbit_degree_apply_orbit = str(self.param_dict['orbit_degree_apply_orbit'])
        self.orbit_degree_ifg = str(self.param_dict['orbit_degree_ifg'])
        self.orbit_degree_tpr = str(self.param_dict['orbit_degree_tpr'])
        if self.sys_argv_int == 0:
            self.base_snap_exe = self.param_dict['base_snap_exe']
        else:
            self.base_snap_exe = self.param_dict['base_snap_no_aux']
        if os.path.isfile(self.DEM_path):
            self.DEM_to_use = self.DEM_path
            self.DEM_param = 'DEM_external'
            self.dem_p = '-PexternalDEMFile=' + self.DEM_path + ' '
        else:
            self.DEM_param = 'DEM_auto'
            self.dem_p = '-PdemName=\'SRTM 1Sec HGT\' '
        self.gdal_exe = self.param_dict['gdal_exe']
        self.snaphu_exe = self.param_dict['snaphu_exe']
        self.swaths = self.param_dict['swaths']
        self.burst_start = self.param_dict['burststart']
        self.burst_end = self.param_dict['burstend']
        self.polarisation = self.param_dict['polarisation']
        self.dem_resample_method = self.param_dict['dem_resample_method']
        self.image_resample_method = self.param_dict['image_resample_method']
        self.fine_win_width = self.param_dict['fine_win_width']
        self.fine_win_height = self.param_dict['fine_win_height']
        self.fine_win_acc_azimuth = self.param_dict['fine_win_acc_azimuth']
        self.fine_win_acc_range = self.param_dict['fine_win_acc_range']
        self.fine_win_oversampling = self.param_dict['fine_win_oversampling']
        self.xCorr_threshold = self.param_dict['xCorr_threshold']
        self.coh_threshold = self.param_dict['coh_threshold']
        self.coh_threshold_gpf = self.param_dict['coh_threshold_gpf']
        self.num_blocks_per_overlap = self.param_dict['num_blocks_per_overlap']
        self.srp_polynomial_degree = self.param_dict['srp_polynomial_degree']
        self.srp_number_points = self.param_dict['srp_number_points']
        self.num_rg_Looks = self.param_dict['num_rg_Looks']
        self.num_az_looks = None
        self.alpha = self.param_dict['alpha']
        self.fft_size = self.param_dict['fft_size']
        self.window_size = self.param_dict['window_size']
        self.coherence_mask = self.param_dict['coherence_mask']
        self.wkt_polygon = self.param_dict['wkt_polygon']
        self.subSamplingX = self.param_dict['subSamplingX']
        self.subSamplingY = self.param_dict['subSamplingY']
        self.stat_cost_mode = self.param_dict['stat_cost_mode']
        self.init_method = self.param_dict['init_method']
        self.number_of_tile_rows = self.param_dict['number_of_tile_rows']
        self.number_of_tile_cols = self.param_dict['number_of_tile_cols']
        self.row_over_lap = self.param_dict['row_over_lap']
        self.col_over_lap = self.param_dict['col_over_lap']
        self.tile_cost_threshold = self.param_dict['tile_cost_threshold']
        self.tile_cache_size = self.param_dict['tile_cache_size']
        self.disable_tile_cache = self.param_dict['disable_tile_cache']
        self.use_file_tile_cache = self.param_dict['use_file_tile_cache']
        self.tie_point_geocoding = self.param_dict['tie_point_geocoding']
        self.pixel_geocoding_fraction = self.param_dict['pixel_geocoding_max_fraction']
        self.slc_pair = None
        self.pair_dates = None
        self.product_output_directory = None
        self.orbit_corrected_file_paths = None
        self.topsar_split_paths = None
        self.bg_paths = None
        self.esd_paths = None
        self.int_paths = None
        self.deb_paths = None
        self.mrg_c_path = None
        self.tpr_c_path = None
        self.mlt_c_path = None
        self.gpf_c_path = None
        self.sbst_c_path = None
        self.sne_c_path = None
        self.snu_c_path = None
        self.sni_c_path = None
        self.p2d_c_path = None
        self.TC_c_path_1 = None
        self.TC_c_path_2 = None
        self.displacement_path_tif = None
        self.dem_name = None
        self.unw_phase_path_tif = None
        self.w_phase_path_tif = None
        self.coherence_path_tif = None
        self.intensity_db_path_tif = None
        self.intensity_path_tif = None
        self.local_incidence_path_tif = None
        self.time_to_run = None
        self.pair_cache_directory_name = None
        self.pair_cache_directory_path_java = None
        self.pair_cache_directory_path_aux = None
        self.pair_cache_directory_path_aux_primary = None
        self.java_cache_command = None
        self.aux_cache_command = None
        self.tile_cache_size_command = None
        self.use_file_tile_cache_command = None
        self.disable_tile_cache_command = None
        self.tie_point_geocoding_command = None
        self.pixel_geocoding_fraction_command = None
        self.gpt_vmoptions_extra = None
        self.gpt_vmoptions_esd = None
        self.list_of_all_aux_caches = None
        self.list_of_all_aux_caches_paths = None

    def num_of_looks(self):
        if not self.data_after_processed_path:
            self.data_after_processed_path.mkdir()
        if not self.time_series_path.exists():
            self.time_series_path.mkdir()
        heading, incidence_mid_swath, range_pixel_spacing, \
        azimuth_pixel_spacing, start_times = ('platformHeading',
                                              'incidenceAngleMidSwath',
                                              'rangePixelSpacing',
                                              'azimuthPixelSpacing',
                                              'startTime')
        lower_case_polarisation = "*" + self.polarisation.lower() + "*"
        headings, incidences, range_d, azimuth_d, start_time = [], [], [], [], []
        for swath in self.swaths:
            headings_t, incidences_t, range_d_t, azimuth_d_t, start_time_t = [], [], [], [], []
            swath_lwcase = swath.lower()
            swath_xml = [x for x in map(str, list(self.meta_data_path.glob(lower_case_polarisation))) if
                         swath_lwcase in x]
            for swath_xml_file in swath_xml:
                tree_xml = ET.parse(swath_xml_file)
                root_xml = tree_xml.getroot()
                swath_xml_file_path = Path(swath_xml_file)
                time_start = swath_xml_file_path.name.split('-')[4].split('t')[1]
                for heading1, incidence1, rgspacing, azspacing in zip(root_xml.iter(heading),
                                                                      root_xml.iter(incidence_mid_swath),
                                                                      root_xml.iter(range_pixel_spacing),
                                                                      root_xml.iter(azimuth_pixel_spacing)):
                    headings_t.append(float(heading1.text))
                    incidences_t.append(float(incidence1.text))
                    range_d_t.append(float(rgspacing.text))
                    azimuth_d_t.append(float(azspacing.text))
                    start_time_t.append(gen_fun.get_sec(time_start))
                del tree_xml, root_xml
            headings_t_a, incidences_t_a, range_d_t_a, azimuth_d_t_a, start_time_d_t_a = (np.mean(headings_t),
                                                                                          np.mean(incidences_t),
                                                                                          np.mean(range_d_t),
                                                                                          np.mean(azimuth_d_t),
                                                                                          np.mean(start_time_t))

            headings.append(headings_t_a)
            incidences.append(incidences_t_a)
            range_d.append(range_d_t_a)
            azimuth_d.append(azimuth_d_t_a)
            start_time.append(start_time_d_t_a)
        headings_mean, incidences_mean, range_sp_mean, azimuth_sp_mean, start_time_mean = (np.mean(headings),
                                                                                           np.mean(incidences),
                                                                                           np.mean(range_d),
                                                                                           np.mean(azimuth_d),
                                                                                           np.mean(start_time))
        start_time_hh_mm_ss = gen_fun.secondsToTime(start_time_mean)
        delta_rg = (range_sp_mean / np.sin(incidences_mean * (np.pi / 180))) * self.num_rg_Looks
        nb_look = int(round(delta_rg / azimuth_sp_mean))
        dict_4_header_params = {'platformHeading': headings_mean, 'incidenceAngleMidSwath': incidences_mean,
                                'rangePixelSpacing': range_sp_mean, 'azimuthPixelSpacing': azimuth_sp_mean,
                                'startTime': start_time_hh_mm_ss}
        dict_4_time_series_name = self.time_series_path / "platform_parameters.json"
        gen_fun.write_json_from_dict(dict_4_header_params, dict_4_time_series_name)
        self.num_az_looks = nb_look

    def write_param_file(self):
        param_file_dir = self.base_path / 'Main/Parameter_File'
        param_file_name = param_file_dir / 'parameter_file.json'
        if not param_file_dir.exists():
            param_file_dir.mkdir()
        gen_fun.write_json_from_dict(self.param_dict, param_file_name)

    def add_pairs_path(self):
        self.slc_pair = list(gen_fun.open_json_file(self.pair_fp).values())[0]
        slc_date_list = [gen_fun.slc_path_2_date(x) for x in self.slc_pair]
        self.pair_dates = (slc_date_list[0] + '_' + slc_date_list[1])
        self.product_output_directory = self.processed_data_path / self.pair_dates

    def add_cache_paths_commands(self):
        self.pair_cache_directory_name = self.pair_dates + '_cache'
        if self.sys_argv_int == 0:
            self.pair_cache_directory_path_aux_primary = self.base_path / ('Main/Cache/' +
                                                                           self.pair_cache_directory_name + '/.esa_snap')
            self.list_of_all_aux_caches = [str(x) for x in self.base_path.glob('Main/Cache/*_cache')]
            self.list_of_all_aux_caches.sort()
            self.list_of_all_aux_caches_paths = [Path(x + '/.esa_snap') for x in self.list_of_all_aux_caches]
            self.list_of_all_aux_caches_paths.remove(self.pair_cache_directory_path_aux_primary)

        self.pair_cache_directory_path_java = self.base_path / ('Main/Cache/' +
                                                                self.pair_cache_directory_name + '/java_tmp')
        self.pair_cache_directory_path_aux = self.base_path / ('Main/Cache/' +
                                                               self.pair_cache_directory_name + '/.esa_snap')
        self.java_cache_command = '-J-Djava.io.tmpdir=' + str(self.pair_cache_directory_path_java) + ' '
        self.aux_cache_command = '-J-Dsnap.userdir=' + str(self.pair_cache_directory_path_aux) + ' '
        self.tile_cache_size_command = '-J-Dsnap.jai.tileCacheSize=' + str(self.tile_cache_size) + 'M' + ' '
        self.use_file_tile_cache_command = '-J-Dsnap.gpf.useFileTileCache=' + self.disable_tile_cache + ' '
        self.disable_tile_cache_command = '-J-Dsnap.gpf.disableTileCache=' + self.use_file_tile_cache + ' '
        self.tie_point_geocoding_command = '-J-Dsnap.tiePointGeoCoding.maxPrecision=' + self.tie_point_geocoding + ' '
        self.pixel_geocoding_fraction_command = ('-J-Dsnap.pixelGeoCoding.fractionAccuracy=' +
                                                 self.pixel_geocoding_fraction + ' ')
        self.gpt_vmoptions_extra = (self.java_cache_command + self.aux_cache_command + self.tile_cache_size_command +
                                    self.tie_point_geocoding_command + self.pixel_geocoding_fraction_command)
        self.gpt_vmoptions_esd = (self.java_cache_command + self.aux_cache_command + self.tile_cache_size_command +
                                  self.use_file_tile_cache_command + self.disable_tile_cache_command +
                                  self.tie_point_geocoding_command + self.pixel_geocoding_fraction_command)

    def write_time_file(self):
        time_file_dir = self.base_path / 'Main/Time_to_run'
        time_file_name = time_file_dir / ('time_to_run_' + self.pair_dates + '.txt')
        if not time_file_dir.exists():
            time_file_dir.mkdir()
        gen_fun.write_txt_file(time_file_name, self.time_to_run)

    def applyOrbit(self):
        ext_name = '_Orb'
        orbit_corrected_file_paths, commands = [], []
        for slc in self.slc_pair:
            slc_date = gen_fun.slc_path_2_date(slc)
            ao_flag = (' Apply-Orbit-File -e -J-DOrbitFiles.sentinel1POEOrbitPath=' + self.aux_poe + ' ' +
                       '-J-DOrbitFiles.sentinel1RESOrbitPath=' + self.res_orb + ' ' +
                       ' ' + self.gpt_vmoptions_extra)
            o_type = ('-PcontinueOnFail=false' + ' ' +
                      '-PpolyDegree=' + self.orbit_degree_apply_orbit + ' ' +
                      '-PorbitType=\'Sentinel Precise (Auto Download)\' ')
            out = '-t ' + str(self.product_output_directory / (slc_date + ext_name)) + ' '
            cmd = self.base_snap_exe + ao_flag + out + self.number_of_processors + o_type + slc
            orbit_corrected_file_path = str(self.product_output_directory / (slc_date + ext_name + '.dim'))
            commands.append(cmd)
            orbit_corrected_file_paths.append(orbit_corrected_file_path)
        print('Applying Precise Orbit file')
        print(commands)
        for command in commands:
            gen_fun.subprocess(command)
        self.orbit_corrected_file_paths = orbit_corrected_file_paths
        print('Applying Precise Orbit file done')

    def topsar_split(self):
        ext_name = '_Split'
        topsar_split_paths, iw_split_command = [], []
        for Orbit_Path in self.orbit_corrected_file_paths:
            slc_date = str(Path(Orbit_Path).stem).split('_')[0]
            ao_flag = ' TOPSAR-Split -e ' + self.gpt_vmoptions_extra
            for swath, burst_start, burst_end in zip(self.swaths, self.burst_start, self.burstend):
                o_type = ('-PfirstBurstIndex=' + str(burst_start) + ' ' +
                          '-PlastBurstIndex=' + str(burst_end) + ' ' +
                          '-PselectedPolarisations=' + self.polarisation + ' ' +
                          '-Psubswath=' + swath + ' ')
                out = '-t ' + str(self.product_output_directory / (slc_date + '_' + swath + ext_name)) + ' '
                cmd = self.base_snap_exe + ao_flag + out + self.number_of_processors + o_type + Orbit_Path
                iw_split_command.append(cmd)
                topsar_split_corrected_file_path = (str(self.product_output_directory / (slc_date + '_' + swath +
                                                                                         ext_name + '.dim')))
                topsar_split_paths.append(topsar_split_corrected_file_path)
        print('Applying TOPSAR_Split')
        print(iw_split_command)
        for command in iw_split_command:
            gen_fun.subprocess(command)
        self.topsar_split_paths = topsar_split_paths
        print('Applying TOPSAR_Split done')

    def Back_Geocoding(self):
        ext_name = '_BG'
        bg_command, bg_paths = [], []
        ao_flag = (' Back-Geocoding -e -x -J-Dsnap.dataio.reader.tileWidth=512 '
                   '-J-Dsnap.dataio.reader.tileHeight=75 ' + self.gpt_vmoptions_extra)
        for swath in self.swaths:
            topsar_split_swath = [x for x in self.topsar_split_paths if swath in x]
            o_type = (self.dem_p +
                      '-PdemResamplingMethod=' + self.dem_resample_method + ' ' +
                      '-PresamplingType=' + self.image_resample_method + ' ')
            out = '-t ' + str(self.product_output_directory / (self.pair_dates + '_' + swath + ext_name)) + ' '
            sources = '-Ssources=' + topsar_split_swath[1] + ' ' + topsar_split_swath[0]
            cmd = self.base_snap_exe + ao_flag + out + self.number_of_processors + o_type + sources
            bg_command.append(cmd)
            bg_split_corrected_file_path = (str(self.product_output_directory / (self.pair_dates + '_' + swath +
                                                                                 ext_name + '.dim')))
            bg_paths.append(bg_split_corrected_file_path)
        print('Applying Back-Geocoding')
        print(bg_command)
        for command in bg_command:
            gen_fun.subprocess(command)
        self.bg_paths = bg_paths
        print('Applying Back-Geocoding done')

    def esd(self):
        ext_name = '_ESD'
        esd_command, esd_paths = [], []
        ao_flag = ' Enhanced-Spectral-Diversity -e -x ' + self.gpt_vmoptions_esd
        num_proc_esd = '-q ' + str(int(int(self.number_of_processors_snaphu) / 2)) + ' '
        for bg_path, swath in zip(self.bg_paths, self.swaths):
            o_type = ('-PfineWinWidthStr=' + str(self.fine_win_width) + ' ' +
                      '-PfineWinHeightStr=' + str(self.fine_win_height) + ' ' +
                      '-PfineWinAccAzimuth=' + str(self.fine_win_acc_azimuth) + ' ' +
                      '-PfineWinAccRange=' + str(self.fine_win_acc_range) + ' ' +
                      '-PfineWinOversampling=' + str(self.fine_win_oversampling) + ' ' +
                      '-PxCorrThreshold=' + str(self.xCorr_threshold) + ' ' +
                      '-PcohThreshold=' + str(self.coh_threshold) + ' ' +
                      '-PnumBlocksPerOverlap=' + str(self.num_blocks_per_overlap) + ' ')
            out = ('-t ' + str(self.product_output_directory / (self.pair_dates + '_' + swath + ext_name)) +
                   ' ')
            cmd = self.base_snap_exe + ao_flag + out + num_proc_esd + o_type + bg_path
            esd_command.append(cmd)
            esd_corrected_file_path = (str(self.product_output_directory / (self.pair_dates + '_' + swath +
                                                                            ext_name + '.dim')))
            esd_paths.append(esd_corrected_file_path)
        print('Applying ESD')
        print(esd_command)
        for command in esd_command:
            gen_fun.subprocess(command)
        self.esd_paths = esd_paths
        print('Applying ESD done')

    def make_interferogram(self):
        ext_name = '_Int'
        int_command, int_paths = [], []
        ao_flag = ' Interferogram -e ' + self.gpt_vmoptions_extra
        for esd_path, swath in zip(self.esd_paths, self.swaths):
            o_type = ('-PsrpPolynomialDegree=' + str(self.srp_polynomial_degree) + ' ' +
                      '-PsrpNumberPoints=' + str(self.srp_number_points) + ' ' +
                      '-PorbitDegree=' + self.orbit_degree_ifg + ' ' +
                      '-PcohWinAz=' + str(self.num_az_looks) + ' ' +
                      '-PcohWinRg=' + str(self.num_rg_Looks) + ' '
                      + self.dem_p + ' ')
            out = ('-t ' + str(self.product_output_directory / (self.pair_dates + '_' + swath +
                                                                ext_name)) + ' ')
            cmd = self.base_snap_exe + ao_flag + out + self.number_of_processors + o_type + esd_path
            int_command.append(cmd)
            int_corrected_file_path = (str(self.product_output_directory / (self.pair_dates + '_' + swath +
                                                                            ext_name + '.dim')))
            int_paths.append(int_corrected_file_path)
        print('Making Interferogram')
        print(int_command)
        for command in int_command:
            gen_fun.subprocess(command)
        self.int_paths = int_paths
        print('Making Interferogram done')

    def topsar_deburst(self):
        ext_name = '_Deb'
        deb_command, deb_paths = [], []
        ao_flag = ' TOPSAR-Deburst -e ' + self.gpt_vmoptions_extra
        for int_path, swath in zip(self.int_paths, self.swaths):
            o_type = ('-PselectedPolarisations=' + self.polarisation + ' ')
            out = ('-t ' + str(self.product_output_directory / (self.pair_dates + '_' + swath +
                                                                ext_name)) + ' ')
            cmd = self.base_snap_exe + ao_flag + out + self.number_of_processors + o_type + int_path
            deb_command.append(cmd)
            deb_corrected_file_path = (str(self.product_output_directory / (self.pair_dates + '_' + swath +
                                                                            ext_name + '.dim')))
            deb_paths.append(deb_corrected_file_path)
        print('Debursting')
        for command in deb_command:
            gen_fun.subprocess(command)
        self.deb_paths = deb_paths
        print('Debursting done')

    def merge_swaths(self):
        if len(self.swaths) > 1:
            ext_name = 'mrg'
            mrg_command = []
            ao_flag = ' TOPSAR-Merge -e ' + self.gpt_vmoptions_extra
            o_type = ('-PselectedPolarisations=' + self.polarisation + ' ')
            out = ('-t ' + str(self.product_output_directory / (self.pair_dates + '_' +
                                                                ext_name)) + ' ')
            if len(self.swaths) == 3:
                sources = '-Ssources=' + self.deb_paths[1] + ' ' + self.deb_paths[0] + ' ' + self.deb_paths[2]
                print(sources)
            else:
                sources = '-Ssources=' + self.deb_paths[1] + ' ' + self.deb_paths[0]
                print(sources)
            cmd = self.base_snap_exe + ao_flag + out + self.number_of_processors + o_type + sources
            print(cmd)
            mrg_c_path = (str(self.product_output_directory / (self.pair_dates + '_' +
                                                               ext_name + '.dim')))
            print(mrg_c_path)
            mrg_command.append(cmd)
            print('Merging Swaths')
            for command in mrg_command:
                gen_fun.subprocess(command)
            self.mrg_c_path = mrg_c_path
            print('Merging Swaths done')
        else:
            pass

    def TopoPhaseRemoval(self):
        tpr_command = []
        ao_flag = ' TopoPhaseRemoval -e ' + self.gpt_vmoptions_extra
        o_type = ('-PorbitDegree=' + self.orbit_degree_tpr + ' ' +
                  self.dem_p + ' ')
        ext_name = 'TPR'
        out = ('-t ' + str(self.product_output_directory / (self.pair_dates + '_' +
                                                            ext_name)) + ' ')
        tpr_c_path = (str(self.product_output_directory / (self.pair_dates + '_' + ext_name + '.dim')))
        if len(self.swaths) > 1:
            cmd = self.base_snap_exe + ao_flag + out + self.number_of_processors + o_type + self.mrg_c_path
        else:
            cmd = self.base_snap_exe + ao_flag + out + self.number_of_processors + o_type + self.deb_paths[0]
        tpr_command.append(cmd)
        print('Topo Phase Removal')
        print(tpr_command)
        for command in tpr_command:
            gen_fun.subprocess(command)
        self.tpr_c_path = tpr_c_path
        print('Topo Phase Removal done')

    def Multilook(self):
        mlt_command = []
        ao_flag = ' Multilook -e ' + self.gpt_vmoptions_extra
        o_type = ('-PnRgLooks=' + str(self.num_rg_Looks) + ' ' +
                  '-PnAzLooks=' + str(self.num_az_looks) + ' ' +
                  '-PgrSquarePixel=false ')
        ext_name = 'mlt'
        out = ('-t ' + str(self.product_output_directory / (self.pair_dates + '_' +
                                                            ext_name)) + ' ')
        mlt_c_path = (str(self.product_output_directory / (self.pair_dates + '_' + ext_name + '.dim')))
        cmd = self.base_snap_exe + ao_flag + out + self.number_of_processors + o_type + self.tpr_c_path
        mlt_command.append(cmd)
        print('Multilooking')
        print(mlt_command)
        for command in mlt_command:
            gen_fun.subprocess(command)
        self.mlt_c_path = mlt_c_path
        print('Multilooking done')

    def Goldstein_Phase_Filtering(self):
        GPF_command = []
        ao_flag = ' GoldsteinPhaseFiltering -e ' + self.gpt_vmoptions_extra
        o_type = ('-Palpha=' + str(self.alpha) + ' ' +
                  '-PFFTSizeString=' + str(self.fft_size) + ' ' +
                  '-PwindowSizeString=' + str(self.window_size) + ' ' +
                  '-PuseCoherenceMask=' + str(self.coherence_mask) + ' ' +
                  '-PcoherenceThreshold=' + str(self.coh_threshold_gpf) + ' ')
        ext_name = 'GPF'
        out = ('-t ' + str(self.product_output_directory / (self.pair_dates + '_' +
                                                            ext_name)) + ' ')
        gpf_c_path = (str(self.product_output_directory / (self.pair_dates + '_' + ext_name + '.dim')))
        cmd = self.base_snap_exe + ao_flag + out + self.number_of_processors + o_type + self.mlt_c_path
        GPF_command.append(cmd)
        print('Goldstein Phase Filtering')
        print(GPF_command)
        for command in GPF_command:
            gen_fun.subprocess(command)
        self.gpf_c_path = gpf_c_path
        print('Goldstein Phase Filtering done')

    def spatial_subset(self):
        if not self.wkt_polygon == ' ':
            subset_command = []
            ao_flag = ' subset -e ' + self.gpt_vmoptions_extra
            o_type = ('-PgeoRegion=' + str(self.wkt_polygon) + ' ' +
                      '-PsubSamplingX=' + str(self.subSamplingX) + ' ' +
                      '-PsubSamplingY=' + str(self.subSamplingY) + ' ')
            source = '-Ssource=' + self.gpf_c_path
            ext_name = 'sbst'
            out = ('-t ' + str(self.product_output_directory / (self.pair_dates + '_' +
                                                                ext_name)) + ' ')
            sbst_c_path = (str(self.product_output_directory / (self.pair_dates + '_' + ext_name + '.dim')))
            cmd = self.base_snap_exe + ao_flag + out + self.number_of_processors + o_type + source
            subset_command.append(cmd)
            print('Subset Image')
            print(subset_command)
            for command in subset_command:
                gen_fun.subprocess(command)
            self.sbst_c_path = sbst_c_path
            print('Subset Image done')
        else:
            pass

    def Snaphu_Export(self):
        ao_flag = ' SnaphuExport -e ' + self.gpt_vmoptions_extra
        ext_name = '_sne'
        snaphu_wd_c = self.product_output_directory / ('Snaphu_Output' + ext_name)
        snaphu_wd = str(self.product_output_directory / ('Snaphu_Output' + ext_name))
        if not snaphu_wd_c.exists():
            snaphu_wd_c.mkdir()
        o_type = ('-PtargetFolder=' + snaphu_wd + ' ' +
                  '-PstatCostMode=' + self.stat_cost_mode + ' ' +
                  '-PinitMethod=' + self.init_method + ' ' +
                  '-PnumberOfTileRows=' + str(self.number_of_tile_rows) + ' ' +
                  '-PnumberOfTileCols=' + str(self.number_of_tile_cols) + ' ' +
                  '-PnumberOfProcessors=' + str(self.number_of_processors_snaphu) + ' ' +
                  '-ProwOverlap=' + str(self.row_over_lap) + ' ' +
                  '-PcolOverlap=' + str(self.col_over_lap) + ' ' +
                  '-PtileCostThreshold=' + str(self.tile_cost_threshold) + ' ')
        if not self.wkt_polygon == ' ':
            source = '-Ssource=' + self.sbst_c_path
        else:
            source = '-Ssource=' + self.gpf_c_path
        snaph_e_command = self.base_snap_exe + ao_flag + o_type + source
        sne_c_path = str(self.product_output_directory / ('Snaphu_Output' + ext_name + '/' + self.pair_dates + '_GPF'))
        print('Snaphu Export')
        print(snaph_e_command)
        gen_fun.subprocess(snaph_e_command)
        self.sne_c_path = sne_c_path
        print('Snaphu Export done')

    def Snaphu_Unwrapping(self):
        cwd = Path.cwd()
        print(cwd)
        snaphu_conf = str(Path(self.sne_c_path) / 'snaphu.conf')
        snu_c_path = str(list(Path(self.sne_c_path).glob('UnwPhase_*.hdr'))[0])
        f = open(snaphu_conf, 'r')
        content = f.readlines()
        parsed_line = content[6].split('#')[1].split(' ')
        parsed_line = [x for x in parsed_line if x]
        snaphu_cmd_parse = (parsed_line[1] + ' ' + parsed_line[2] + ' ' + parsed_line[3] + ' ' +
                            parsed_line[4].rstrip("\n"))
        f.close()
        print(snaphu_cmd_parse)
        snu_command = self.snaphu_exe + ' ' + snaphu_cmd_parse
        print(snu_command)
        print('Unwrapping')
        path_2_unwrap = Path(self.sne_c_path)
        os.chdir(path_2_unwrap)
        gen_fun.subprocess(snu_command)
        os.chdir(cwd)
        self.snu_c_path = snu_c_path
        print('Unwrapping done')

    def Snaphu_Import(self):
        ao_flag = ' SnaphuImport -e ' + self.gpt_vmoptions_extra
        o_type = ('-PdoNotKeepWrapped=false ')
        ext_name = '_unw'
        if not self.wkt_polygon == ' ':
            source = '-Ssources=' + self.sbst_c_path + ' ' + self.snu_c_path
        else:
            source = '-Ssources=' + self.gpf_c_path + ' ' + self.snu_c_path
        out = '-t ' + str(Path(self.sne_c_path) / (self.pair_dates + ext_name)) + ' '
        sni_command = self.base_snap_exe + ao_flag + out + o_type + source
        sni_c_path = str(Path(self.sne_c_path) / (self.pair_dates + ext_name + '.dim'))
        print('Snaphu Import')
        print(sni_command)
        gen_fun.subprocess(sni_command)
        self.sni_c_path = sni_c_path
        print('Snaphu Import done')

    def Phase_2_Displacement(self):
        ao_flag = ' PhaseToDisplacement -e ' + self.gpt_vmoptions_extra
        source = '-SsourceProduct=' + self.sni_c_path
        ext_name = '_P2D'
        out = '-t ' + str(Path(self.sne_c_path) / (self.pair_dates + ext_name)) + ' '
        p2d_command = self.base_snap_exe + ao_flag + out + source
        p2d_c_path = str(Path(self.sne_c_path) / (self.pair_dates + ext_name + '.dim'))
        print('Phase to displacement')
        print(p2d_command)
        gen_fun.subprocess(p2d_command)
        self.p2d_c_path = p2d_c_path
        print('Phase to displacement done')

    def Terrain_Correction(self):
        ao_flag = ' Terrain-Correction -e ' + self.gpt_vmoptions_extra
        if self.DEM_param == 'DEM_external':
            o_type1 = (self.dem_p +
                       '-PdemName=\'SRTM 1Sec HGT\' ' +
                       '-PdemResamplingMethod=' + self.dem_resample_method + ' ' +
                       '-PimgResamplingMethod=' + self.image_resample_method + ' ')
        else:
            o_type1 = ('-PdemName=\'SRTM 1Sec HGT\' ' +
                       '-PdemResamplingMethod=' + self.dem_resample_method + ' ' +
                       '-PimgResamplingMethod=' + self.image_resample_method + ' ')
        source1 = '-SsourceProduct=' + self.p2d_c_path
        source2 = '-SsourceProduct=' + self.sni_c_path
        ext_name_1 = '_TC1'
        ext_name_2 = '_TC2'
        out1 = '-t ' + str(Path(self.sne_c_path) / (self.pair_dates + ext_name_1)) + ' '
        out2 = '-t ' + str(Path(self.sne_c_path) / (self.pair_dates + ext_name_2)) + ' '
        tc_command1 = self.base_snap_exe + ao_flag + out1 + self.number_of_processors + o_type1 + source1
        tc_command2 = self.base_snap_exe + ao_flag + out2 + self.number_of_processors + o_type1 + source2
        tc_c_path_1 = str(Path(self.sne_c_path) / (self.pair_dates + ext_name_1 + '.dim'))
        tc_c_path_2 = str(Path(self.sne_c_path) / (self.pair_dates + ext_name_2 + '.dim'))
        print('Terrain Correction')
        print(tc_command1)
        print(tc_command2)
        tc_command = [tc_command1, tc_command2]
        for command in tc_command:
            gen_fun.subprocess(command)
        self.TC_c_path_1 = tc_c_path_1
        self.TC_c_path_2 = tc_c_path_2
        print('Terrain Correction done')

    def img_2_geotiff(self):
        print("converting img files to geotiff")
        tc_paths_1_parse = list(Path(self.sne_c_path).glob("*_TC1.data"))[0]
        tc_paths_2_parse = list(Path(self.sne_c_path).glob("*_TC2.data"))[0]

        displacement_image = str(list(tc_paths_1_parse.glob("*displacement_??.img*"))[0])
        displacement_tif_name = str(self.product_output_directory / (self.pair_dates + '_displacement_' +
                                                                     self.polarisation + '.tif'))
        self.displacement_path_tif = displacement_tif_name

        # dem_image = str(list(tc_paths_1_parse.glob("*elevation.img*"))[0])
        # dem_name = str(self.product_output_directory / (self.pair_dates + '_elevation_' +
        #                                                 self.polarisation + '.tif'))
        # self.dem_name = dem_name

        # local_incidence_image = str(list(tc_paths_1_parse.glob("*localIncidenceAngle.img*"))[0])
        # local_incidence_tif_name = str(self.product_output_directory / (self.pair_dates +
        #                                                                 '_local_incidence_' + self.polarisation + '.tif'))
        # self.local_incidence_path_tif = local_incidence_tif_name

        coh_image = str(list(tc_paths_2_parse.glob("*coh_*.img"))[0])
        coh_tif_name = str(self.product_output_directory / (self.pair_dates +
                                                            '_coherence_' + self.polarisation + '.tif'))
        self.coherence_path_tif = coh_tif_name

        intensity_image = str(list(tc_paths_2_parse.glob("*Intensity_*.img"))[0])
        intensity_tif_name = str(self.product_output_directory / (self.pair_dates +
                                                                  '_intensity_' + self.polarisation + '.tif'))
        self.intensity_path_tif = intensity_tif_name

        intensity_image_db = str(list(tc_paths_2_parse.glob("*Intensity_*.img"))[1])
        intensity_tif_name_db = str(self.product_output_directory / (self.pair_dates +
                                                                     '_intensity_db_' + self.polarisation + '.tif'))
        self.intensity_db_path_tif = intensity_tif_name_db

        phase_image = str(list(tc_paths_2_parse.glob("Phase_*.img"))[0])
        phase_tif_name = str(self.product_output_directory / (self.pair_dates +
                                                              '_w_phase_' + self.polarisation + '.tif'))
        self.w_phase_path_tif = phase_tif_name

        unw_phase_image = str(list(tc_paths_2_parse.glob("*Unw_*.img"))[0])
        unw_phase_tif_name = str(self.product_output_directory / (self.pair_dates +
                                                                  '_unw_phase_' + self.polarisation + '.tif'))
        self.unw_phase_path_tif = unw_phase_tif_name

        displacement_2_tiff = self.gdal_exe + ' -of GTiff ' + displacement_image + ' ' + displacement_tif_name
        # dem_2_tiff = self.gdal_exe + ' -of GTiff ' + dem_image + ' ' + dem_name
        # local_incidence_2_tiff = self.gdal_exe + ' -of GTiff ' + local_incidence_image + ' ' + local_incidence_tif_name
        coh_image_2_tiff = self.gdal_exe + ' -of GTiff ' + coh_image + ' ' + coh_tif_name
        intensity_image_2_tiff = self.gdal_exe + ' -of GTiff ' + intensity_image + ' ' + intensity_tif_name
        intensity_db_2_tiff = self.gdal_exe + ' -of GTiff ' + intensity_image_db + ' ' + intensity_tif_name_db
        phase_image_2_tiff = self.gdal_exe + ' -of GTiff ' + phase_image + ' ' + phase_tif_name
        unw_phase_image_2_tiff = self.gdal_exe + ' -of GTiff ' + unw_phase_image + ' ' + unw_phase_tif_name

        images_2_tiff_list = [displacement_2_tiff, coh_image_2_tiff,
                              intensity_image_2_tiff, intensity_db_2_tiff, phase_image_2_tiff, unw_phase_image_2_tiff]
        for image in images_2_tiff_list:
            gen_fun.subprocess(image)
        print("converting img files to geotiff done")

    def clean_tif_borders(self):
        print("cleaning borders")
        rc.coherence_mask(self.coherence_path_tif, self.unw_phase_path_tif)
        rc.coherence_mask(self.coherence_path_tif, self.displacement_path_tif)
        print("cleaning borders done")

    def detrend(self):
        print("detrending tif")
        rc.detrend_function(self.unw_phase_path_tif)
        rc.detrend_function(self.unw_phase_path_tif)
        print("detrending done")

    def copy_aux_from_main(self):
        for aux_path in self.list_of_all_aux_caches_paths:
            if not aux_path.exists():
                shutil.copytree(self.pair_cache_directory_path_aux_primary, aux_path)

    def clean_useless_data(self):
        print("cleaning useless data")
        if self.sys_argv_int == 0:
            if not self.snaphu_logs_path.exists():
                self.snaphu_logs_path.mkdir()
        snaphu_conf = list(Path(self.sne_c_path).glob("snaphu.conf"))[0]
        snaphu_log = list(Path(self.sne_c_path).glob("snaphu.log"))[0]
        snaphu_conf_new_name = self.snaphu_logs_path / ("snaphu_" + self.pair_dates + ".conf")
        snaphu_log_new_name = self.snaphu_logs_path / ("snaphu_" + self.pair_dates + ".log")
        shutil.move(snaphu_conf, snaphu_conf_new_name)
        shutil.move(snaphu_log, snaphu_log_new_name)

        snaphu_directory = self.product_output_directory / 'Snaphu_Output_sne'
        data_list = self.product_output_directory.glob('*.data')
        dim_list = self.product_output_directory.glob('*.dim')
        for data_f, dim_f in zip(data_list, dim_list):
            try:
                os.remove(dim_f)
                shutil.rmtree(data_f)
            except OSError:
                pass
        try:
            shutil.rmtree(snaphu_directory)
            # shutil.rmtree(self.pair_cache_directory_path_java)
            # shutil.rmtree(self.pair_cache_directory_path_aux)
        except OSError:
            pass
        print("cleaning useless data done")

    def run_insar_processing(self):
        start_time = datetime.datetime.now()
        if self.sys_argv_int == 0:
            set_up_object = setup_code.set_up_for_insar_run(self.json_pairs_list_path, self.base_path, self.max_delta_t,
                                                            self.sys_argv_int)
            set_up_object.run_functions()
        else:
            pass
        self.num_of_looks()
        self.write_param_file()
        self.add_pairs_path()
        self.add_cache_paths_commands()
        self.applyOrbit()
        self.topsar_split()
        self.Back_Geocoding()
        self.esd()
        self.make_interferogram()
        self.topsar_deburst()
        self.merge_swaths()
        self.TopoPhaseRemoval()
        self.Multilook()
        self.Goldstein_Phase_Filtering()
        self.spatial_subset()
        self.Snaphu_Export()
        self.Snaphu_Unwrapping()
        self.Snaphu_Import()
        self.Phase_2_Displacement()
        self.Terrain_Correction()
        self.img_2_geotiff()
        self.clean_tif_borders()
        self.detrend()
        if self.sys_argv_int == 0:
            self.copy_aux_from_main()
        else:
            pass
        self.clean_useless_data()
        end_time = datetime.datetime.now()
        total_time = gen_fun.timestamp(end_time) - gen_fun.timestamp(start_time)
        self.time_to_run = total_time
        self.write_time_file()


if __name__ == '__main__':
    sys_index_var = sys.argv[1]
    run_insar_processing_object = insar_processing_snap(sys_index_var)
    run_insar_processing_object.run_insar_processing()
