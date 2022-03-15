import sys

sys.path.append(".")
from Main.general_set_up_functions import rasterio_basic_functions as rb
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from scipy import interpolate


def incidence_tif_file_creation(base_meta_path, base_file_rcut, base_file_rcut_mcoh, save_path, dsc_asc):
    inc_save_1 = save_path / ('incidence_' + dsc_asc + '.tif')
    inc_save_2 = save_path / ('incidence_rcut_' + dsc_asc + '.tif')
    inc_save_3 = save_path / ('incidence_rcut_mcoh_' + dsc_asc + '.tif')
    geo_grid_point = 'geolocationGridPoint'
    incidence = 'incidenceAngle'
    pixel_col = 'pixel'
    line_row = 'line'
    number_of_samples = 'numberOfSamples'
    number_of_lines = 'numberOfLines'
    base_kwargs_rcut = rb.get_kwargs(base_file_rcut)
    base_width = rb.get_tif_width(base_file_rcut)
    base_height = rb.get_tif_height(base_file_rcut)
    meta_data_paths = [str(x) for x in base_meta_path.glob('s1a*xml')]
    meta_data_paths.sort()
    base_date = base_file_rcut_mcoh.name.split('_')[0]
    iws = [('iw1-slc-vv-' + base_date), ('iw2-slc-vv-' + base_date), ('iw3-slc-vv-' + base_date)]
    base_meta_list = []
    for iw in iws:
        match = [s for s in meta_data_paths if iw in s][0]
        base_meta_list.append(Path(match))

    # xml parse
    max_incidence_angles = []
    min_incidence_angles = []
    for swath_xml_file in base_meta_list:
        tree_xml = ET.parse(swath_xml_file)
        root_xml = tree_xml.getroot()

        list_inc = []
        list_col = []
        list_row = []
        for geo_grid_ in root_xml.iter(geo_grid_point):
            incidence_angle = float(geo_grid_.find(incidence).text)
            col_ = int(geo_grid_.find(pixel_col).text)
            row_ = int(geo_grid_.find(line_row).text)
            list_inc.append(incidence_angle)
            list_col.append(col_)
            list_row.append(row_)
        for num_sample in root_xml.iter(number_of_samples):
            num_samples = int(num_sample.text)
        for num_line in root_xml.iter(number_of_lines):
            num_lines = int(num_line.text)
        del tree_xml, root_xml
        max_incidence = np.max(list_inc)
        min_incidence = np.min(list_inc)
        max_incidence_angles.append(max_incidence)
        min_incidence_angles.append(min_incidence)

    max_overall_incidence_angle = np.max(max_incidence_angles)
    min_overall_incidence_angle = np.min(min_incidence_angles)
    cols_all = np.arange(0, base_width)

    if dsc_asc == 'asc':
        f = interpolate.interp1d([0, base_width - 1], [min_overall_incidence_angle, max_overall_incidence_angle])
    else:
        f = interpolate.interp1d([0, base_width - 1], [max_overall_incidence_angle, min_overall_incidence_angle])
    f_int = list(f(cols_all))
    final_int_array = np.tile(f_int, (base_height, 1)).astype('float32')
    rb.write_reprojected_array_2_tif(final_int_array, inc_save_1, base_kwargs_rcut)

    incidence_array = rb.tif_2_array(inc_save_1)
    rcut_array = rb.tif_2_array(base_file_rcut)
    rcut_kwargs_inc = rb.get_kwargs(inc_save_1)
    rcut_mask_array = np.where(rcut_array != 0,
                               1,
                               0)
    incidence_cut_array = np.multiply(incidence_array, rcut_mask_array).astype('float32')
    rb.write_reprojected_array_2_tif(incidence_cut_array, inc_save_2, rcut_kwargs_inc)

    window = rb.get_data_window_f(inc_save_2)
    kwargs_windowed = rb.windowed_kwargs(inc_save_2, window)
    coherence_mask_windowed_use = rb.tif_2_array_window(inc_save_2, window)
    mcoh_array = rb.tif_2_array(base_file_rcut_mcoh)
    mcoh_mask_array = np.where(mcoh_array != 0,
                               1,
                               0)
    project_array_2_mcoh_incidence_mask = np.multiply(coherence_mask_windowed_use, mcoh_mask_array).astype('float32')
    rb.write_reprojected_array_2_tif(project_array_2_mcoh_incidence_mask, inc_save_3, kwargs_windowed)
