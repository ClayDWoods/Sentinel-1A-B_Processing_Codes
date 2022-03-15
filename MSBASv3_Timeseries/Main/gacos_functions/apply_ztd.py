import sys

sys.path.append(".")

import numpy as np
from pathlib import Path
from numpy.linalg import inv
from Main.general_set_up_functions import rasterio_basic_functions as rb


def detrend_function(array_2_detrend):
    # (ATA)-1ATB
    tif_width = array_2_detrend.shape[1]
    tif_height = array_2_detrend.shape[0]
    x_vector = np.array([i for i in range(0, tif_width, 1)])
    y_vector = np.array([i for i in range(0, tif_height, 1)])

    xx, yy = np.meshgrid(x_vector, y_vector)
    xx_flat = xx.flatten()
    yy_flat = yy.flatten()
    ones_ = np.ones(len(xx_flat))
    b = array_2_detrend.flatten()

    a = np.column_stack((xx_flat, yy_flat, ones_))
    part_1 = inv(np.matmul(a.T, a))
    part_2 = np.matmul(part_1, a.T)
    plane_fit = np.matmul(part_2, b)
    z = xx * plane_fit[0] + yy * plane_fit[1] + plane_fit[2]
    tif_detrend_array = np.where(array_2_detrend == 0.,
                                 0,
                                 array_2_detrend - z).astype('float32')
    return tif_detrend_array

def apply_gacos_correction(base_rcut_mcoh_displacement_path, incidence_file_path, base_ztd_path):
    base_path_2_displacement = base_rcut_mcoh_displacement_path.absolute().parent
    base_kwargs = rb.get_kwargs(base_rcut_mcoh_displacement_path)
    incidence_array = np.deg2rad(rb.tif_2_array(incidence_file_path))
    mask_4_dtzd = np.where(incidence_array != 0,
                           1,
                           0)
    displacement_array = rb.tif_2_array(base_rcut_mcoh_displacement_path)
    rcut_mcoh_date_1 = base_rcut_mcoh_displacement_path.name.split('_')[0]
    rcut_mcoh_date_2 = base_rcut_mcoh_displacement_path.name.split('_')[1]
    dtzd_save_name = base_path_2_displacement / (rcut_mcoh_date_1 + '_' + rcut_mcoh_date_2 + '_dtzd.tif')
    base_gacos_correct_name = (base_path_2_displacement / (base_rcut_mcoh_displacement_path.stem +
                                                           '_gacos.tif'))
    gacos_1 = [Path(str(x)) for x in base_ztd_path.glob(rcut_mcoh_date_1 + '.ztd.tif')][0]
    gacos_2 = [Path(str(x)) for x in base_ztd_path.glob(rcut_mcoh_date_2 + '.ztd.tif')][0]
    gacos_array_1_reprojected = rb.reproject_tif_array_ztd(base_rcut_mcoh_displacement_path, gacos_1)
    gacos_array_2_reprojected = rb.reproject_tif_array_ztd(base_rcut_mcoh_displacement_path, gacos_2)
    gacos_dtzd = np.multiply((gacos_array_2_reprojected - gacos_array_1_reprojected), (1 / np.cos(incidence_array)))
    gacos_dtzd_masked = np.multiply(gacos_dtzd, mask_4_dtzd).astype('float32')
    corrected_array = displacement_array + gacos_dtzd_masked
    corrected_detrended_array = detrend_function(corrected_array)
    rb.write_reprojected_array_2_tif(gacos_dtzd_masked, dtzd_save_name, base_kwargs)
    rb.write_reprojected_array_2_tif(corrected_detrended_array, base_gacos_correct_name, base_kwargs)
