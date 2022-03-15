import sys
sys.path.append(".")

from Main.helper_functions import rasterio_basic_functions as rb
import numpy as np
from numpy.linalg import inv


def coherence_mask(coh_tif, tif_2_mask):
    coherence_array = rb.tif_2_array(coh_tif)
    coherence_mask_array = np.where(coherence_array > 0,
                                    1,
                                    0)
    array_masked = np.multiply(coherence_mask_array, rb.tif_2_array(tif_2_mask)).astype('float32')
    kwargs_tif = rb.get_kwargs(tif_2_mask)
    rb.write_reprojected_array_2_tif(array_masked, tif_2_mask, kwargs_tif)


def detrend_function(tif_2_detrend):
    tif_2_detrend_array = rb.tif_2_array(tif_2_detrend)
    kwargs = rb.get_kwargs(tif_2_detrend)

    # (ATA)-1ATB
    tif_width = rb.get_tif_width(tif_2_detrend)
    tif_height = rb.get_tif_height(tif_2_detrend)
    x_vector = np.array([i for i in range(0, tif_width, 1)])
    y_vector = np.array([i for i in range(0, tif_height, 1)])

    xx, yy = np.meshgrid(x_vector, y_vector)
    xx_flat = xx.flatten()
    yy_flat = yy.flatten()
    ones_ = np.ones(len(xx_flat))
    B = tif_2_detrend_array.flatten()

    A = np.column_stack((xx_flat, yy_flat, ones_))
    part_1 = inv(np.matmul(A.T, A))
    part_2 = np.matmul(part_1, A.T)
    plane_fit = np.matmul(part_2, B)
    z = xx * plane_fit[0] + yy * plane_fit[1] + plane_fit[2]
    tif_detrend_array = np.where(tif_2_detrend_array == 0.,
                                 0,
                                 tif_2_detrend_array - z).astype('float32')
    rb.write_reprojected_array_2_tif(tif_detrend_array, tif_2_detrend, kwargs)
