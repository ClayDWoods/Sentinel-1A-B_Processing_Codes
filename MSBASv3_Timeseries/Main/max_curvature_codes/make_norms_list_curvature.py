import sys, os
sys.path.append(".")

import numpy as np
from scipy.interpolate import CubicSpline


def make_norms_list_f(msbas_log_list, list_factor):
    list_factor_int = int(list_factor)
    log_list_norms_x = []
    log_list_norms_axy = []
    lam_lis = []
    string_norm = 'computed ||x|| and ||Ax-Y|| norms:'
    for header in msbas_log_list:
        if header is not None:
            lines = [line.rstrip('\n') for line in open(header)]
        else:
            lines = None
        if lines is not None:
            for log1 in lines:
                if log1 is not None:
                    if log1.startswith(string_norm):
                        norms_string_1 = log1
                        x_norm_1 = float(norms_string_1.split(" ")[-2])
                        axy_norm_1 = float(norms_string_1.split(" ")[-1])
                        log_list_norms_x.append(x_norm_1)
                        log_list_norms_axy.append(axy_norm_1)
        t_ = '.'.join(header.name.split('_')[2:4])
        lam_lis.append(t_)
    lam_axy = lam_lis.copy()
    log_array_norms_x = np.array(log_list_norms_x)
    log_array_norms_axy = np.array(log_list_norms_axy)
    lam_lis, log_array_norms_x = zip(*sorted(zip(lam_lis, log_array_norms_x)))
    lam_axy, log_array_norms_axy = zip(*sorted(zip(lam_axy, log_array_norms_axy)))
    lam_lis = lam_lis[list_factor_int:]
    log_array_norms_x = log_array_norms_x[list_factor_int:]
    lam_axy = lam_axy[list_factor_int:]
    log_array_norms_axy = log_array_norms_axy[list_factor_int:]
    return log_array_norms_x, log_array_norms_axy, lam_lis, lam_axy


def find_curvatures(log_array_norms_x, log_array_norms_axy, lam_lis, lam_axy):
    cs_2 = CubicSpline(np.array(lam_lis).astype('float64'),
                       np.log10(np.square(np.array(log_array_norms_x)).astype('float64')))
    cs_1 = CubicSpline(np.array(lam_axy).astype('float64'),
                       np.log10(np.square(np.array(log_array_norms_axy).astype('float64'))))

    # first order derivative
    dcs_11_eval = cs_1(lam_lis, 1)
    dcs_12_eval = cs_1(lam_lis, 2)
    # second order derivative
    dcs_21_eval = cs_2(lam_axy, 1)
    dcs_22_eval = cs_2(lam_axy, 2)

    curv = []
    for dr1_11, dr1_12, dr1_21, dr1_22 in zip(dcs_11_eval, dcs_12_eval, dcs_21_eval, dcs_22_eval):
        k = 2*(((dr1_21 * dr1_12) - (dr1_22 * dr1_11)) / ((dr1_21 ** 2 + dr1_11 ** 2) ** (3 / 2)))
        curv.append(k)
    curv = np.array(curv)
    return curv


def find_max_curv_lambda(curv, lam_list):
    curv_max_index = np.argmax(curv)
    lam_max_value = np.array(lam_list)[curv_max_index]
    return lam_max_value
