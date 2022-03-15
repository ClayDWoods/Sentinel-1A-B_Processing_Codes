import sys

import numpy as np

sys.path.append(".")
import pandas as pd


def filter_months_years_regression(years_1, months_1, years_2, months_2,
                                   dates_yyyy_mm_list, dates_yyyy_yyy_list, msbas_np_names, regression_mode):
    df = pd.DataFrame(list(zip(dates_yyyy_mm_list, dates_yyyy_yyy_list, msbas_np_names)))
    str_filter = []
    if regression_mode == 'months':
        for year in years_1:
            for month in months_1:
                str_filter_u = year + month
                str_filter.append(str_filter_u)
        for str_filter_use in str_filter:
            df = df[~df[0].str.contains(str_filter_use)]
    else:
        for year in years_2:
            for month in months_2:
                str_filter_u = year + month
                str_filter.append(str_filter_u)
        for str_filter_use in str_filter:
            df = df[~df[0].str.contains(str_filter_use)]
    dates_yyyy_mm_list_updated = list(df[0])
    dates_yyyy_yyy_list_updated = list(df[1])
    msbas_np_names_updated = list(df[2])
    return dates_yyyy_mm_list_updated, dates_yyyy_yyy_list_updated, msbas_np_names_updated


def make_regressions_parallel(np_array, time_array):
    reg_coeff = np.polyfit(time_array, np_array, 1)
    slope_ = reg_coeff[0]
    return slope_
