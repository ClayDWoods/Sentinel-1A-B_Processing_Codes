import os
import sys

sys.path.append(".")

import json
import pandas as pd
import numpy as np
import math


def open_txt_file(txt_file_name):
    f_open = open(txt_file_name)
    return f_open


def write_txt_file(fname, txt_str):
    if isinstance(txt_str, str):
        conv_2_string = str(txt_str)
        with open(fname, mode="w+") as text_file:
            text_file.write(conv_2_string)
    else:
        with open(fname, mode='w+') as myfile:
            myfile.write('\n'.join(txt_str))


def df_from_lists(list_1, list_2, list_3, list_4):
    df = pd.DataFrame(list(zip(list_1, list_2, list_3, list_4)))
    df.sort_values(by=2)
    return df


def save_df(df, txt_file_name):
    df.to_csv(txt_file_name, index=False, header=False, sep=' ')


def read_df(txt_file_name):
    df_open = pd.read_csv(txt_file_name)
    return df_open


def write_json_from_dict(dict_, json_save_name):
    with open(json_save_name, "w+") as outfile:
        json.dump(dict_, outfile, indent=4)


def open_json_file(json_file_path):
    with open(json_file_path) as json_file:
        slc_pairs_dict = json.load(json_file)
    return slc_pairs_dict


def write_array_2_npy(array, save_name):
    if save_name.exists():
        os.remove(save_name)
    np.save(save_name, array, allow_pickle=False)


def open_numpy_file(save_name):
    np_array = np.load(save_name)
    return np_array


def nearest_bound_z(bound_float, base):
    if bound_float < 0:
        new = base * math.floor(bound_float / base)
    else:
        new = base * math.ceil(bound_float / base)
    return new


def nearest_bound_coord(bound_float, base, u_d):
    if u_d == 'up':
        new = base * math.ceil(bound_float / base)
    else:
        new = base * math.floor(bound_float / base)
    return new


def set_bounds(mean_, std_, num_of_stds, base):
    lower_bound = mean_ - num_of_stds * std_
    upper_bound = mean_ + num_of_stds * std_
    lower_bound_u = nearest_bound_z(lower_bound, base)
    upper_bound_u = nearest_bound_z(upper_bound, base)
    return lower_bound_u, upper_bound_u


def lon_lat_2_col_row_using_transform(transform, lons, lats):
    a, c = transform[0], transform[2]
    e, f = transform[4], transform[5]
    cols = []
    rows = []
    for lon, lat in zip(lons, lats):
        col = int(abs((lon - c) / a))
        row = int(abs((lat - f) / e))
        cols.append(col)
        rows.append(row)
    return cols, rows
