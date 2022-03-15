import sys

sys.path.append(".")

import json
import pandas as pd


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


def write_txt_file_gacos(fname, txt_str):
    with open(fname, mode='w+') as myfile:
        for lines in txt_str:
            myfile.write('\n'.join(str(line) for line in lines))
            myfile.write('\n')


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


def HHMMSS_2_HHMM(time_2_convert):
    hour_ = time_2_convert[0:2]
    minute_ = time_2_convert[2:4]
    seconds_ = time_2_convert[4:]
    seconds_2_minutes = round((int(seconds_) / 60))
    minutes = str(int(minute_) + seconds_2_minutes)
    return hour_, minutes


def unique_dates(processed_files_list):
    dates_list = []
    for d in processed_files_list:
        split_name_path = d.name.split('_')
        date_1, date_2 = split_name_path[0], split_name_path[1]
        dates_list.append(date_1)
        dates_list.append(date_2)
    dates_list_un = list(set(dates_list))
    dates_list_un.sort()
    return dates_list_un
