from pathlib import Path
import json
from subprocess import Popen, PIPE, STDOUT
import time


def slc_path_2_date(slc_path):
    if isinstance(slc_path, str):
        slc_name_from_path = Path(slc_path).name
        slc_date = slc_name_from_path.split('_')[5].split('T')[0]
    else:
        slc_name_from_path = slc_path.name
        slc_date = slc_name_from_path.split('_')[5].split('T')[0]
    return slc_date


def write_json_from_dict(dict, json_save_name):
    with open(json_save_name, "w+") as outfile:
        json.dump(dict, outfile, indent=4)


def timestamp(date):
    return time.mktime(date.timetuple())


def write_txt_file(fname, txt_str):
    conv_2_string = str(txt_str)
    with open(fname, "w+") as text_file:
        text_file.write(conv_2_string)


def open_json_file(json_file_path):
    with open(json_file_path) as json_file:
        slc_pairs_dict = json.load(json_file)
    return slc_pairs_dict


def subprocess(command):
    p = Popen(command, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT)
    output = p.stdout.read()
    print(output)


def make_dictionary_from_params(param_keys_all, param_vals_all):
    param_dict = {}
    for param_names_set, param_vals in zip(param_keys_all, param_vals_all):
        for param_name, param_val in zip(param_names_set, param_vals):
            param_dict[param_name] = param_val
    return param_dict


def get_sec(time_str):
    h, m, s = time_str[0:2], time_str[2:4], time_str[4:6]
    return int(h) * 3600 + int(m) * 60 + float(s)


def secondsToTime(seconds):
    minutes, seconds = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    time_hh_mm_ss = ("%02d:%02d:%02d" % (hours, minutes, seconds))
    hh, mm, ss = time_hh_mm_ss.split(':')
    time_string = hh + mm + ss
    return time_string
