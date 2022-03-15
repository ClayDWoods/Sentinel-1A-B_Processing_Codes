import sys, os
sys.path.append(".")

import pathlib
import shutil

def rename_msbas_log_f(msbas_log_paths):
    for msbas_log_path in msbas_log_paths:
        msbas_log_path_base = msbas_log_path.parent
        header_name = msbas_log_path_base.name
        new_msbas_name = msbas_log_path_base / (header_name + '_MSBAS_LOG.txt')
        shutil.copy2(msbas_log_path, new_msbas_name)

