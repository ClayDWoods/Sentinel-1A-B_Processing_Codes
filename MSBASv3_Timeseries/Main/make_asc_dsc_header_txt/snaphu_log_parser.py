import sys
sys.path.append(".")

from Main.general_set_up_functions import general_functions

def snaphu_log_parser_f(snaphu_txt_file):
    snaphu_log_open = general_functions.open_txt_file(snaphu_txt_file)
    base = 'BASELINE '
    ans = []
    for line in snaphu_log_open:
        line = line.strip()
        if line.startswith(base):
            ans.append(line)
    base_line_length = float(ans[0].split(' ')[1])
    snaphu_log_open.close()
    return base_line_length