from Main.general_set_up_functions import general_functions


def header_txt_file_creation_f(lambda_, l_, master_width, master_height,
                               asc_time, asc_heading, asc_incidence,
                               dsc_time, dsc_heading, dsc_incidence,
                               msbasv3_path,
                               c_size, c_flag_col, c_flag_row, num_of_ref,
                               ref_type):
    for L_ in l_:
        for lambda__ in lambda_:
            lambda_string = '{:.8f}'.format(lambda__)
            format_string = 'FORMAT=2, 2'
            file_size_string = 'FILE_SIZE=' + str(int(master_width)) + ', ' + str(int(master_height))
            window_size_string = 'WINDOW_SIZE=0, ' + str(int(master_width - 1)) + ', 0, ' + str(int(master_height - 1))
            base_str = ''
            if c_flag_col and c_flag_row is not None:
                if ref_type == 'gps':
                    base_str = base_str + str(c_flag_col) + ', ' + str(c_flag_row) + ', '
                    c_flag_string = ('C_FLAG=' + str(num_of_ref) + ', ' + base_str + str(c_size) + ', ' + str(c_size))
                elif ref_type == 'auto':
                    for col, row in zip(c_flag_col, c_flag_row):
                        base_str = base_str + str(col) + ', ' + str(row) + ', '
                    if len(c_flag_col) != num_of_ref:
                        num_of_ref = len(c_flag_col)
                    c_flag_string = ('C_FLAG=' + str(num_of_ref) + ', ' + base_str + str(c_size) + ', ' + str(c_size))
                else:
                    base_str = base_str + str(c_flag_col) + ', ' + str(c_flag_row) + ', '
                    c_flag_string = ('C_FLAG=' + str(num_of_ref) + ', ' + base_str + str(c_size) + ', ' + str(c_size))
            else:
                c_flag_string = str('C_FLAG=0')
            r_flag_string = 'R_FLAG=' + str(L_) + ', ' + lambda_string
            i_flag_string = 'I_FLAG=0'
            set_asc_string = 'SET=' + str(asc_time) + ', ' + str(asc_heading) + ', ' + str(asc_incidence) + ', ' + \
                             ' asc.txt'
            set_dsc_string = 'SET=' + str(dsc_time) + ', ' + str(dsc_heading) + ', ' + str(dsc_incidence) + ', ' + \
                             ' dsc.txt'
            header_save_name = "header_" + str(L_) + "_" + lambda_string.split(".")[0] + "_" + \
                               lambda_string.split(".")[1] + '.txt'
            header_list = [format_string, file_size_string, window_size_string, c_flag_string, r_flag_string,
                           i_flag_string, set_asc_string, set_dsc_string]
            header_txt_save_name = msbasv3_path / header_save_name
            general_functions.write_txt_file(header_txt_save_name, header_list)
