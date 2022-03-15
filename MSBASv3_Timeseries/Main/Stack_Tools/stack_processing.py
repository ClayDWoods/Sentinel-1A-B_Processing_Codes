import sys
sys.path.append(".")

from Main.general_set_up_functions import rasterio_basic_functions, general_functions
from Main.Stack_Tools import stack_base
from pathlib import Path


class make_stack_asc_dsc():
    def __init__(self, stack_path, asc_list_mcoh, dsc_list_mcoh, asc_incidence,
                 dsc_incidence, c_flag_col, c_flag_row, c_size, coherence_mask_tif, num_of_ref, coh_mask_val):
        self.stack_path = stack_path
        self.stacked_txt_name_asc = self.stack_path / ('asc_stacked.txt')
        self.stacked_txt_name_dsc = self.stack_path / ('dsc_stacked.txt')
        self.coherence_mask_tif_path = coherence_mask_tif
        self.num_of_ref = num_of_ref
        self.coherence_threshold = coh_mask_val
        self.stack_path_name = None
        temp_asc_list = [str(x) for x in asc_list_mcoh]
        temp_dsc_list = [str(x) for x in dsc_list_mcoh]
        self.asc_list_mcoh = [Path(x) for x in temp_asc_list]
        self.dsc_list_mcoh = [Path(x) for x in temp_dsc_list]
        self.asc_list_mcoh_ = [Path(x) for x in temp_asc_list]
        self.dsc_list_mcoh_ = [Path(x) for x in temp_dsc_list]
        self.asc_list_mcoh__ = [Path(x) for x in temp_asc_list]
        self.dsc_list_mcoh__ = [Path(x) for x in temp_dsc_list]
        self.asc_list_mcoh___ = [Path(x) for x in temp_asc_list]
        self.dsc_list_mcoh___ = [Path(x) for x in temp_dsc_list]
        self.c_size = c_size
        self.asc_list_mcoh_filtered = None
        self.dsc_list_mcoh_filtered = None
        self.asc_stack_array = None
        self.dsc_stack_array = None
        self.stack_array = None
        self.asc_incidence = asc_incidence
        self.dsc_incidence = dsc_incidence
        self.base_kwargs = None
        self.c_flag_col = c_flag_col
        self.c_flag_row = c_flag_row

    def make_dirs(self):
        if not self.stack_path.exists():
            self.stack_path.mkdir()

    def make_new_lists(self):
        self.asc_list_mcoh_filtered = stack_base.filter_asc_dsc_list(self.asc_list_mcoh, self.asc_list_mcoh_)
        self.dsc_list_mcoh_filtered = stack_base.filter_asc_dsc_list(self.dsc_list_mcoh, self.dsc_list_mcoh_)
        asc_list_mcoh_filtered = stack_base.filter_asc_dsc_list(self.asc_list_mcoh__, self.asc_list_mcoh___)
        dsc_list_mcoh_filtered = stack_base.filter_asc_dsc_list(self.dsc_list_mcoh__, self.dsc_list_mcoh___)
        asc_filtered_write = [str(x) for x in asc_list_mcoh_filtered]
        dsc_filtered_write = [str(x) for x in dsc_list_mcoh_filtered]
        general_functions.write_txt_file(self.stacked_txt_name_asc, asc_filtered_write)
        general_functions.write_txt_file(self.stacked_txt_name_dsc, dsc_filtered_write)
        base_dir = self.asc_list_mcoh_filtered[0]
        self.base_kwargs = rasterio_basic_functions.get_kwargs(base_dir)

    def make_stack_arrays(self):
        self.asc_stack_array = stack_base.make_stack_base_asc_dsc(self.asc_list_mcoh_filtered, self.asc_incidence,
                                                                  self.c_flag_col, self.c_flag_row,
                                                                  self.c_size,
                                                                  self.coherence_mask_tif_path,
                                                                  self.num_of_ref,
                                                                  self.coherence_threshold)
        self.dsc_stack_array = stack_base.make_stack_base_asc_dsc(self.dsc_list_mcoh_filtered, self.dsc_incidence,
                                                                  self.c_flag_col, self.c_flag_row,
                                                                  self.c_size,
                                                                  self.coherence_mask_tif_path,
                                                                  self.num_of_ref,
                                                                  self.coherence_threshold)
        self.stack_array = ((self.asc_stack_array + self.dsc_stack_array)/2).astype('float32')

    def write_stack_array(self):
        self.stack_path_name = self.stack_path / ('Stack_Rate.tif')
        rasterio_basic_functions.write_reprojected_array_2_tif(self.stack_array, self.stack_path_name, self.base_kwargs)

    def run_functions(self):
        self.make_dirs()
        self.make_new_lists()
        self.make_stack_arrays()
        self.write_stack_array()













