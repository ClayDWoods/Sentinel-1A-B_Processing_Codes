import sys
sys.path.append(".")

from Main.helper_functions import general_functions
from pathlib import Path
import pathlib
import zipfile
import shutil
import datetime
from datetime import timedelta
from datetime import datetime


class set_up_for_insar_run():
    def __init__(self, json_file_path, base_path, max_delta_t, sys_arg):
        self.sys_arg = sys_arg
        self.json_file_path = json_file_path
        self.base_path = base_path
        self.data_path = self.base_path / 'Main/Data'
        self.max_delta_t = int(max_delta_t)
        self.slc_pairs_dict = None
        self.kml_frame_coords_kml_path = None

    def make_pair_list(self):
        slc_list = self.data_path.glob('*.zip')
        slc_list__ = list(self.data_path.glob('*.zip'))
        slc_date_list = []
        max_delta_t_ob = timedelta(self.max_delta_t)
        for SLC in slc_list:
            slc_date = general_functions.slc_path_2_date(SLC)
            slc_date_list.append(slc_date)
        slc_list_ = [x for y, x in sorted(zip(slc_date_list, slc_list__))]
        slc_date_list_sorted = sorted(slc_date_list)
        slc_len = len(slc_date_list_sorted)
        j, k = 0, 0
        slc_pair_dict = {}
        for slc_date in slc_date_list_sorted:
            for i in range(j, slc_len):
                slc_pair_iset_dict = {}
                if i + 2 <= slc_len:
                    t1, t2 = (datetime.strptime(slc_date, '%Y%m%d'),
                              datetime.strptime(slc_date_list_sorted[i + 1], '%Y%m%d'))
                    t_delta = t2 - t1
                    if t_delta <= max_delta_t_ob:
                        pair_name = 'Pair' + '_' + str(k)
                        slc_list_temp = [str(slc_list_[j]), str(slc_list_[i + 1])]
                        slc_pair_dict[pair_name] = slc_list_temp
                        slc_pair_iset_dict[pair_name] = slc_list_temp
                        pairs_dir = self.data_path.parent / "Pairs/"
                        if not pairs_dir.exists():
                            pairs_dir.mkdir()
                        slc_pair_iset_name = self.data_path.parent / ("Pairs/" + pair_name + ".json")
                        general_functions.write_json_from_dict(slc_pair_iset_dict, slc_pair_iset_name)
                        k = k + 1
            j = j + 1
        json_save_name = self.data_path / "Pairs_List.json"
        general_functions.write_json_from_dict(slc_pair_dict, json_save_name)
        self.slc_pairs_dict = general_functions.open_json_file(self.json_file_path)

    def create_pair_directories(self):
        jj = 0
        for slc_list in self.slc_pairs_dict.values():
            slc_1 = general_functions.slc_path_2_date(slc_list[0])
            slc_2 = general_functions.slc_path_2_date(slc_list[1])
            path_processed_data = self.base_path / 'Main/Processed_Data'
            path_cache_directory = self.base_path / 'Main/Cache'
            pair_directory_name = slc_1 + '_' + slc_2
            pair_cache_directory_name = pair_directory_name + '_cache'
            pair_cache_directory_path = self.base_path / ('Main/Cache/' + pair_cache_directory_name)
            pair_cache_directory_path_java = self.base_path / ('Main/Cache/' + pair_cache_directory_name + '/java_tmp')
            self.kml_frame_coords_kml_path = self.data_path / 'kml_frame'
            pair_directory_path = path_processed_data / pair_directory_name
            if not path_processed_data.exists():
                path_processed_data.mkdir()
            if not path_cache_directory.exists():
                path_cache_directory.mkdir()
            if not pair_cache_directory_path.exists():
                pair_cache_directory_path.mkdir()
            if not pair_cache_directory_path_java.exists():
                pair_cache_directory_path_java.mkdir()
            if not self.kml_frame_coords_kml_path.exists():
                self.kml_frame_coords_kml_path.mkdir()
            if not pair_directory_path.exists():
                pair_directory_path.mkdir()
            if jj == 0 and self.sys_arg == 0:
                pair_cache_directory_path_aux = (self.base_path /
                                                 ('Main/Cache/' + pair_cache_directory_name + '/.esa_snap'))
                if not pair_cache_directory_path_aux.exists():
                    pair_cache_directory_path_aux.mkdir()
                jj = jj + 1

    def extract_meta_data(self):
        meta_data_directory = self.data_path / 'Meta_Data'
        if not meta_data_directory.exists():
            meta_data_directory.mkdir()
        s1_data_dir_names = self.data_path.glob('S1*.zip')
        ii = 0
        for zip_safe_dir in s1_data_dir_names:
            with zipfile.ZipFile(zip_safe_dir, "r") as zf:
                if ii == 0:
                    preview_path = [s for s in zf.namelist() if 'preview' and '.kml' in s][0]
                    parent_direct_kml = self.kml_frame_coords_kml_path / Path(preview_path).parent.parent
                    zf.extract(preview_path, path=self.kml_frame_coords_kml_path)
                    path_to_extracted = self.kml_frame_coords_kml_path / preview_path
                    shutil.copy2(path_to_extracted, self.kml_frame_coords_kml_path)
                    shutil.rmtree(parent_direct_kml)
                ii = ii + 1

                list_annotation_xml = \
                    [s for s in zf.namelist() if "annotation" and "xml" in s and "annotation" and "/s1a-" in s]
                paths_list = []
                parent_direct_meta = meta_data_directory / Path(list_annotation_xml[0]).parent.parent
                for annotation_xml_path in list_annotation_xml:
                    zf.extract(annotation_xml_path, path=meta_data_directory)
                    path1 = meta_data_directory / annotation_xml_path
                    paths_list.append(path1)
                for item in paths_list:
                    item_path = Path(item)
                    shutil.copy2(item_path, meta_data_directory)
                shutil.rmtree(parent_direct_meta)
            zf.close()

    def run_functions(self):
        self.make_pair_list()
        self.create_pair_directories()
        self.extract_meta_data()



