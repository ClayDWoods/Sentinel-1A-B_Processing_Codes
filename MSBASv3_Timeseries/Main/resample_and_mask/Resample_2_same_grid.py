from Main.general_set_up_functions import rasterio_basic_functions
from statistics import mean
import math
import numpy as np


class resample_2_common_grid():
    def __init__(self, asc_dsc_list_displacement, asc_dsc_list_coh, base_path, coherence_threshold):
        self.asc_dsc_list_displacement = asc_dsc_list_displacement
        self.asc_dsc_list_coh = asc_dsc_list_coh
        self.base_path = base_path
        self.coherence_mask_path = self.base_path / ('coherence_mask')
        self.coh_thresh = coherence_threshold
        self.master_width = None
        self.master_height = None
        self.master_affine = None
        self.general_mask = None
        self.new_kwargs = None
        self.window = None
        self.window_kwargs = None

    def define_master_resample_parameters(self):
        longitude_spacing_list, upper_left_corner_longitudes_list = [], []
        lattitude_spacing_list, upper_left_corner_lattitude_list = [], []
        lower_right_corner_longitude_list = []
        lower_right_corner_lattitude_list = []

        for asc_dsc in self.asc_dsc_list_displacement:
            Affine_transform = rasterio_basic_functions.get_tif_transform(asc_dsc)
            tif_height, tif_width = (rasterio_basic_functions.get_tif_height(asc_dsc),
                                     rasterio_basic_functions.get_tif_width(asc_dsc))
            lower_right_corner_longitude = (Affine_transform * (tif_width, tif_height))[0]
            lower_right_corner_lattitude = (Affine_transform * (tif_width, tif_height))[1]
            longitude_spacing, upper_left_corner_longitude = Affine_transform.a, Affine_transform.c
            lattitude_spacing, upper_left_corner_lattitude = Affine_transform.e, Affine_transform.f
            upper_left_corner_longitudes_list.append(upper_left_corner_longitude)
            upper_left_corner_lattitude_list.append(upper_left_corner_lattitude)
            lower_right_corner_longitude_list.append(lower_right_corner_longitude)
            lower_right_corner_lattitude_list.append(lower_right_corner_lattitude)
            longitude_spacing_list.append(longitude_spacing)
            lattitude_spacing_list.append(lattitude_spacing)

        max_longitude = max(lower_right_corner_longitude_list) + 0.05
        min_longitude = min(upper_left_corner_longitudes_list) - 0.05
        max_lattitude = max(upper_left_corner_lattitude_list) + 0.05
        min_lattitude = min(lower_right_corner_lattitude_list) - 0.05

        mean_longitude_spacing = mean(longitude_spacing_list)
        mean_lattitude_spacing = mean(lattitude_spacing_list)

        self.master_width = math.ceil((abs(max_longitude - min_longitude))
                                      / abs(mean_longitude_spacing))
        self.master_height = math.ceil((abs(max_lattitude - min_lattitude))
                                       / abs(mean_lattitude_spacing))
        self.master_affine = rasterio_basic_functions.make_affine(mean_longitude_spacing,
                                                                  0.0,
                                                                  min_longitude,
                                                                  0.0,
                                                                  mean_lattitude_spacing,
                                                                  max_lattitude)

    def resample_displace_and_coh(self):
        mask_array_base = np.zeros((self.master_height, self.master_width))
        num_of_images = len(self.asc_dsc_list_displacement)
        for displacement_tif in self.asc_dsc_list_displacement:
            displacement_array_old = rasterio_basic_functions.tif_2_array(displacement_tif)
            displacement_transform_old = rasterio_basic_functions.get_tif_transform(displacement_tif)
            base_crs = rasterio_basic_functions.get_tif_crs(displacement_tif)
            displacement_reprojected_array = rasterio_basic_functions.reproject_tif_array(displacement_tif,
                                                                                          self.master_height,
                                                                                          self.master_width,
                                                                                          displacement_array_old,
                                                                                          displacement_transform_old,
                                                                                          base_crs,
                                                                                          self.master_affine)
            general_mask = np.where(displacement_reprojected_array != 0,
                                    1,
                                    0)
            mask_array_base = mask_array_base + general_mask
        self.general_mask = np.where(mask_array_base == num_of_images,
                                     1,
                                     0)
        coherence_mask_array = np.zeros((self.master_height, self.master_width))
        ii = 0
        displacement_tif_rcut_list_names = []
        for displacement_tif, coh_tif in zip(self.asc_dsc_list_displacement, self.asc_dsc_list_coh):
            displacement_array_old = rasterio_basic_functions.tif_2_array(displacement_tif)
            coh_array_old = rasterio_basic_functions.tif_2_array(coh_tif)
            displacement_transform_old = rasterio_basic_functions.get_tif_transform(displacement_tif)
            coh_transform_old = rasterio_basic_functions.get_tif_transform(coh_tif)
            base_crs = rasterio_basic_functions.get_tif_crs(displacement_tif)
            base_kwargs = rasterio_basic_functions.get_kwargs(displacement_tif)
            base_kwargs.update({'nodata': 0.,
                                'width': self.master_width,
                                'height': self.master_height,
                                'transform': self.master_affine
                                })
            if ii == 0:
                self.new_kwargs = base_kwargs
            ii = ii + 1

            displacement_tif_new_name = displacement_tif.parent / (displacement_tif.stem + '_rcut.tif')
            displacement_tif_rcut_list_names.append(displacement_tif_new_name)
            coherence_tif_new_name = coh_tif.parent / (coh_tif.stem + '_rcut.tif')
            displacement_reprojected_array = rasterio_basic_functions.reproject_tif_array(displacement_tif,
                                                                                          self.master_height,
                                                                                          self.master_width,
                                                                                          displacement_array_old,
                                                                                          displacement_transform_old,
                                                                                          base_crs,
                                                                                          self.master_affine)
            coherence_reprojected_array = rasterio_basic_functions.reproject_tif_array(coh_tif,
                                                                                       self.master_height,
                                                                                       self.master_width,
                                                                                       coh_array_old,
                                                                                       coh_transform_old,
                                                                                       base_crs,
                                                                                       self.master_affine)
            displacement_masked_array = np.multiply(displacement_reprojected_array, self.general_mask).astype('float32')
            coherence_masked_array = np.multiply(coherence_reprojected_array, self.general_mask).astype('float32')
            coherence_mask_array = coherence_mask_array + coherence_masked_array
            rasterio_basic_functions.write_reprojected_array_2_tif(displacement_masked_array,
                                                                   displacement_tif_new_name,
                                                                   base_kwargs)
            rasterio_basic_functions.write_reprojected_array_2_tif(coherence_masked_array,
                                                                   coherence_tif_new_name,
                                                                   base_kwargs)
        coherence_mask_array_use = (coherence_mask_array / num_of_images).astype('float32')
        if not self.coherence_mask_path.exists():
            self.coherence_mask_path.mkdir()
        coherence_mask_array_tif_name = self.coherence_mask_path / ('coherence_mask.tif')
        rasterio_basic_functions.write_reprojected_array_2_tif(coherence_mask_array_use,
                                                               coherence_mask_array_tif_name,
                                                               self.new_kwargs)
        window = rasterio_basic_functions.get_data_window_f(displacement_tif_rcut_list_names[0])
        kwargs_windowed = rasterio_basic_functions.windowed_kwargs(displacement_tif_rcut_list_names[0], window)
        self.window = window
        self.window_kwargs = kwargs_windowed
        coherence_mask_windowed_use = rasterio_basic_functions.tif_2_array_window(coherence_mask_array_tif_name,
                                                                                  self.window)
        rasterio_basic_functions.write_reprojected_array_2_tif(coherence_mask_windowed_use,
                                                               coherence_mask_array_tif_name,
                                                               self.window_kwargs)
        for displacement_tif_rcut in displacement_tif_rcut_list_names:
            displacement_tif_rcut_array = (rasterio_basic_functions.
                                           tif_2_array_window(displacement_tif_rcut, self.window))
            coherence_mask_apply = np.where(coherence_mask_windowed_use < self.coh_thresh,
                                            0,
                                            1)
            displacement_tif_rcut_coh_masked_name = (displacement_tif_rcut.parent
                                                     / (displacement_tif_rcut.stem + '_mcoh.tif'))
            displacement_tif_coh_masked_array = (np.multiply(displacement_tif_rcut_array,
                                                             coherence_mask_apply).astype('float32'))

            rasterio_basic_functions.write_reprojected_array_2_tif(displacement_tif_coh_masked_array,
                                                                   displacement_tif_rcut_coh_masked_name,
                                                                   self.window_kwargs)

    def run_resample(self):
        self.define_master_resample_parameters()
        self.resample_displace_and_coh()
