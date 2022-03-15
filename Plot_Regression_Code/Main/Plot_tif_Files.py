import sys, os

sys.path.append("..")
os.environ['PROJ_LIB'] = '/projects/clwo4142/.conda_pkgss/proj4-4.9.3-hc8507d1_7/share/proj'

from general_functions_all import rasterio_basic_functions, general_functions
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt

plt.style.use('seaborn-whitegrid')
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
from matplotlib import cm, colors
import numpy as np


class plot_tifs_png():
    def __init__(self, linear_or_coh, tif_name):
        self.base_path = Path.cwd().absolute().parent
        self.path_to_base = self.base_path / 'Plots_Regression'
        self.path_to_png = self.path_to_base / 'Linear_Rate_Spatial_Plots'
        self.linear_or_coh = linear_or_coh
        self.num_of_std = 2
        self.base_linear = 5
        self.base_coh = 0.1
        if self.linear_or_coh == 'linear':
            self.paths_to_use = self.path_to_png / tif_name
            tif_name_split = tif_name.split('.')[0]
            self.path_save_name = self.path_to_png / (tif_name_split + '.png')
            self.path_save_name_hist = self.path_to_png / (tif_name_split + '_hist.png')
        else:
            self.path_to_coh = self.path_to_png
            self.paths_to_use = self.path_to_coh / tif_name
            tif_name_split = tif_name.split('.')[0]
            self.path_save_name = self.path_to_png / (tif_name_split + '.png')
            self.path_save_name_hist = self.path_to_png / (tif_name_split + '_hist.png')
        self.save_name_title = tif_name_split

    def make_plots(self):
        base_array = rasterio_basic_functions.tif_2_array(self.paths_to_use)
        base_array[base_array == 0.] = np.nan
        if self.linear_or_coh == 'linear':
            base_array = base_array * 1000
            mean_tif = np.nanmean(base_array)
            std_tif = np.nanstd(base_array)
            cmap = mpl.cm.rainbow
            lower_bound_u_z, upper_bound_u_z = general_functions.set_bounds(mean_tif,
                                                                            std_tif,
                                                                            self.num_of_std,
                                                                            self.base_linear)
        else:
            cmap = mpl.cm.gray
            min_coh = np.nanmin(base_array)
            max_coh = np.nanmax(base_array)
            lower_bound_u_z, upper_bound_u_z = (general_functions.nearest_bound_coord(min_coh, self.base_coh, 'down'),
                                                general_functions.nearest_bound_coord(max_coh, self.base_coh, 'up'))
        norm = mpl.colors.Normalize(vmin=lower_bound_u_z, vmax=upper_bound_u_z)
        left, bottom, right, top = rasterio_basic_functions.return_bounds(self.paths_to_use)
        lon_spacing = round(abs(right - left)/5, 1)
        lat_spacing = round(abs(top - bottom)/5, 1)
        if 0.25 < lon_spacing < 0.5:
            lon_spacing = 0.25
        elif 0.5 < lon_spacing < 1:
            lon_spacing = 0.5
        else:
            lon_spacing = 0.1
        if 0.25 < lat_spacing < 0.5:
            lat_spacing = 0.25
        elif 0.5 < lat_spacing < 1:
            lat_spacing = 0.5
        else:
            lat_spacing = 0.1
        base_degree_lon = lon_spacing
        base_degree_lat = lat_spacing
        fig, axs1 = plt.subplots(1, 1, figsize=(20, 22), facecolor='w', edgecolor='k')
        left_use, bottom_use, right_use, top_use = (general_functions.nearest_bound_coord(left,
                                                                                          base_degree_lat, 'down'),
                                                    general_functions.nearest_bound_coord(bottom,
                                                                                          base_degree_lon, 'down'),
                                                    general_functions.nearest_bound_coord(right,
                                                                                          base_degree_lat, 'up'),
                                                    general_functions.nearest_bound_coord(top,
                                                                                          base_degree_lon, 'up'))
        lon_mid = (left + right) / 2
        m = Basemap(epsg=4326, llcrnrlat=bottom, urcrnrlat=top,
                    llcrnrlon=left, urcrnrlon=right, resolution='f', lon_0=lon_mid, ax=axs1)
        m.arcgisimage(service='ESRI_Imagery_World_2D', xpixels=2500, verbose=True, dpi=300)
        im1 = m.imshow(base_array, vmin=lower_bound_u_z, vmax=upper_bound_u_z, cmap=cmap, norm=norm, origin="upper")
        parallels = np.arange(bottom_use, top_use + 1, base_degree_lon)
        meridians = np.arange(left_use, right_use + 1, base_degree_lat)
        m.drawparallels(parallels, labels=[True, False, True, False], linewidth=3.0, fontsize=14)
        m.drawmeridians(meridians, labels=[True, False, True, False], linewidth=3.0, fontsize=14)
        divider1 = make_axes_locatable(axs1)
        cax1 = divider1.append_axes("right", size="5%", pad=0.05)
        if self.linear_or_coh == 'linear':
            cb_ticks = list(np.arange(lower_bound_u_z, upper_bound_u_z + 1, self.base_linear))
        else:
            cb_ticks = list(np.arange(lower_bound_u_z, upper_bound_u_z + 0.1, self.base_coh))
        cb = fig.colorbar(im1, ax=axs1, cax=cax1, extend='both', ticks=cb_ticks)
        if self.linear_or_coh == 'linear':
            cb.set_label('mm/yr')
        else:
            cb.set_label('coherence')
        axs1.set_title(self.save_name_title)
        fig.tight_layout()
        if self.path_save_name.exists():
            os.remove(self.path_save_name)
        fig.savefig(self.path_save_name, dpi=600, facecolor='w', edgecolor='k',
                    orientation='portrait', bbox_inches='tight',
                    pad_inches=0.3)

    def make_histograms(self):
        base_array = rasterio_basic_functions.tif_2_array(self.paths_to_use)
        base_array[base_array == 0.] = np.nan
        if self.linear_or_coh == 'linear':
            base_array = base_array * 1000
        h_a = base_array.flatten()[~np.isnan(base_array.flatten())]
        h, b = np.histogram(h_a, bins=100, density=True)
        fig, axs = plt.subplots(1, 1, figsize=(15, 10), facecolor='w', edgecolor='k')
        axs.hist(h_a, bins=b)
        if self.linear_or_coh == 'linear':
            axs.set_xlabel('mm/yr')
        else:
            axs.set_xlabel('coherence')
        if self.path_save_name_hist.exists():
            os.remove(self.path_save_name_hist)
        axs.set_title(self.save_name_title)
        fig.savefig(self.path_save_name_hist, dpi=200, facecolor='w', edgecolor='k',
                    orientation='portrait', bbox_inches='tight', pad_inches=0.3)

    def run_all(self):
        self.make_plots()
        self.make_histograms()


if __name__ == '__main__':
    sys_index_var_1 = sys.argv[1]
    sys_index_var_2 = sys.argv[2]
    make_pngs_from_tif = plot_tifs_png(sys_index_var_1, sys_index_var_2)
    make_pngs_from_tif.run_all()
