import sys, csv

from IPython.display import clear_output
from pathlib import Path
import pathlib
import matplotlib.colors as colors
import os, os.path
import rasterio
import tempfile, shutil
import time
import pandas as pd
import glob
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from rasterio.plot import show
from mpl_toolkits.axes_grid1 import make_axes_locatable
import itertools
import matplotlib as mpl

from Main.general_set_up_functions import rasterio_basic_functions as rb
from Main.general_set_up_functions import general_functions as genf


def return_overall_mean(file_list):
    mean_list = []
    for file in file_list:
        test_array = rb.tif_2_array(file).flatten()
        test_array[test_array == 0.] = np.nan
        test_array = test_array * 1000
        test_array_mean = np.nanmean(test_array)
        mean_list.append(test_array_mean)
    overall_mean = round(np.mean(mean_list), 4)
    return overall_mean, mean_list


def overall_array_hist(file_list):
    overall_array = np.array([])
    for file in file_list:
        test_array = rb.tif_2_array(file).flatten()
        test_array[test_array == 0.] = np.nan
        test_array = test_array * 1000
        overall_array = np.append(overall_array, test_array)
    return overall_array


def make_histograms(base_array, axs, color_):
    h_a = base_array.flatten()[~np.isnan(base_array.flatten())]
    h, b = np.histogram(h_a, bins=100, density=True)
    axs.hist(h_a, bins=b, alpha=0.65, color=color_, density=True)


class visualize_interferograms():
    def __init__(self, path_2_interferograms, asc_dsc, start_new):
        self.asc_dsc = str(asc_dsc)
        path_2_intf_object = Path(path_2_interferograms)
        self.full_path_2_ints = path_2_intf_object / 'Main/Processed_Data'
        self.full_path_2_msbas = path_2_intf_object.parent / 'MSBASv3_Timeseries'
        self.all_gacos_intfs = ([Path(x) for x in
                                 self.full_path_2_ints.glob("20*/*displacement_VV_rcut_mcoh_gacos.tif")])
        self.all_gacos_intfs.sort()
        self.all_non_gacos_intfs = ([Path(x) for x in
                                     self.full_path_2_ints.glob("20*/*displacement_VV_rcut_mcoh.tif")])
        self.all_non_gacos_intfs.sort()
        self.all_coherence_files = ([Path(x) for x in
                                     self.full_path_2_ints.glob("20*/*coherence_VV_rcut.tif")])
        self.all_coherence_files.sort()
        self.num_of_std = 2
        self.base_linear = 5
        self.base_coh = 0.1
        self.overall_std_gacos = None
        self.overall_std_non_gacos = None
        self.start_new = start_new
        self.temp_decis_name = self.full_path_2_msbas / ('temp_decis_' + self.asc_dsc + '.json')
        self.overall_means_temp = self.full_path_2_msbas / ('temp_means_' + self.asc_dsc + '.json')

    def make_plots(self):
        # overall_array_gacos = overall_array_hist(self.all_gacos_intfs)
        # overall_array_non_gacos = overall_array_hist(self.all_non_gacos_intfs)

        if self.overall_means_temp.exists():
            overall_means_open = genf.open_json_file(self.overall_means_temp)
            overall_mean_gacos = float(overall_means_open['mean_gacos'])
            overall_mean_non_gacos = float(overall_means_open['mean_non_gacos'])
            mean_list_gacos = [float(x) for x in overall_means_open['means_list_gacos']]
            mean_list_non_gacos = [float(x) for x in overall_means_open['means_list_non_gacos']]
        else:
            overall_mean_gacos, mean_list_gacos = return_overall_mean(self.all_gacos_intfs)
            overall_mean_non_gacos, mean_list_non_gacos = return_overall_mean(self.all_non_gacos_intfs)
        temp_means = {'mean_gacos': str(overall_mean_gacos), 'mean_non_gacos': str(overall_mean_non_gacos),
                      'means_list_gacos': [str(x) for x in mean_list_gacos],
                      'means_list_non_gacos': [str(x) for x in mean_list_non_gacos]}
        genf.write_json_from_dict(temp_means, self.overall_means_temp)

        overall_std_means_gacos = np.std(mean_list_gacos)
        overall_std_means_non_gacos = np.std(mean_list_non_gacos)
        overall_2nd_std_lower_gacos = overall_mean_gacos - 2*overall_std_means_gacos
        overall_2nd_std_upper_gacos = overall_mean_gacos + 2*overall_std_means_gacos

        overall_2nd_std_lower_non_gacos = overall_mean_non_gacos - 2*overall_std_means_non_gacos
        overall_2nd_std_upper_non_gacos = overall_mean_non_gacos + 2*overall_std_means_non_gacos

        if self.start_new == 'start_new':
            p = 0
            decision_list = [i for i in range(len(self.all_gacos_intfs))]
        else:
            decision_list_open_dict = genf.open_json_file(self.temp_decis_name)
            decision_list = [str(x) for x in decision_list_open_dict['decision']]
            p_list = [int(x) for x in decision_list if x.isdigit()]
            if not p_list:
                p = np.max([int(x.split('_')[1]) for x in decision_list])
            else:
                p = min(p_list)

        while p < len(self.all_gacos_intfs):
            gacos_intf = self.all_gacos_intfs[p]
            non_gacos_intf = self.all_non_gacos_intfs[p]
            # coherence_file = self.all_coherence_files[p]
            date = gacos_intf.parent.stem
            gacos_intf_array = rb.tif_2_array(gacos_intf)
            non_gacos_intf_array = rb.tif_2_array(non_gacos_intf)
            # coh_array = rb.tif_2_array(coherence_file)
            gacos_non_gacos_cmap = mpl.cm.rainbow
            # coh_cmap = mpl.cm.Blues.reversed()
            # fig, ((axs1, axs2), (axs3, axs4)) = plt.subplots(2, 2, figsize=(27, 18), facecolor='w', edgecolor='k')
            fig, (axs1, axs2) = plt.subplots(1, 2, figsize=(27, 18), facecolor='w', edgecolor='k')
            gacos_intf_array[gacos_intf_array == 0.] = np.nan
            gacos_intf_array = gacos_intf_array * 1000
            mean_tif_gacos_intf = np.nanmean(gacos_intf_array)
            std_tif_gacos_intf = np.nanstd(gacos_intf_array)
            lower_bound_u_z_gacos_intf, upper_bound_u_z_gacos_intf = genf.set_bounds(mean_tif_gacos_intf,
                                                                                     std_tif_gacos_intf,
                                                                                     self.num_of_std,
                                                                                     self.base_linear)
            norm_gacos_intf = mpl.colors.Normalize(vmin=lower_bound_u_z_gacos_intf,
                                                    vmax=upper_bound_u_z_gacos_intf)

            non_gacos_intf_array[non_gacos_intf_array == 0.] = np.nan
            non_gacos_intf_array = non_gacos_intf_array * 1000
            mean_tif_non_gacos_intf = np.nanmean(non_gacos_intf_array)
            std_tif_non_gacos_intf = np.nanstd(non_gacos_intf_array)

            lower_bound_u_z_non_gacos_intf, upper_bound_u_z_non_gacos_intf = genf.set_bounds(mean_tif_non_gacos_intf,
                                                                                             std_tif_non_gacos_intf,
                                                                                             self.num_of_std,
                                                                                             self.base_linear)
            # if lower_bound_u_z_gacos_intf < lower_bound_u_z_non_gacos_intf:
            #     lower_bound_all = lower_bound_u_z_gacos_intf
            # else:
            #     lower_bound_all = lower_bound_u_z_non_gacos_intf
            # if upper_bound_u_z_gacos_intf > upper_bound_u_z_non_gacos_intf:
            #     upper_bound_all = upper_bound_u_z_gacos_intf
            # else:
            #     upper_bound_all = upper_bound_u_z_non_gacos_intf

            lower_bound_overall_gacos, upper_bound_overall_gacos = genf.set_bounds(overall_mean_gacos,
                                                                                   overall_std_means_gacos,
                                                                                   2.3,
                                                                                   20)
            lower_bound_overall_non_gacos, upper_bound_overall_non_gacos = genf.set_bounds(overall_mean_non_gacos,
                                                                                           overall_std_means_non_gacos,
                                                                                           2.3,
                                                                                           20)
            if lower_bound_overall_gacos < lower_bound_overall_non_gacos:
                overall_lower_bound = lower_bound_overall_gacos
            else:
                overall_lower_bound = lower_bound_overall_non_gacos
            if upper_bound_overall_gacos > upper_bound_overall_non_gacos:
                overall_upper_bound = upper_bound_overall_gacos
            else:
                overall_upper_bound = upper_bound_overall_non_gacos

            norm_non_gacos_intf = (mpl.colors.Normalize(vmin=lower_bound_u_z_non_gacos_intf,
                                                        vmax=upper_bound_u_z_non_gacos_intf))

            # norm_all = (mpl.colors.Normalize(vmin=lower_bound_all,
            #                                   vmax=upper_bound_all))

            # coh_array[coh_array == 0.] = np.nan
            # coh_array = coh_array * 1000
            # mean_tif_coh = np.nanmean(coh_array)
            # std_tif_coh = np.nanstd(coh_array)
            # lower_bound_u_z_coh, upper_bound_u_z_coh = genf.set_bounds(mean_tif_coh,
            #                                                            std_tif_coh,
            #                                                            self.num_of_std,
            #                                                            self.base_coh)
            # norm_coh = (mpl.colors.Normalize(vmin=lower_bound_u_z_coh,
            #                                  vmax=upper_bound_u_z_coh))

            im1 = (axs1.imshow(gacos_intf_array,
                               vmin=lower_bound_u_z_gacos_intf, vmax=upper_bound_u_z_gacos_intf,
                               cmap=gacos_non_gacos_cmap,
                               norm=norm_gacos_intf, origin="upper"))
            im2 = (axs2.imshow(non_gacos_intf_array,
                               vmin=lower_bound_u_z_non_gacos_intf, vmax=upper_bound_u_z_non_gacos_intf,
                               cmap=gacos_non_gacos_cmap,
                               norm=norm_non_gacos_intf, origin="upper"))
            # im3 = (axs3.imshow(coh_array,
            #                    vmin=lower_bound_u_z_coh, vmax=upper_bound_u_z_coh,
            #                    cmap=coh_cmap,
            #                    norm=norm_coh, origin="upper"))

                # make_histograms(overall_array_gacos, axs3, 'b')
                # make_histograms(overall_array_non_gacos, axs3, 'r')

                # make_histograms(gacos_intf_array, axs4, 'g')
                # make_histograms(non_gacos_intf_array, axs4, 'c', )
                #
                # axs3.set_xlim(left=overall_lower_bound, right=overall_upper_bound)
                #
                # axs4.set_xlim(left=lower_bound_all, right=upper_bound_all)
                # axs4.legend(['g_s', 'ng_s'])
                #
                # axs3.axvline(x=overall_mean_gacos, color='b')
                # axs3.axvline(x=overall_mean_non_gacos, color='r')
                # axs3.legend(['g_o', 'ng_o'])
                # axs3.axvline(x=overall_2nd_std_lower_gacos, color='b', linestyle='--')
                # axs3.axvline(x=overall_2nd_std_upper_gacos, color='b', linestyle='--')
                # axs3.axvline(x=overall_2nd_std_lower_non_gacos, color='r', linestyle='--')
                # axs3.axvline(x=overall_2nd_std_upper_non_gacos, color='r', linestyle='--')
                # axs3.axvline(x=mean_tif_gacos_intf, color='g')
                # axs3.axvline(x=mean_tif_non_gacos_intf, color='c')
                # axs4.axvline(x=overall_mean_gacos, color='b')
                # axs4.axvline(x=overall_mean_non_gacos, color='r')
                # axs4.axvline(x=mean_tif_gacos_intf, color='g')
                # axs4.axvline(x=mean_tif_non_gacos_intf, color='c')


            divider1 = make_axes_locatable(axs1)
            divider2 = make_axes_locatable(axs2)
            # divider3 = make_axes_locatable(axs3)
            cax1 = divider1.append_axes("right", size="5%", pad=0.05)
            cax2 = divider2.append_axes("right", size="5%", pad=0.05)
            # cax3 = divider3.append_axes("right", size="5%", pad=0.05)
            fig.colorbar(im1, ax=axs1, cax=cax1)
            fig.colorbar(im2, ax=axs2, cax=cax2)
            # fig.colorbar(im3, ax=axs3, cax=cax3)
            axs1.set_title('gacos_intf_' + date)
            axs2.set_title('non_gacos_intf_' + date)
            # axs3.set_title('coh_' + date)
            axs1.set_xticklabels([])
            axs1.set_yticklabels([])
            axs2.set_xticklabels([])
            axs2.set_yticklabels([])
            axs1.text(0.01, 0.97, ('gacos_std=' + str(std_tif_gacos_intf)), fontsize=12,
                      transform=axs1.transAxes)
            axs1.text(0.23, 0.97, ('gacos_overall_mean=' + str(overall_mean_gacos)), fontsize=12,
                      transform=axs1.transAxes)
            axs1.text(0.23, 0.91, ('gacos_mean=' + str(mean_tif_gacos_intf)), fontsize=12,
                      transform=axs1.transAxes)
            axs1.text(0.01, 0.91, ('gacos_min=' + str(lower_bound_u_z_gacos_intf)), fontsize=12,
                      transform=axs1.transAxes)
            axs1.text(0.01, 0.85, ('gacos_max=' + str(upper_bound_u_z_gacos_intf)), fontsize=12,
                      transform=axs1.transAxes)

            axs2.text(0.01, 0.97, ('non_gacos_std=' + str(std_tif_non_gacos_intf)), fontsize=12,
                      transform=axs2.transAxes)
            axs2.text(0.23, 0.97, ('non_gacos_overall_mean=' + str(overall_mean_non_gacos)), fontsize=12,
                      transform=axs2.transAxes)
            axs2.text(0.23, 0.91, ('non_gacos_mean=' + str(mean_tif_non_gacos_intf)), fontsize=12,
                      transform=axs2.transAxes)
            axs2.text(0.01, 0.91, ('non_gacos_min=' + str(lower_bound_u_z_non_gacos_intf)), fontsize=12,
                      transform=axs2.transAxes)
            axs2.text(0.01, 0.85, ('non_gacos_max=' + str(upper_bound_u_z_non_gacos_intf)), fontsize=12,
                      transform=axs2.transAxes)
            # axs3.set_xticklabels([])
            # axs3.set_yticklabels([])
            fig.tight_layout()
            plt.show()
            print('Move Forward or backward: F = forward : B = backward')
            print("Current location in loop:", p)
            choice_1 = input()
            print("current decision in storage:", decision_list[p])
            print('Keep or throw away: 1 = keep gacos: 2 = keep non gacos: 3 = Trash')
            choice_2 = input()
            decision = str(choice_2) + "_" + str(p)
            decision_list[p] = decision
            if choice_1 == 'F':
                p = p + 1
            else:
                p = p - 1
            if p < 0:
                print('p is below 0 you cannot go backward any further')
                p = 0
            clear_output(wait=True)

            temp_dict = {'decision': decision_list}
            genf.write_json_from_dict(temp_dict, self.temp_decis_name)

        gacos_trash_list = []
        non_gacos_trash_list = []
        gacos_or_non_gacos_list = []
        for decis in decision_list:
            keep_or_trash = decis.split('_')[0]
            number = int(decis.split('_')[1])
            if keep_or_trash == '1' or keep_or_trash == '2':
                if keep_or_trash == '1':
                    save_var = 'g'
                else:
                    save_var = 'ng'
                gacos_or_non_gacos_list.append(save_var)
            else:
                gacos_trash = self.all_gacos_intfs[number]
                non_gacos_trash = self.all_non_gacos_intfs[number]
                gacos_trash_list.append(str(gacos_trash))
                non_gacos_trash_list.append(str(non_gacos_trash))
        gacos_trash_name = self.full_path_2_msbas / ('visual_remove_gacos_' + self.asc_dsc + '.json')
        # non_gacos_trash_name = self.full_path_2_msbas / ('visual_remove_non_gacos_' + self.asc_dsc + '.json')
        gacos_or_non_gacos_name = self.full_path_2_msbas / ('gacos_or_non_gacos_use_' + self.asc_dsc + '.json')
        gacos_trash_dict = {'trash_files': gacos_trash_list}
        # non_gacos_trash_dict = {'trash_files': non_gacos_trash_list}
        gacos_or_non_gacos_dict = {'g_ng': gacos_or_non_gacos_list}
        genf.write_json_from_dict(gacos_trash_dict, gacos_trash_name)
        # genf.write_json_from_dict(non_gacos_trash_dict, non_gacos_trash_name)
        genf.write_json_from_dict(gacos_or_non_gacos_dict, gacos_or_non_gacos_name)

    def run_all(self):
        self.make_plots()

# if __name__ == '__main__':
#     sys_index_var_1 = sys.argv[1]
#     sys_index_var_2 = sys.argv[2]
#     visualize_tifs = visualize_interferograms(sys_index_var_1, sys_index_var_2)
#     visualize_tifs.run_all()
