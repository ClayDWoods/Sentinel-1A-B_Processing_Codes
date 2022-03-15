import sys, os

sys.path.append(".")

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


def make_max_curv_plot(lam_list, curv, plot_directory, l_value):
    fig1 = plt.figure()
    ax1 = plt.axes()
    ax1.scatter(lam_list, curv, color='b')
    plot_name = plot_directory / ('Max_Curve_Plot_' + str(l_value) + '.png')
    x_name = 'lambda'
    y_name = 'curvature'
    title = 'Curvature of L_Curve'
    ax1.set_xlabel(x_name)
    ax1.set_ylabel(y_name)
    ax1.set_title(title)
    if plot_name.exists():
        os.remove(plot_name)
    fig1.savefig(plot_name)


def make_L_curv_plot(log_x_1, log_x_2, log_x_3, log_axy_1, log_axy_2, log_axy_3, plot_directory):
    plot_name = plot_directory / 'L_Curve_Plot'
    fig1 = plt.figure()
    ax1 = plt.axes()
    x_name = 'log(||Ax-Y||)'
    y_name = 'log(||x||)'
    title = 'L_Curve'
    if log_axy_1 and log_x_1:
        log_axy_1_plot, log_x_1_plot = np.log10(np.square(log_axy_1))/2, np.log10(np.square(log_x_1))/2
        ax1.scatter(log_axy_1_plot, log_x_1_plot, color='b', label='L=1')
    if log_axy_2 and log_x_2:
        log_axy_2_plot, log_x_2_plot = np.log10(np.square(log_axy_2))/2, np.log10(np.square(log_x_2))/2
        ax1.scatter(log_axy_2_plot, log_x_2_plot, color='r', label='L=2')
    if log_axy_3 and log_x_3:
        log_axy_3_plot, log_x_3_plot = np.log10(np.square(log_axy_3))/2, np.log10(np.square(log_x_3))/2
        ax1.scatter(log_axy_3_plot, log_x_3_plot, color='g', label='L=3')
    ax1.set_xlabel(x_name)
    ax1.set_ylabel(y_name)
    ax1.set_title(title)
    ax1.legend()
    if plot_name.exists():
        os.remove(plot_name)
    fig1.savefig(plot_name)
