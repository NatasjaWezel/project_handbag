# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This script is part of the quantification pipeline of 3D experimental data of crystal structures that I wrote for my
# thesis in the Master Computational Science, University of Amsterdam, 2021.
#
# `density_slider`
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import math
import sys

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.widgets import Button, Slider
from mpl_toolkits.mplot3d import Axes3D
from constants.paths import WORKDIR
from classes.Settings import Settings
from helpers.plot_functions import plot_fragment_colored, plot_density

from constants.constants import STANDARD_THRESHOLD, STANDARD_RES


def main():

    if len(sys.argv) != 3:
        print("Usage: python analyze_density.py <path/to/inputfile> <atom or center to count>")
        sys.exit(1)

    settings = Settings(WORKDIR, sys.argv[1])

    # resolution of the bins, in Angstrom
    settings.set_resolution(STANDARD_RES)
    settings.set_threshold(STANDARD_THRESHOLD)
    settings.set_contact_reference_point(sys.argv[2])

    avg_fragment = pd.read_csv(settings.get_avg_frag_filename())

    make_density_slider_plot(avg_fragment, settings)


def make_density_slider_plot(avg_fragment, settings):
    df = pd.read_hdf(settings.get_density_df_filename(), settings.get_density_df_key())
    df[settings.contact_rp] = df[settings.contact_rp] / df[settings.contact_rp].sum()

    maximum = df[settings.contact_rp].max()

    fig = plt.figure(figsize=(8, 5))

    plt.subplots_adjust(bottom=0.30)
    ax: Axes3D = fig.add_subplot(111, projection='3d')

    ax = plot_fragment_colored(ax, avg_fragment)

    global p
    p, ax = plot_density(ax=ax, df=df, settings=settings)

    title = f"{settings.central_name}--{settings.contact_name} ({settings.contact_rp}) density\n"
    title += f"Resolution: {settings.resolution :.2f}, fraction: {settings.threshold :.2f}"

    ax.set_title(title)

    xlim, ylim, zlim = ax.get_xlim(), ax.get_ylim(), ax.get_zlim()
    xlim, ylim, zlim = (math.floor(xlim[0]), math.ceil(xlim[1])), (math.floor(ylim[0]), math.ceil(ylim[1])),\
                       (math.floor(zlim[0]), math.ceil(zlim[1]))

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_zlim(zlim)

    cax = fig.add_axes([0.8, 0.25, 0.03, 0.6])
    cbar = fig.colorbar(p, cax=cax, pad=0.2)

    # format colorbar ticks into percentages
    cbar.set_ticks(cbar.ax.get_yticks())
    cbar.ax.set_yticklabels(['{:.2f}%'.format(x * 100) for x in cbar.get_ticks()])

    axcolor = 'lightgoldenrodyellow'
    ax_resolution = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
    resolution = Slider(ax_resolution, 'Res', 0.2, 1, valinit=settings.resolution, valstep=0.05)

    ax_lowerlim = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    lowerlim = Slider(ax_lowerlim, 'Lim', 0, 1, valinit=settings.threshold, valstep=0.01)

    percentage = round(df[df[settings.contact_rp] >=
                          settings.threshold * maximum][settings.contact_rp].sum() * 100, 2)
    text_holder = fig.text(.25, .05, f"Showing {percentage}% of all data")

    # update everything
    def update_res(val):
        print("\nChanged resolution to:", round(val, 2))
        settings.set_resolution(round(val, 2))
        print(f"Threshold: {settings.threshold}")
        df = pd.read_hdf(settings.get_density_df_filename(), settings.get_density_df_key())
        df[settings.contact_rp] = df[settings.contact_rp] / df[settings.contact_rp].sum()
        maximum = df[settings.contact_rp].max()

        global p

        ax.cla()
        cax.cla()

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_zlim(zlim)

        title = f"{settings.central_name}--{settings.contact_name} ({settings.contact_rp}) density\n"
        title += f"Resolution: {settings.resolution :.2f}, fraction: {settings.threshold :.2f}"
        ax.set_title(title)

        p, ax1 = plot_density(ax=ax, df=df, settings=settings)
        ax1 = plot_fragment_colored(ax1, avg_fragment)

        percentage = round(df[df[settings.contact_rp] >=
                              settings.threshold * maximum][settings.contact_rp].sum() * 100, 2)
        text_holder.set_text(f"Showing {percentage}% of all data")

        cbar = fig.colorbar(p, cax=cax, pad=0.2)

        # format colorbar ticks into percentages
        cbar.set_ticks(cbar.ax.get_yticks())
        cbar.ax.set_yticklabels(['{:.2f}%'.format(x * 100) for x in cbar.get_ticks()])

    def update_lim(val):
        print("\nChanged threshold too:", round(val, 2))
        print("Resolution:", settings.resolution)
        global p

        df = pd.read_hdf(settings.get_density_df_filename(), settings.get_density_df_key())
        df[settings.contact_rp] = df[settings.contact_rp] / df[settings.contact_rp].sum()
        maximum = df[settings.contact_rp].max()

        settings.set_threshold(round(val, 2))

        ax.cla()
        cax.cla()

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_zlim(zlim)

        _, ax1 = plot_density(ax=ax, df=df, settings=settings)
        ax1 = plot_fragment_colored(ax1, avg_fragment)

        title = f"{settings.central_name}--{settings.contact_name} ({settings.contact_rp}) density\n"
        title += f"Resolution: {settings.resolution :.2f}, fraction: {settings.threshold :.2f}"
        ax.set_title(title)

        no_bins_cluster = len(df[df[settings.contact_rp] >= settings.threshold * maximum])
        bins_cluster = df[df[settings.contact_rp] >= settings.threshold * maximum]

        percentage = round(bins_cluster[settings.contact_rp].sum() * 100, 2)

        text_holder.set_text(f"Showing {percentage}% of all data")

        print(f"Volume: {no_bins_cluster*settings.resolution**3 :.2f}")

        cbar = fig.colorbar(p, cax=cax, pad=0.2)

        # format colorbar ticks into percentages
        cbar.set_ticks(cbar.ax.get_yticks())
        cbar.ax.set_yticklabels(['{:.2f}%'.format(x * 100) for x in cbar.get_ticks()])

    resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

    def reset(event):
        resolution.reset()
        lowerlim.reset()

    button.on_clicked(reset)
    resolution.on_changed(update_res)
    lowerlim.on_changed(update_lim)

    plt.show()


if __name__ == "__main__":
    main()
