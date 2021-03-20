# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This is a script that I wrote for my master thesis
# It loads the coordinates of the aligned fragments. It then divides the
# surrounding space into a number of bins, depending on which resolution is
# set. It counts how many of the contact atoms/ centers of contact groups are
# are in each bin and normalizes that by the total amount of contact atoms or
# groups. Then a plot is made that shows the density of the contacts in "4D".
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import sys

import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D

from classes.Settings import Settings
from helpers.plot_functions import plot_density, plot_fragment_colored, plot_vdw_spheres

from constants.paths import WORKDIR


def main():

    if len(sys.argv) != 5:
        print("Usage: python analyze_density.py <path/to/inputfile> <resolution> <threshold> <atom or center to count>")
        sys.exit(1)

    settings = Settings(WORKDIR, sys.argv[1])

    # resolution of the bins, in Angstrom
    settings.set_resolution(float(sys.argv[2]))
    settings.set_threshold(float(sys.argv[3]))

    settings.set_contact_reference_point(sys.argv[4])

    try:
        avg_fragment = pd.read_csv(settings.get_avg_frag_filename())
        density_df = pd.read_hdf(settings.get_density_df_filename(), settings.get_density_df_key())
    except (FileNotFoundError, KeyError) as exception:
        print(exception)
        print("Run avg_frag and calc_density first")
        sys.exit(1)

    make_density_plot(avg_fragment, density_df, settings)


def make_density_plot(avg_fragment, density_df, settings):
    plot_spheres = True
    plotname = settings.get_density_plotname()
    fig = plt.figure()
    ax: Axes3D = fig.add_subplot(111, projection='3d')

    ax = plot_fragment_colored(ax, avg_fragment)
    p, ax = plot_density(ax=ax, df=density_df, settings=settings)

    if plot_spheres:
        ax, _ = plot_vdw_spheres(avg_fragment, ax)

    title = f"{settings.central_name}--{settings.contact_name} ({settings.contact_rp}) density\n"
    title += f"Resolution: {settings.resolution :.2f}, fraction: {settings.threshold :.2f}"

    ax.set_title(title)

    xlim, ylim, zlim = list(ax.get_xlim()), list(ax.get_ylim()), list(ax.get_zlim())
    minn = min([xlim[0], ylim[0], zlim[0]])
    maxx = max([xlim[1], ylim[1], zlim[1]])

    ax.set_xlim((minn, maxx))
    ax.set_ylim((minn, maxx))
    ax.set_zlim((minn, maxx))

    ax.view_init(elev=90, azim=90)

    cbar = plt.colorbar(p)

    # format colorbar ticks into percentages
    cbar.set_ticks(cbar.ax.get_yticks())
    cbar.ax.set_yticklabels(['{:.2f}%'.format(x * 100) for x in cbar.get_ticks()])

    plt.savefig(plotname)
    plt.show()


if __name__ == "__main__":
    main()
