# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This script is part of the quantification pipeline of 3D experimental data of crystal structures that I wrote for my
# thesis in the Master Computational Science, University of Amsterdam, 2021.
#
# `plot_avg_fragment`
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import sys

import pandas as pd

from classes.Settings import Settings
from helpers.plot_functions import plot_fragment_colored

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from constants.paths import WORKDIR


def main():

    if len(sys.argv) != 2:
        print("Usage: python plot_density.py <path/to/inputfile>")
        sys.exit(1)

    inputfilename = sys.argv[1]

    settings = Settings(WORKDIR, inputfilename)
    plot_avg_fragment(settings)


def plot_avg_fragment(settings):
    fragment = pd.read_csv(settings.get_avg_frag_filename())

    fig = plt.figure()
    ax: Axes3D = fig.add_subplot(111, projection='3d')

    ax.set_title(f"Average {settings.central_name}--{settings.contact_name} Fragment")

    # plot the (average of the) central group
    ax = plot_fragment_colored(ax, fragment)

    xlim, ylim, zlim = list(ax.get_xlim()), list(ax.get_ylim()), list(ax.get_zlim())
    minn = min([xlim[0], ylim[0], zlim[0]]) - 1
    maxx = max([xlim[1], ylim[1], zlim[1]]) + 1

    ax.set_xlim((minn, maxx))
    ax.set_ylim((minn, maxx))
    ax.set_zlim((minn, maxx))

    ax.set_xlabel('X axis ($\\AA$)')
    ax.set_ylabel('Y axis ($\\AA$)')
    ax.set_zlabel('Z axis ($\\AA$)')

    plt.show()


if __name__ == "__main__":
    main()
