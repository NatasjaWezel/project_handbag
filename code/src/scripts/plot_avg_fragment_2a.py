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

    fragment = pd.read_csv(settings.get_avg_frag_filename())

    fig = plt.figure()
    ax: Axes3D = fig.add_subplot(111, projection='3d')

    ax.set_title(f"Average {settings.central_group_name}--{settings.contact_group_name} Fragment")

    # plot the (average of the) central group
    ax = plot_fragment_colored(ax, fragment)

    xlim, ylim, zlim = list(ax.get_xlim()), list(ax.get_ylim()), list(ax.get_zlim())
    minn = min([xlim[0], ylim[0], zlim[0]])
    maxx = max([xlim[1], ylim[1], zlim[1]])

    ax.set_xlim((minn, maxx))
    ax.set_ylim((minn, maxx))
    ax.set_zlim((minn, maxx))

    ax.set_xlabel('X axis ($\\AA$)')
    ax.set_ylabel('Y axis ($\\AA$)')
    ax.set_zlabel('Z axis ($\\AA$)')

    plt.show()


if __name__ == "__main__":
    main()
