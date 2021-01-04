import sys

import pandas as pd

from classes.Settings import Settings
from helpers.plot_functions import plot_fragment_colored

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def main():

    if len(sys.argv) != 2:
        print("Usage: python plot_density.py <path/to/inputfile>")
        sys.exit(1)

    inputfilename = sys.argv[1]

    settings = Settings(inputfilename)

    fragment = pd.read_csv(settings.get_avg_frag_filename())

    fig = plt.figure()
    ax: Axes3D = fig.add_subplot(111, projection='3d')

    # plot the (average of the) central group
    ax = plot_fragment_colored(ax, fragment)

    zlim = list(ax.get_zlim())

    if zlim[0] > -0.1:
        zlim[0] = -0.1
    if zlim[1] < 0.1:
        zlim[1] = 0.1

    ax.set_zlim(zlim)

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    plt.show()


if __name__ == "__main__":
    main()
