import sys

from classes.Settings import Settings
from helpers.geometry_helpers import make_avg_fragment_if_not_exists
from helpers.helpers import read_results_alignment
from helpers.plot_functions import plot_fragment_colored

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def main():

    if len(sys.argv) != 2:
        print("Usage: python plot_density.py <path/to/inputfile>")
        sys.exit(1)

    inputfilename = sys.argv[1]

    settings = Settings(inputfilename)

    aligned_fragments_df = read_results_alignment(settings.get_aligned_csv_filename())

    fragment = make_avg_fragment_if_not_exists(settings, aligned_fragments_df)

    fig = plt.figure()
    ax: Axes3D = fig.add_subplot(111, projection='3d')

    # plot the (average of the) central group
    ax = plot_fragment_colored(ax, fragment)

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    plt.show()


if __name__ == "__main__":
    main()
