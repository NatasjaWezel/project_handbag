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
from helpers.geometry_helpers import make_avg_fragment_if_not_exists
from helpers.helpers import read_results_alignment
from helpers.plot_functions import plot_density, plot_fragment_colored

from sklearn.cluster import KMeans
import numpy as np


def main():

    if len(sys.argv) != 4:
        print("Usage: python analyze_density.py <path/to/inputfile> <resolution> <atom or center to count>")
        sys.exit(1)

    settings = Settings(sys.argv[1])

    # resolution of the bins, in Angstrom
    settings.set_resolution(float(sys.argv[2]))

    settings.set_atom_to_count(sys.argv[3])

    df = read_results_alignment(settings.get_aligned_csv_filename())
    avg_fragment = make_avg_fragment_if_not_exists(settings, df)

    df = df[df.in_central_group == False]

    try:
        density_df = pd.read_hdf(settings.get_density_df_filename(), settings.get_density_df_key())
    except (FileNotFoundError, KeyError):
        print("Run calc_density first")
        sys.exit(1)

    print(df.columns)
    print(density_df.columns)

    threshold = density_df[settings.to_count_contact].max() * 5 / 100

    print(len(density_df[density_df[settings.to_count_contact] > float(threshold)]))

    indices = np.where(np.array(density_df[settings.to_count_contact]) >= threshold)

    amount_of_clusters = 2
    calc_clusters(density_df, indices, amount_of_clusters)

    # normalize
    density_df["centroid_normalized"] = density_df["centroid"] / density_df.centroid.max()
    
    densities = sorted(list(density_df["centroid_normalized"]))

    plt.plot(range(0, len(densities)), densities)
    plt.yscale("log")
    plt.show()

    make_cluster_plot(avg_fragment, density_df[density_df[settings.to_count_contact] > float(threshold)], settings)
    make_plot(avg_fragment, density_df, settings)


def calc_clusters(density_df, indices, amount_of_clusters):
    density_df["cluster"] = - 1
    df = density_df.iloc[indices]

    X = np.transpose(np.array([df.xstart, df.ystart, df.zstart]))

    kmeans = KMeans(n_clusters=amount_of_clusters, random_state=1)
    kmeans.fit(X)

    for index, label in zip(indices[0], kmeans.labels_):
        density_df.loc[index, "cluster"] = label

    return density_df


def make_plot(avg_fragment, density_df, settings):
    plotname = settings.get_density_plotname()
    fig = plt.figure()
    ax: Axes3D = fig.add_subplot(111, projection='3d')

    ax = plot_fragment_colored(ax, avg_fragment)

    p, ax = plot_density(ax=ax, df=density_df, settings=settings)

    ax.set_title("4D density plot\n Resolution: " + str(settings.resolution))

    fig.colorbar(p)
    plt.savefig(plotname)
    plt.show()


def make_cluster_plot(avg_fragment, density_df, settings):
    colors = ["grey", "red", "green"]

    density_df[settings.to_count_contact] = density_df[settings.to_count_contact] / density_df[settings.to_count_contact].sum()

    density_df = density_df[density_df[settings.to_count_contact] > 0.00001]

    density_df["cluster_color"] = [colors[i + 1] for i in list(density_df.cluster)]

    fig = plt.figure()
    ax: Axes3D = fig.add_subplot(111, projection='3d')

    ax = plot_fragment_colored(ax, avg_fragment)

    ax.scatter(list(density_df.xstart), list(density_df.ystart), list(density_df.zstart),
               color=list(density_df.cluster_color))

    ax.set_title("4D density plot\n Resolution: " + str(settings.resolution))

    plt.show()


if __name__ == "__main__":
    main()
