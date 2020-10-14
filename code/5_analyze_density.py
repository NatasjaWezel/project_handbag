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

import csv
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm

from classes.Settings import Settings
from helpers.geometry_helpers import make_avg_fragment_if_not_exists, get_vdw_distance_contact
from helpers.helpers import read_results_alignment
from helpers.plot_functions import plot_density, plot_fragment_colored


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

    df = df[df.in_central_group is False]

    vdw_distance_contact = get_vdw_distance_contact(df, settings)

    try:
        density_df = pd.read_hdf(settings.get_density_df_filename(), settings.get_density_df_key())
    except (FileNotFoundError, KeyError):
        print("Run calc_density first")
        sys.exit(1)

    print(density_df.describe())

    count_bins_in_vdw(density_df, avg_fragment, settings, vdw_distance_contact)
    calculate_80_percent(density_df, settings)
    make_plot(avg_fragment, density_df, settings)


def count_bins_in_vdw(density_df, avg_fragment, settings, vdw_distance_contact):
    bin_coordinates = np.array([density_df.xstart, density_df.ystart, density_df.zstart])

    index_x, index_y, index_z = np.where(bin_coordinates[0] < 0),\
        np.where(bin_coordinates[1] < 0),\
        np.where(bin_coordinates[2] < 0)

    bin_coordinates[0][index_x] += settings.resolution
    bin_coordinates[1][index_y] += settings.resolution
    bin_coordinates[2][index_z] += settings.resolution

    bin_coordinates = bin_coordinates.T

    in_vdw_volume = np.zeros(len(bin_coordinates))

    points_avg_f = np.array([avg_fragment.atom_x, avg_fragment.atom_y, avg_fragment.atom_z, avg_fragment.vdw_radius]).T

    for i in tqdm(range(len(points_avg_f))):
        indices = np.transpose(np.where(in_vdw_volume == 0))

        avg_f_p = points_avg_f[i]

        calc_distances(in_vdw_volume, bin_coordinates, avg_f_p, indices, vdw_distance_contact)

    total = np.sum(in_vdw_volume)

    print(total, '/', len(in_vdw_volume), 'bins in overlapping vdw volume + 0.5')
    print('Thats a fraction of:', total/len(in_vdw_volume) * 100)

    return total


def calc_distances(in_vdw_volume, bin_coordinates, avg_f_p, indices, vdw_distance_contact):

    for idx in indices:
        bin_point = bin_coordinates[idx[0]]
        distance = np.sum((bin_point - avg_f_p[:3])**2, axis=0)**0.5

        # TODO: this hardcoded 0.5
        if distance < avg_f_p[3] + vdw_distance_contact + 0.5:
            in_vdw_volume[idx] = 1

    return in_vdw_volume


def make_plot(avg_fragment, density_df, settings):
    plotname = settings.get_density_plotname()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax = plot_fragment_colored(ax, avg_fragment)

    p, ax = plot_density(ax=ax, df=density_df, settings=settings)

    ax.set_title("4D density plot\n Resolution: " + str(settings.resolution))

    fig.colorbar(p)
    plt.savefig(plotname)
    plt.show()


def calculate_80_percent(df, settings):
    total_atoms = df[settings.to_count_contact].sum()

    # take non-empty bins and sort them from highest fraction to lowest fraction
    population = list(df[df[settings.to_count_contact] > 0][settings.to_count_contact])
    population.sort(reverse=True)

    fraction = 0
    i = 0

    while fraction < 0.8:
        fraction += population[i]/total_atoms
        i += 1

    non_empty_bins = len(population)
    bins = len(df)

    with open(settings.get_directionality_results_filename(), 'a', newline='') as resultsfile:
        writer = csv.writer(resultsfile)

        writer.writerow([settings.central_group_name, settings.contact_group_name, settings.to_count_contact,
                        fraction, i, non_empty_bins])

    print(str(round(fraction * 100, 2)) + "% of contact group atoms is in " + str(round(i/non_empty_bins*100, 2)) +
          "% of the non-empty bins")

    print(str(round(fraction * 100, 2)) + "% of contact group atoms is in " + str(round(i/bins*100, 2)) +
          "% of the total bins")


if __name__ == "__main__":
    main()
