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
import os

import numpy as np
import pandas as pd
from tqdm import tqdm

import time

from classes.Settings import Settings
from helpers.geometry_helpers import (get_vdw_distance_contact,
                                      make_avg_fragment_if_not_exists)
from helpers.helpers import read_results_alignment


def main():

    if len(sys.argv) != 4:
        print("Usage: python analyze_density.py <path/to/inputfile> <resolution> <atom or center to count>")
        sys.exit(1)

    settings = Settings(sys.argv[1])

    # resolution of the bins, in Angstrom
    settings.set_resolution(float(sys.argv[2]))

    settings.set_atom_to_count(sys.argv[3])

    try:
        density_df = pd.read_hdf(settings.get_density_df_filename(), settings.get_density_df_key())
    except (FileNotFoundError, KeyError):
        print("Run calc_density first")
        sys.exit(1)

    df = read_results_alignment(settings.get_aligned_csv_filename())
    avg_fragment = make_avg_fragment_if_not_exists(settings, df)

    # grab only the atoms that are in the contact groups
    df = df[df.in_central_group == False]

    empty_bins, bins_80 = calculate_80_percent(density_df, settings)

    vdw_distance_contact = get_vdw_distance_contact(df, settings)
    bins_in_vdw = count_bins_in_vdw(density_df, avg_fragment, settings, vdw_distance_contact)

    if not os.path.exists(settings.get_directionality_results_filename()):
        with open(settings.get_directionality_results_filename(), 'w', newline='') as resultsfile:
            writer = csv.writer(resultsfile)
            writer.writerow(["centralgroup", "contactgroup", "to_count_contact", "resolution", "datapoints", "bins",
                             "empty_bins", "bins_in_vdw", "bins_80"])

    with open(settings.get_directionality_results_filename(), 'a', newline='') as resultsfile:
        writer = csv.writer(resultsfile)

        writer.writerow([settings.central_group_name, settings.contact_group_name, settings.to_count_contact,
                         settings.resolution, density_df[settings.to_count_contact].sum(), len(density_df),
                         empty_bins, bins_in_vdw, bins_80])


def count_bins_in_vdw(density_df, avg_fragment, settings, vdw_distance_contact, extra_radius=0.5):
    # put the centroids instead of the boundaries in the bin_coordinates
    bin_coordinates = np.array([density_df.xstart + 0.5 * settings.resolution,
                                density_df.ystart + 0.5 * settings.resolution,
                                density_df.zstart + 0.5 * settings.resolution])

    bin_coordinates = bin_coordinates.T

    in_vdw_volume = np.zeros(len(bin_coordinates))

    points_avg_f = np.array([avg_fragment.atom_x, avg_fragment.atom_y, avg_fragment.atom_z, avg_fragment.vdw_radius]).T

    print("Calculating if bins are in vdw volume central group for resolution ", settings.resolution)
    time.sleep(0.5)
    for i in tqdm(range(len(points_avg_f))):
        indices = np.transpose(np.where(in_vdw_volume == 0))

        avg_f_p = points_avg_f[i]

        in_vdw_volume = calc_distances(in_vdw_volume, bin_coordinates, avg_f_p, indices, vdw_distance_contact, extra_radius)

    total = np.sum(in_vdw_volume)

    return total


def calc_distances(in_vdw_volume, bin_coordinates, avg_f_p, indices, vdw_distance_contact, extra_radius):

    for idx in indices:
        bin_point = bin_coordinates[idx[0]]
        distance = np.sum((bin_point - avg_f_p[:3])**2, axis=0)**0.5

        # TODO: this hardcoded 0.5
        if distance < avg_f_p[3] + vdw_distance_contact + extra_radius:
            in_vdw_volume[idx] = 1

    return in_vdw_volume


def calculate_80_percent(df, settings):
    total_atoms = df[settings.to_count_contact].sum()

    # take non-empty bins and sort them from highest fraction to lowest fraction
    population = list(df[df[settings.to_count_contact] > 0][settings.to_count_contact])
    population.sort(reverse=True)

    empty_bins = len(df) - len(population)

    fraction = 0
    i = 0

    print("Calculating fraction of bins containing 80% of data")
    while fraction < 0.8:
        fraction += population[i]/total_atoms
        i += 1

    return empty_bins, i


if __name__ == "__main__":
    main()
