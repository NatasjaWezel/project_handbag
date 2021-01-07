# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This is a script that I wrote for my master thesis
# It loads the coordinates of the aligned fragments. It then divides the
# surrounding space into a number of bins, depending on which resolution is
# set. It counts how many of the contact atoms/ centers of contact groups are
# are in each bin and normalizes that by the total amount of contact atoms or
# groups.
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import sys
import time

from numba import jit

import numpy as np
import pandas as pd

from classes.Settings import Settings
from helpers.density_helpers import find_available_volume, prepare_df
from helpers.geometry_helpers import (get_vdw_distance_contact,
                                      make_coordinate_df)


def main():

    if len(sys.argv) != 4:
        print("Usage: python plot_density.py <path/to/inputfile> <resolution> <atom or center to count>")
        sys.exit(1)

    t0 = time.time()

    settings = Settings(sys.argv[1])
    settings.set_atom_to_count(sys.argv[3])

    # resolution of the bins, in Angstrom
    settings.set_resolution(float(sys.argv[2]))

    try:
        df = pd.read_csv(settings.get_kabsch_aligned_csv_filename(), header=0)
        avg_frag = pd.read_csv(settings.outputfile_prefix + "_avg_fragment.csv", header=0)
    except FileNotFoundError:
        print('First align and calculate average fragment.')
        sys.exit(2)

    # grab only the atoms that are in the contact groups
    df_central = df[df['label'] == '-']
    coordinate_df = make_coordinate_df(df_central, settings, avg_frag)

    # find the volume of the central group
    tolerance = 0.5
    contact_group_radius = get_vdw_distance_contact(df, settings)
    volume = find_available_volume(avg_fragment=avg_frag, extra=(tolerance + contact_group_radius))
    print('Available volume:', volume)

    try:
        pd.read_hdf(settings.get_density_df_filename(), settings.get_density_df_key())
    except (FileNotFoundError, KeyError):
        empty_density_df = prepare_df(df=coordinate_df, settings=settings)

        density_df = count_points_per_square(df=empty_density_df, contact_points_df=coordinate_df, settings=settings)

        # save so we can use the data but only change the plot - saves time :)
        density_df.to_hdf(settings.get_density_df_filename(), settings.get_density_df_key())

    t1 = time.time() - t0
    print("Duration: %.2f s." % t1)


def count_points_per_square(df, contact_points_df, settings):
    contact_points_df = contact_points_df

    print("Counting points per bin: ")
    # prepare vector that will contain the amount
    amount = np.zeros(len(df))

    bin_coordinates = np.array([df.xstart, df.ystart, df.zstart])
    contact_coordinates = np.transpose(np.array([contact_points_df.x,
                                                 contact_points_df.y,
                                                 contact_points_df.z]))

    amount = fill_bins(amount, bin_coordinates, contact_coordinates, settings.resolution)

    df[settings.to_count_contact] = amount

    return df


@jit(nopython=True)
def fill_bins(amount, bin_coordinates, contact_coordinates, resolution):
    x, y, z = 0, 1, 2
    i = 0
    total = len(contact_coordinates)

    for i in range(total):
        cor = contact_coordinates[i]
        idx = np.where((bin_coordinates[x] <= cor[x]) & (bin_coordinates[x] + resolution >= cor[x]) &
                       (bin_coordinates[y] <= cor[y]) & (bin_coordinates[y] + resolution >= cor[y]) &
                       (bin_coordinates[z] <= cor[z]) & (bin_coordinates[z] + resolution >= cor[z]))

        amount[idx] += 1

        if i % 5000 == 0:
            print(i, "/", total)

    return amount


if __name__ == "__main__":
    main()
