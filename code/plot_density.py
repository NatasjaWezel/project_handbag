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

from helpers.density_helpers import prepare_df, add_one_to_bin
from helpers.plot_functions import plot_fragment_colored, plot_density
from helpers.geometry_helpers import average_molecule, calculate_center

import pandas as pd
import numpy as np

import sys
import time

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def main():

    if len(sys.argv) != 3:
        print("Usage: python plot_density.py <path/to/inputfile> <resolution>")
        sys.exit(1)
    
    inputfilename = sys.argv[1]

    # TODO: make this more general
    to_count = ["center"]

    # resolution of the bins, in Angstrong
    resolution = float(sys.argv[2])

    intermediate_hdf_file = inputfilename.rsplit("/\\", 1)[-1].rsplit(".", 1)[0] + "_" + str(resolution) + ".hdf"
    plotname = inputfilename.rsplit("/\\", 1)[-1].rsplit(".", 1)[0] + "_" + str(resolution) + "_density.png"

    fragments_df = pd.read_csv(inputfilename, header=None)
    fragments_df.columns = ["entry_id", "fragment_id", "atom_label", "fragment_or_contact", "atom_x", "atom_y", "atom_z"]

    avg_fragment = average_molecule(fragments_df)

    try:
        density_df = pd.read_hdf(intermediate_hdf_file, 'key')
    except FileNotFoundError:
        starttime = time.time()

        empty_density_df = prepare_df(fragments_df=fragments_df, 
                                    resolution=resolution, 
                                    to_count=to_count)

        density_df = count_points_per_square(df=empty_density_df, points_df=fragments_df)

        # save so we can use the data but only change the plot - saves time :)
        density_df.to_hdf(intermediate_hdf_file, 'key')

        print("It took me", time.time() - starttime, "s to calculate the df for a resolution of", resolution)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax = plot_fragment_colored(ax, avg_fragment)
    p, ax = plot_density(ax=ax, df=density_df)

    ax.set_title("4D density plot\n Resolution: " + str(resolution))

    fig.colorbar(p)
    plt.savefig(plotname)
    plt.show()


def count_points_per_square(df, points_df):
    columns = list(df.columns)
    column_name = [i.split("_") for i in columns if "amount" in i]

    small_points_df = points_df[points_df.fragment_or_contact == "f"]
    unique_entries = small_points_df.entry_id.unique()
    total_entries = len(unique_entries)

    for i, entry_id in enumerate(unique_entries):
        entry_df = small_points_df[small_points_df.entry_id == entry_id]

        for fragment_id in entry_df.fragment_id.unique():
            fragment_df = entry_df[entry_df.fragment_id == fragment_id]
                           
            # if center, calculate per fragment instead of per atom
            if "center" in column_name:
                # TODO: can it say just C here?
                coordinates = calculate_center(fragment_df=fragment_df, atoms=["C"])

            else:
                point = fragment_df[fragment_df.atom_label.str.contains(column_name.split("_")[1])]

                assert (len(point) == 1), " atom label is not unique, can't count per bin"

                coordinates = [float(point.atom_x), float(point.atom_y), float(point.atom_z)]
                
            df = add_one_to_bin(df, column_name, coordinates)
        
        if i % 100 == 0:
            print(str(i) + "/" + str(total_entries) + " done")

    # TODO: see if counting is right
    # test_count(df, small_points_df)
    return df


if __name__ == "__main__":
    main()