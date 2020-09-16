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

from helpers.density_helpers import count_points_per_square, plot_density, prepare_df, average_molecule

import pandas as pd
import numpy as np

import sys
import time

def main():

    if len(sys.argv) != 3:
        print("Usage: python plot_density.py <path/to/inputfile> <resolution>")
        sys.exit(1)
    
    inputfilename = sys.argv[1]

    # resolution of the bins, in Angstrong
    resolution = float(sys.argv[2])

    

    intermediate_hdf_file = inputfilename.rsplit("/\\", 1)[-1].rsplit(".", 1)[0] + "_" + str(resolution) + ".hdf"
    plotname = inputfilename.rsplit("/\\", 1)[-1].rsplit(".", 1)[0] + "_" + str(resolution) + "_density.png"

    # TODO: make this more general:
    to_count = ["O"]
    # to_count = ["center"]

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

        density_df = count_points_per_square(df=empty_density_df, points_df=fragments_df, to_count=to_count)

        # save so we can use the data but only change the plot - saves time :)
        density_df.to_hdf(intermediate_hdf_file, 'key')

        print("It took me", time.time() - starttime, "s to calculate the df for a resolution of", resolution)

    plot_density(plotname, to_count, avg_fragment, density_df, resolution)


if __name__ == "__main__":
    main()