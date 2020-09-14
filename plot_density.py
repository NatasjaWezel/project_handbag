from helpers.density_helpers import count_points_per_square, plot_density, prepare_df, average_molecule

import pandas as pd
import numpy as np

import time 

def main():
    # TODO: if total isnt a cube, it doens make sense to have the same amount of bins
    # in each direction
    # Hardcoding these makes that the volumes aren't always the same in every plot
    bins_x = 9
    bins_y = 9
    bins_z = 5

    resultsdir = "results/no3_c6h5r/"
    inputfile = resultsdir + "coord_test_no3_c6h5r.csv"
    intermediate_hdf_file = resultsdir + "density_df_" + str(bins_x) + "-" + str(bins_y) + "-" + str(bins_z) + ".hdf"
    plotname = resultsdir + "density" + str(time.time()) + ".png"

    # TODO: make this more general:
    # to_count = ["O"]
    to_count = ["center"]

    fragments_df = pd.read_csv(inputfile, header=None)
    fragments_df.columns = ["entry_id", "fragment_id", "atom_label", "fragment_or_contact", "atom_x", "atom_y", "atom_z"]

    amount_bins = bins_x * bins_y * bins_z

    avg_fragment = average_molecule(fragments_df)

    try:
        density_df = pd.read_hdf(intermediate_hdf_file, 'no3_co')
    except FileNotFoundError:
        density_df = prepare_df(fragments_df=fragments_df, 
                                    amount_bins=amount_bins, 
                                    no_bins_x=bins_x, 
                                    no_bins_y=bins_y, 
                                    no_bins_z=bins_z, 
                                    to_count=to_count)

        density_df = count_points_per_square(df=density_df, points_df=fragments_df, to_count=to_count)

        # save so we can use the data but only change the plot - saves time :)
        density_df.to_hdf(intermediate_hdf_file, 'no3_co')

    plot_density(plotname, to_count, avg_fragment, density_df, amount_bins)


if __name__ == "__main__":
    main()