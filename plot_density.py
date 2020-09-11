from helpers.density_helpers import find_maximum, count_points_per_square, plot_density, prepare_df

import pandas as pd
import numpy as np

def main():
    plotname = "density.png"

    df = pd.read_csv("results/NO3_CO_vdw5.csv", header=None)
    df.columns = ["entry_id", "fragment_id", "atom_label", "atom_x", "atom_y", "atom_z"]
    print(df.head())

    maxx, maxy, maxz = find_maximum(df)
    minx, miny, minz = -maxx, -maxy, -maxz

    # TODO: if total isnt a cube, it doens make sense to have the same amount of bins
    # in each direction
    amount_bins = 10
    bins = [np.linspace(minx, maxx, num=amount_bins+1),
                np.linspace(miny, maxy, num=amount_bins+1),
                np.linspace(minz, maxz, num=amount_bins+1)] 

    indices = [i for i in range(0, amount_bins * amount_bins * amount_bins)]

    density_df = prepare_df(amount_bins=amount_bins, bins=bins, indices=indices)
    density_df = count_points_per_square(df=density_df, points_df=df)
    density_df.to_hdf("test2.hdf", 'no3_co')

    # density_df = pd.read_hdf("test2.hdf", 'no3_co')
    print(density_df.head())

    plot_density(plotname, density_df, amount_bins, minx, maxx, miny, maxy, minz, maxz)


if __name__ == "__main__":
    main()