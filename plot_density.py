from density_helpers import find_maximum, count_points_per_square, plot_density, prepare_df

import pandas as pd
import numpy as np

def main():
    df = pd.read_csv("results/NO3_CO_vdw5.csv", header=None)
    df.columns = ["entry_id", "fragment_id", "atom_label", "atom_x", "atom_y", "atom_z"]
    print(df.head())

    maximum = find_maximum(df)
    minimum = -maximum

    amount_bins = 10
    bins = np.linspace(minimum, maximum, num=amount_bins+1)
    indices = [i for i in range(0, amount_bins * amount_bins * amount_bins)]

    # density_df = prepare_df(amount_bins=amount_bins, bins=bins, indices=indices)
    # density_df = count_points_per_square(df=density_df, points_df=df)
    # density_df.to_hdf("test2.hdf", 'no3_co')

    density_df = pd.read_hdf("test2.hdf", 'no3_co')
    print(density_df.head())

    plot_density(density_df, amount_bins, minimum, maximum)


if __name__ == "__main__":
    main()