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
from helpers.geometry_helpers import average_fragment, calculate_center
from helpers.helpers import read_results_alignment

import pandas as pd
import numpy as np

import sys
import time

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from classes.Settings import Settings

from tqdm import tqdm

def main():

    if len(sys.argv) != 4:
        print("Usage: python plot_density.py <path/to/inputfile> <resolution> <atom or center to count>")
        sys.exit(1)
    
    inputfilename = sys.argv[1]
    
    # resolution of the bins, in Angstrong
    resolution = float(sys.argv[2])

    to_count = sys.argv[3]

    settings = Settings(inputfilename)

    aligned_fragments_df = read_results_alignment(settings.get_aligned_csv_filename())
    
    avg_fragment = average_fragment(settings.get_avg_fragment_filename(), aligned_fragments_df)

    try:
        density_df = pd.read_hdf(settings.get_density_df_filename(resolution), settings.get_density_df_key(resolution))
    except FileNotFoundError:
        starttime = time.time()

        empty_density_df = prepare_df(fragments_df=aligned_fragments_df, resolution=resolution)

        empty_density_df["amount_" + to_count] = 0

        density_df = count_points_per_square(df=empty_density_df, points_df=aligned_fragments_df, resolution=resolution)

        # save so we can use the data but only change the plot - saves time :)
        density_df.to_hdf(settings.get_density_df_filename(resolution), settings.get_density_df_key(resolution))

        print("It took me", time.time() - starttime, "s to calculate the df for a resolution of", resolution)


    calculate_80_percent(density_df, to_count)
    make_plot(avg_fragment, density_df, resolution, settings.get_density_plotname(resolution))

def make_plot(avg_fragment, density_df, resolution, plotname):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax = plot_fragment_colored(ax, avg_fragment)
    p, ax = plot_density(ax=ax, df=density_df, resolution=resolution)

    ax.set_title("4D density plot\n Resolution: " + str(resolution))

    fig.colorbar(p)
    plt.savefig(plotname)
    plt.show()


def count_points_per_square(df, points_df, resolution):
    columns = list(df.columns)

    # TODO: check out this [0]
    column_name = [i for i in columns if "amount" in i][0]

    small_points_df = points_df[points_df.in_central_group == False]
    unique_fragments = small_points_df.id.unique()

    print("Counting points per bin: ")
    for fragment_id in tqdm(unique_fragments):
        fragment_df = small_points_df[small_points_df.id == fragment_id]
       
        # if center, calculate per fragment instead of per atom
        if "center" in column_name:
            # TODO: can it say just C here?
            coordinates = calculate_center(fragment_df=fragment_df, atoms=["C"])

        else:
            point = fragment_df[fragment_df.atom_label.str.contains(column_name.split("_")[1])]
           
            assert (len(point) == 1), " atom label is not unique, can't count per bin"

            coordinates = [float(point.atom_x), float(point.atom_y), float(point.atom_z)]
            
        df = add_one_to_bin(df=df, column_name=column_name, resolution=resolution, coordinates=coordinates)
        
    return df

def calculate_80_percent(df, to_count):
    total_atoms = df["amount_" + to_count].sum()
    population = list(df[df["amount_" + to_count] > 0]["amount_" + to_count])
    population.sort(reverse=True)

    fraction = 0
    i = 0

    while fraction < 0.8:
        fraction += population[i]/total_atoms
        i+=1
    
    non_empty_bins = len(population)
    bins = len(df)
    print(str(round(fraction * 100, 2)) + "% of the contact group atoms is in " + str(round(i/non_empty_bins*100, 2)) + "% of the non-empty bins")
    print(str(round(fraction * 100, 2)) + "% of the contact group atoms is in " + str(round(i/bins*100, 2)) + "% of the total bins")


if __name__ == "__main__":
    main()