

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
from helpers.geometry_helpers import make_avg_fragment_if_not_exists, calculate_center, calculate_longest_vdw_radius_contact
from helpers.helpers import read_results_alignment

import math

import pandas as pd
import numpy as np

import sys
import time

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from classes.Settings import Settings

from tqdm import tqdm

import csv

def main():

    if len(sys.argv) != 4:
        print("Usage: python analyze_density.py <path/to/inputfile> <resolution> <atom or center to count>")
        sys.exit(1)
    
    settings = Settings(sys.argv[1])
    settings.set_central_group()

    # resolution of the bins, in Angstrom
    settings.set_resolution(float(sys.argv[2]))

    aligned_fragments_df = read_results_alignment(settings.get_aligned_csv_filename())

    try:
        density_df = pd.read_hdf(settings.get_density_df_filename(), settings.get_density_df_key())
    except (FileNotFoundError, KeyError):
        print("Run calc_density first")
        sys.exit(1)


    calculate_80_percent(density_df, settings)
    make_plot(first_fragment_df, density_df, settings)

def make_plot(avg_fragment, density_df, settings):
    resolution = settings.resolution
    plotname = settings.get_density_plotname()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax = plot_fragment_colored(ax, avg_fragment)

    p, ax = plot_density(ax=ax, df=density_df, resolution=resolution)

    ax.set_title("4D density plot\n Resolution: " + str(resolution))

    fig.colorbar(p)
    plt.savefig(plotname)
    plt.show()


def calculate_80_percent(df, settings):
    column_name = "amount_" + settings.to_count_contact
    total_atoms = df[column_name].sum()

    # take bins that aren't empty and sort them from highest fraction to lowest fraction
    population = list(df[df[column_name] > 0][column_name])
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

        writer.writerow([settings.central_group_name, settings.contact_group_name, settings.to_count_contact, fraction, i, non_empty_bins,])

    print(str(round(fraction * 100, 2)) + "% of the contact group atoms is in " + str(round(i/non_empty_bins*100, 2)) + "% of the non-empty bins")
    print(str(round(fraction * 100, 2)) + "% of the contact group atoms is in " + str(round(i/bins*100, 2)) + "% of the total bins")


if __name__ == "__main__":
    main()