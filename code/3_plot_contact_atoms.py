# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This is a script that I wrote for my master thesis
# It loads the coordinates of the aligned fragments, and then plots the central
# group and all contact atoms/the centers of the contact groups around it.
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from helpers.geometry_helpers import make_avg_fragment_if_not_exists, calculate_center, calculate_longest_vdw_radius_contact
from helpers.plot_functions import plot_fragment_colored, plot_vdw_spheres
from helpers.helpers import read_results_alignment

from classes.Settings import Settings

from matplotlib.widgets import Slider

import math

import pandas as pd
import numpy as np

import time 

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from helpers.headers import AXCOLOR, RADII_CSV

import sys

from tqdm import tqdm

def main():

    if len(sys.argv) != 3:
        print("Usage: python plot_all_contact_atoms.py <path/to/inputfile> <atom/center to count>")
        sys.exit(1)
    
    settings = Settings(sys.argv[1])
    settings.set_atom_to_count(sys.argv[2])

    df = read_results_alignment(settings.get_aligned_csv_filename())

    avg_fragment = make_avg_fragment_if_not_exists(settings, df)

    # grab only the atoms that are in the contact groups
    df = df[df.in_central_group == False]

    first_fragment_df = df[df.id == df.id.unique()[0]]

    if settings.to_count_contact == "centroid":
        # plot centroids of all contact fragments
        coordinate_df = df.groupby("id").mean()
        vdw_distance_contact = calculate_longest_vdw_radius_contact(first_fragment_df, settings)
    elif len(first_fragment_df[first_fragment_df["atom_symbol"] == settings.to_count_contact]) == 1:
        # atom is unique, plot all of them
        coordinate_df = df[df.atom_symbol == settings.to_count_contact]
        vdw_distance_contact = settings.get_vdw_radius(settings.to_count_contact)
    else:
        # TODO: atom is not unique, find closest
        pass

    coordinate_df = make_coordinate_df(coordinate_df, settings, avg_fragment)

    make_plot(avg_fragment, coordinate_df, vdw_distance_contact)


def make_coordinate_df(df, settings, avg_fragment):

    try:
        coordinate_df = pd.read_hdf(settings.get_coordinate_df_filename(), settings.get_coordinate_df_key())

        return coordinate_df
    except FileNotFoundError:
        coordinate_df = distances_closest_vdw_central(df, avg_fragment, settings)

        coordinate_df.to_hdf(settings.get_coordinate_df_filename(), settings.get_coordinate_df_key())

    return coordinate_df


def distances_closest_vdw_central(coordinate_df, avg_fragment, settings):
    closest_distances = []
    closest_atoms_vdw = []

    points_avg_f = np.array([avg_fragment.atom_x, avg_fragment.atom_y, avg_fragment.atom_z]).T
    
    print("Searching for nearest atom from contact group...")
    print(len(coordinate_df.atom_x))
    for x, y, z in tqdm(zip(coordinate_df.atom_x, coordinate_df.atom_y, coordinate_df.atom_z)):
        
        p2 = np.array([x,y,z])

        dist = np.sqrt([np.sum((f - p2)**2, axis=0) for f in points_avg_f])
        
        min_dist_idx = dist.argmin()
        min_dist = dist[min_dist_idx]
        
        min_atom_vdw = avg_fragment.iloc[min_dist_idx]['vdw_radius']

        closest_distances.append(min_dist)
        closest_atoms_vdw.append(min_atom_vdw)

    coordinate_df["distance"] = closest_distances
    coordinate_df["vdw_closest_atom"] = closest_atoms_vdw

    return coordinate_df


def make_plot(avg_fragment, coordinate_df, longest_vdw_contact):
    """ Plot all the surrounding contact groups around the central group. """

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # plot the (average of the) central group 
    ax = plot_fragment_colored(ax, avg_fragment)
    ax, _ = plot_vdw_spheres(avg_fragment, ax, color='pink')

    points = ax.scatter(list(coordinate_df[coordinate_df.distance <= coordinate_df.vdw_closest_atom + longest_vdw_contact].atom_x), 
                        list(coordinate_df[coordinate_df.distance <= coordinate_df.vdw_closest_atom + longest_vdw_contact].atom_y), 
                        list(coordinate_df[coordinate_df.distance <= coordinate_df.vdw_closest_atom + longest_vdw_contact].atom_z), s=1, c="red")

    vdw_slider_ax = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=AXCOLOR)
    vdw_slider = Slider(vdw_slider_ax, 'VDW radius + ', -2, 3, valinit=0, valstep=0.1)

    def update(val):
        val = round(val, 2)
       
        show_df = coordinate_df[coordinate_df.distance <= coordinate_df.vdw_closest_atom + longest_vdw_contact + val]

        print(len(coordinate_df) , len(show_df), val)    
        points._offsets3d = (list(show_df.atom_x), list(show_df.atom_y), list(show_df.atom_z))
        fig.canvas.draw_idle()

    vdw_slider.on_changed(update)
    
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    plt.show()


if __name__ == "__main__":
    main()