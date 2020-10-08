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

def main():

    if len(sys.argv) != 3:
        print("Usage: python plot_all_contact_atoms.py <path/to/inputfile> <atom/center to count>")
        sys.exit(1)
    
    settings = Settings(sys.argv[1])
    settings.set_central_group()
    settings.set_atom_to_count(sys.argv[2])

    aligned_fragments_df = read_results_alignment(settings.get_aligned_csv_filename())
    avg_fragment = make_avg_fragment_if_not_exists(settings, aligned_fragments_df)

    for atom in avg_fragment.atoms.values():
        atom.set_vdw_radius(settings.get_vdw_radius(atom.symbol))

    coordinate_df = count_contact_atoms(aligned_fragments_df, settings)
    coordinate_df = distances_closest_vdw_central(coordinate_df, avg_fragment, settings)

    first_fragment_df = aligned_fragments_df[aligned_fragments_df.id == aligned_fragments_df.id.unique()[0]]
    vdw_distance_contact = calculate_longest_vdw_radius_contact(first_fragment_df, settings)

    make_plot(avg_fragment, coordinate_df, vdw_distance_contact)

def distances_closest_vdw_central(coordinate_df, avg_fragment, settings):
    closest_distances = []
    closest_atoms_vdw = []
    
    for x, y, z in zip(coordinate_df.x, coordinate_df.y, coordinate_df.z):
        closest_distance = math.inf

        for atom in avg_fragment.atoms.values():
            distance = np.sqrt((x - atom.x)**2 + (y - atom.y)**2 + (z - atom.z)**2)
            
            if distance < closest_distance:
                closest_distance = distance
                closest_atom_vdw = settings.get_vdw_radius(atom.symbol)

        # assert closest_distance <5, "closest distance can't be smaller then 3"

        closest_distances.append(closest_distance)
        closest_atoms_vdw.append(closest_atom_vdw)

    coordinate_df["distance"] = closest_distances
    coordinate_df["vdw_closest_atom"] = closest_atoms_vdw

    return coordinate_df

def count_contact_atoms(fragments_df, settings):
    """ This is a function that counts the contact atoms or centers of contact groups near the
        central group. """

    contact_group_df = fragments_df[fragments_df.in_central_group == False]
    ids = contact_group_df.id.unique()

    coordinate_df = pd.DataFrame(columns=["x", "y", "z"], index=ids)

    for _id in ids:
        fragment_df = contact_group_df[contact_group_df.id == _id]
            
        coordinates = get_coordinates(fragment_df, settings)
        coordinate_df.loc[coordinate_df.index == _id, ["x", "y", "z"]] = coordinates[0], coordinates[1], coordinates[2]

    return coordinate_df

def get_coordinates(fragment_df, settings):
    # if center, calculate per fragment instead of per atom
    if settings.to_count_contact == "center":
        return calculate_center(fragment_df=fragment_df)
    else:
        point = fragment_df[fragment_df.atom_symbol == settings.to_count_contact]

        assert (len(point) == 1), " atom label is not unique, can't count per bin"
        
        return [float(point.atom_x), float(point.atom_y), float(point.atom_z)]

def make_plot(avg_fragment, coordinate_df, longest_vdw_contact):
    """ Plot all the surrounding contact groups around the central group. """

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # plot the (average of the) central group 
    ax = plot_fragment_colored(ax, avg_fragment)
    ax, _ = plot_vdw_spheres(avg_fragment, ax, color='pink')

    points = ax.scatter(list(coordinate_df[coordinate_df.distance <= coordinate_df.vdw_closest_atom + longest_vdw_contact].x), 
                        list(coordinate_df[coordinate_df.distance <= coordinate_df.vdw_closest_atom + longest_vdw_contact].y), 
                        list(coordinate_df[coordinate_df.distance <= coordinate_df.vdw_closest_atom + longest_vdw_contact].z), s=1, c="red")

    vdw_slider_ax = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=AXCOLOR)
    vdw_slider = Slider(vdw_slider_ax, 'VDW radius + ', -2, 3, valinit=0, valstep=0.1)

    def update(val):
        val = round(val, 2)
       
        show_df = coordinate_df[coordinate_df.distance <= coordinate_df.vdw_closest_atom + longest_vdw_contact + val]

        print(len(coordinate_df) , len(show_df), val)    
        points._offsets3d = (list(show_df.x), list(show_df.y), list(show_df.z))
        fig.canvas.draw_idle()

    vdw_slider.on_changed(update)
    
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    plt.show()


if __name__ == "__main__":
    main()