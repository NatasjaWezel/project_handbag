# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This is a script that I wrote for my master thesis
# It loads the coordinates of the aligned fragments, and then plots the central
# group and all contact atoms/the centers of the contact groups around it.
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from helpers.geometry_helpers import average_fragment, calculate_center
from helpers.plot_functions import plot_fragment_colored, plot_vdw_spheres
from helpers.helpers import read_results_alignment

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
    
    inputfilename = sys.argv[1]
    to_count = sys.argv[2]

    prefix = inputfilename.rsplit("/\\", 1)[-1].rsplit(".", 1)[0] 
    avg_fragment_name = prefix + "_avg_fragment.pkl"

    aligned_fragments_df = read_results_alignment(inputfilename)

    avg_fragment = average_fragment(avg_fragment_name, aligned_fragments_df)

    coordinate_df = count_contact_atoms(aligned_fragments_df, to_count)

    # TODO: this only works for single atom, not for center
    radii_df = pd.read_csv(RADII_CSV)

    vdw_radius = float(radii_df[radii_df.symbol == to_count].radius)

    coordinate_df["vdw"] = vdw_radius

    coordinate_df = distances_closest_vdw_central(coordinate_df, avg_fragment)

    make_plot(avg_fragment, coordinate_df)

def distances_closest_vdw_central(coordinate_df, avg_fragment):
    closest_distances = []
    closest_atoms_vdw = []
    
    for x, y, z, vdw in zip(coordinate_df.x, coordinate_df.y, coordinate_df.z, coordinate_df.vdw):
        closest_distance = math.inf

        for atom in avg_fragment.atoms.values():
            distance = np.sqrt((x - atom.x)**2 + (y - atom.y)**2 + (z - atom.z)**2) - vdw - atom.vdw_radius
            
            if distance < closest_distance:
                closest_distance = distance
                closest_atom_vdw = atom.vdw_radius
                
        # assert closest_distance <5, "closest distance can't be smaller then 3"

        closest_distances.append(closest_distance)
        closest_atoms_vdw.append(closest_atom_vdw)

    coordinate_df["distance"] = closest_distances
    coordinate_df["vdw_closest_atom"] = closest_atoms_vdw

    return coordinate_df

def count_contact_atoms(fragments_df, to_count):
    """ This is a function that counts the contact atoms or centers of contact groups near the
        central group. """

    contact_group_df = fragments_df[fragments_df.fragment_or_contact == "f"]

    unique_fragments = contact_group_df.unique_fragment
    coordinate_df = pd.DataFrame(columns=["x", "y", "z"], index=unique_fragments)

    for unique_fragment_id in unique_fragments:
        fragment_df = contact_group_df[contact_group_df.unique_fragment == unique_fragment_id]
            
        coordinates = get_coordinates(fragment_df, to_count)
            
        coordinate_df.loc[coordinate_df.index == unique_fragment_id, ["x", "y", "z"]] = coordinates[0], coordinates[1], coordinates[2]

    return coordinate_df

def get_coordinates(fragment_df, to_count):
    # if center, calculate per fragment instead of per atom
    if to_count == "center":
        return calculate_center(fragment_df=fragment_df, atoms=["C"])
    else:
        point = fragment_df[fragment_df['atom_label'].str.contains(to_count)]

        assert (len(point) == 1), " atom label is not unique, can't count per bin"
        
        return [float(point.atom_x), float(point.atom_y), float(point.atom_z)]

def make_plot(avg_fragment, coordinate_df):
    """ Plot all the surrounding contact groups around the central group. """

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # plot the (average of the) central group 
    ax = plot_fragment_colored(ax, avg_fragment)
    ax, _ = plot_vdw_spheres(avg_fragment, ax)

    points = ax.scatter(list(coordinate_df.x), list(coordinate_df.y), list(coordinate_df.z), s=1, c="red")

    print(coordinate_df.distance.min(), coordinate_df.distance.max())

    vdw_slider_ax = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=AXCOLOR)
    vdw_slider = Slider(vdw_slider_ax, 'VDW radius + ', -0.5, 2, valinit=0, valstep=0.1)

    def update(val):
        val = round(val, 2)
       
        show_df = coordinate_df[coordinate_df.distance <= 2 * val]

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