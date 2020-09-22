# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This is a script that I wrote for my master thesis
# It loads the coordinates of the aligned fragments, and then plots the central
# group and all contact atoms/the centers of the contact groups around it.
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from helpers.geometry_helpers import average_fragment, calculate_center
from helpers.plot_functions import plot_fragment_colored

import math

import pandas as pd
import numpy as np

import time 

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import sys

def main():

    if len(sys.argv) != 3:
        print("Usage: python plot_all_contact_atoms.py <path/to/inputfile> <atom to count or center>")
        sys.exit(1)
    
    inputfilename = sys.argv[1]
    to_count = sys.argv[2]

    prefix = inputfilename.rsplit("/\\", 1)[-1].rsplit(".", 1)[0] 
    avg_fragment_name = prefix + "_avg_fragment.pkl"

    fragments_df = pd.read_csv(inputfilename, header=None)
    fragments_df.columns = ["entry_id", "fragment_id", "atom_label", "atom_symbol", "fragment_or_contact", "atom_x", "atom_y", "atom_z"]

    avg_fragment = average_fragment(avg_fragment_name, fragments_df)

    coordinate_df = count_contact_atoms(fragments_df, to_count)

    coordinate_df = distances_closest_atom_central(coordinate_df, avg_fragment)

    make_plot(avg_fragment, coordinate_df)

def distances_closest_atom_central(coordinate_df, avg_fragment):
    closest_distances = []
    closest_atoms_vdw = []
    
    for x, y, z in zip(coordinate_df.x, coordinate_df.y, coordinate_df.z):
        closest_distance = math.inf

        for atom in avg_fragment.atoms.values():
            distance = np.sqrt((x - atom.x)**2 + (y - atom.y)**2 + (z - atom.z)**2)
            
            if distance < closest_distance:
                closest_distance = distance
                closest_atom_vdw = atom.vdw_radius
                
        closest_distances.append(closest_distance)
        closest_atoms_vdw.append(closest_atom_vdw)

    coordinate_df["distance"] = closest_distances
    coordinate_df["vdw_closest_atom"] = closest_atoms_vdw

    return coordinate_df

def count_contact_atoms(fragments_df, to_count):
    """ This is a function that counts the contact atoms or centers of contact groups near the
        central group. """

    contact_group_df = fragments_df[fragments_df.fragment_or_contact == "f"]

    labels = contact_group_df.entry_id.unique()
    coordinate_df = pd.DataFrame(columns=["x", "y", "z", "distance", "vdw_closest_atom"], index=labels)

    for entry_id in labels:
        entry_df = contact_group_df[contact_group_df.entry_id == entry_id]

        for fragment_id in entry_df.fragment_id.unique():
            fragment_df = entry_df[entry_df.fragment_id == fragment_id]
            
            coordinates = get_coordinates(fragment_df, to_count)
                
            coordinate_df.loc[coordinate_df.index == entry_id, ["x", "y", "z"]] = coordinates[0], coordinates[1], coordinates[2]

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
    
    ax.scatter(coordinate_df.x, coordinate_df.y, coordinate_df.z, s=1, c="red")

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    plt.show()


if __name__ == "__main__":
    main()