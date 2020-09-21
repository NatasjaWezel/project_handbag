# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This is a script that I wrote for my master thesis
# It loads the coordinates of the aligned fragments, and then plots the central
# group and all contact atoms/the centers of the contact groups around it.
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from helpers.geometry_helpers import average_fragment, calculate_center

import pandas as pd
import numpy as np

import time 

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import sys

def main():

    if len(sys.argv) != 2:
        print("Usage: python plot_all_contact_atoms.py <path/to/inputfile>")
        sys.exit(1)
    
    inputfilename = sys.argv[1]

    prefix = inputfilename.rsplit("/\\", 1)[-1].rsplit(".", 1)[0] 
    avg_fragment_name = prefix + "_avg_fragment.pkl"

    # TODO: make this more general:
    # to_count = ["O"]
    to_count = ["center"]

    fragments_df = pd.read_csv(inputfilename, header=None)
    fragments_df.columns = ["entry_id", "fragment_id", "atom_label", "atom_symbol", "fragment_or_contact", "atom_x", "atom_y", "atom_z"]

    avg_fragment = average_fragment(avg_fragment_name, fragments_df)

    coordinate_lists = count_contact_atoms(fragments_df, to_count)

    make_plot(avg_fragment, coordinate_lists)

def count_contact_atoms(fragments_df, to_count):
    """ This is a function that counts the contact atoms or centers of contact groups near the
        central group. """

    Xs, Ys, Zs = [], [], []

    contact_group_df = fragments_df[fragments_df.fragment_or_contact == "f"]

    for entry_id in contact_group_df.entry_id.unique():
        entry_df = contact_group_df[contact_group_df.entry_id == entry_id]

        for fragment_id in entry_df.fragment_id.unique():
            fragment_df = entry_df[entry_df.fragment_id == fragment_id]
            
            # TODO: this will break if to_count is more than a single atom
            for column in to_count:
                # if center, calculate per fragment instead of per atom
                if "center" == column:
                    coordinates = calculate_center(fragment_df=fragment_df, atoms=["C"])
                else:
                    point = fragment_df[fragment_df['atom_label'].str.contains(column)]

                    assert (len(point) == 1), " atom label is not unique, can't count per bin"

                    coordinates = [point.atom_x, point.atom_y, point.atom_z]

                Xs.append(coordinates[0])
                Ys.append(coordinates[1])
                Zs.append(coordinates[2])

    return [Xs, Ys, Zs]


def make_plot(avg_fragment, coordinate_lists):
    """ Plot all the surrounding contact groups around the central group. """

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # plot the (average of the) central group 
    for atom in avg_fragment.atoms.values():
        if "O" in atom.label:
            ax.scatter(atom.x, atom.y, atom.z, c="red", s=20, edgecolor="black")
        elif "N" in atom.label:
            ax.scatter(atom.x, atom.y, atom.z, c="blue", s=20, edgecolor="black")
    
    ax.scatter(coordinate_lists[0], coordinate_lists[1], coordinate_lists[2], s=1, c="red")

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    plt.title("Central group and all contact atoms")
    plt.show()


if __name__ == "__main__":
    main()