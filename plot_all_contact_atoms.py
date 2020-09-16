# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This is a script that I wrote for my master thesis
# It loads the coordinates of the aligned fragments, and then plots the central
# group and all contact atoms/the centers of the contact groups around it.
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from helpers.density_helpers import count_points_per_square, plot_density, prepare_df, average_molecule, calculate_center

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

    # TODO: make this more general:
    to_count = ["O"]
    # to_count = ["center"]

    fragments_df = pd.read_csv(inputfilename, header=None)
    fragments_df.columns = ["entry_id", "fragment_id", "atom_label", "fragment_or_contact", "atom_x", "atom_y", "atom_z"]

    avg_fragment = average_molecule(fragments_df)

    Xs, Ys, Zs = [], [], []

    contact_group_df = fragments_df[fragments_df.fragment_or_contact == "f"]

    for entry_id in contact_group_df.entry_id.unique():
        entry_df = contact_group_df[contact_group_df.entry_id == entry_id]

        for fragment_id in entry_df.fragment_id.unique():
            fragment_df = entry_df[entry_df.fragment_id == fragment_id]
            
            for column in to_count:
                # if center, calculate per fragment instead of per atom
                if "center" == column:
                    x, y, z = calculate_center(fragment_df=fragment_df, atoms=["C"])
                    Xs.append(x)
                    Ys.append(y)
                    Zs.append(z)
                else:
                    point = fragment_df[fragment_df['atom_label'].str.contains(column)]

                    assert (len(point) == 1), " atom label is not unique, can't count per bin"

                    x, y, z = point.atom_x, point.atom_y, point.atom_z
                    Xs.append(x)
                    Ys.append(y)
                    Zs.append(z)

    make_plot(avg_fragment, Xs, Ys, Zs)


def make_plot(avg_fragment, Xs, Ys, Zs):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # plot the (average of the) central group 
    for atom in avg_fragment.atoms.values():
        if "O" in atom.label:
            ax.scatter(atom.x, atom.y, atom.z, c="red", s=20, edgecolor="black")
        elif "N" in atom.label:
            ax.scatter(atom.x, atom.y, atom.z, c="blue", s=20, edgecolor="black")
    
    ax.scatter(Xs, Ys, Zs, s=1, c="red")

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    plt.show("Central group and all contact atoms")
    plt.show()





if __name__ == "__main__":
    main()