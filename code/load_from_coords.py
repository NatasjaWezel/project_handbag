# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This is a script that I wrote for my master thesis
# It loads the coordinates of the fragments exported from a conquest query and
# aligns the central groups by using rotation matrices and other linear algebra.
# It then saves the new coordinates in a .csv file.
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import csv
import sys
import os

import time

import numpy as np

import pandas as pd
from tqdm import tqdm


def main():

    if len(sys.argv) != 2:
        print("Usage: python load_from_coords.py <path/to/inputfile>")
        sys.exit(1)

    t0 = time.time()

    filename = sys.argv[1]

    lines = read_coord_file(filename)
    filename_csv = filename.rsplit(".", 1)[0] + ".csv"
    labels_df = read_parameter_csv(filename_csv)

    # alignment = {}
    # df = pd.read_csv(".\\data\\central_groups.csv")
    # df = df[df.name == "RC6H5"]

    # alignment['center'] = df.center_label.max()
    # alignment['yaxis'] = df.y_axis_label.max()
    # alignment['xyplane'] = df.xy_plane_label.max()

    # alignment['R'] = df.R.max()

    # if alignment["R"] == "-":
    #     alignment["R"] = None


    lines = load_fragments_from_coords(labels_df, lines)

    
    open("bla.csv", "w").writelines([line for line in lines])

    t1 = time.time() - t0

    print("duration: %.2f s." % t1)


def read_parameter_csv(filename):
    labels_contact_df = pd.read_csv(filename)

    # clean column names
    columns = labels_contact_df.columns
    columns = [i.strip() for i in columns]
    labels_contact_df.columns = columns

    labels_contact_df = labels_contact_df.drop(columns=["Query", "Refcode"])
    # labels_contact_df = labels_contact_df.transpose()
    return labels_contact_df


def read_coord_file(filename):
    """ Reads the file and saves its lines as a list. """

    with open(filename) as inputfile:
        lines = inputfile.readlines()

    return lines


def load_fragments_from_coords(labels_df, lines):
    """ Reads part of the coordinate file and returns an entire fragment. """

    print(labels_df)

    idx = -1
    row = None

    atom_labels = []
    new_lines = []

    print("Reading fragments from .cor file")
    for line in tqdm(lines):

        if "FRAG" in line:
            idx += 1
            row = labels_df[labels_df.index == idx].transpose().reset_index()
            row.columns = ["label", "atom"]

            atom_labels = []
            information = line.split("**")
            entry = information[0].strip()
            fragment_no = information[2].strip()

            # if idx == 1:
            #     sys.exit()
        else:
            information = line.split()
            atom_label = information[0].strip("%*")

            # print(row, atom_label)
            label = row[row["atom"] == atom_label]["label"].max()
            # if type(label) != float:
            #     print(label)

            atom_label = check_if_label_exists(atom_label, atom_labels)
            atom_labels.append(atom_label)

            x, y, z = information[1].split("(")[0], information[2].split("(")[0], information[3].split("(")[0]

            line = entry + fragment_no + "," + atom_label + "," + x + "," + y + "," + z + "," + '\n'
            new_lines.append(line)

    return new_lines


def check_if_label_exists(label, labels):
    """ Checks if label already exists in the fragment. Adds the letter 'a' to it if it does. """

    if label in labels:
        label += 'a'
        label = check_if_label_exists(label, labels)

    return label



if __name__ == "__main__":
    main()
