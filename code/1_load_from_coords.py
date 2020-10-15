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

import pandas as pd
from tqdm import tqdm

from classes.Atom import Atom
from classes.Fragment import Fragment
from classes.Settings import Settings
from helpers.helpers import check_if_label_exists
from helpers.rotation_helpers import perform_rotations


def main():

    if len(sys.argv) != 2:
        print("Usage: python load_from_coords.py <path/to/inputfile>")
        sys.exit(1)

    filename = sys.argv[1]

    settings = Settings(filename)

    coordinate_lines = read_coord_file(filename=filename)
    fragments = load_fragments_from_coords(coordinate_lines)

    labels_contact_df = pd.read_csv(settings.parameter_csv)

    # clean column names
    columns = labels_contact_df.columns
    columns = [i.strip() for i in columns]
    labels_contact_df.columns = columns

    labels = [i for i in columns if "LAB" in i]

    alignment_labels = settings.alignment_labels()
    print(labels, alignment_labels)

    outputfile = open(settings.get_aligned_csv_filename(), 'w', newline='')
    writer = csv.writer(outputfile)

    print("Aligning fragments and writing result to csv")
    for i, fragment in enumerate(tqdm(fragments)):
        # get contact group and relabel it
        row = labels_contact_df[labels_contact_df.index == i]

        for j, label in enumerate(labels):
            # use .max to get the string out of the row
            atom_label = row[label].max()
            atom = fragment.atoms[atom_label]
            atom.add_to_central_group()

            if label == alignment_labels['center']:
                center_atom = atom.label
            elif label == alignment_labels['yaxis']:
                y_axis_atom = atom.label
            elif label == alignment_labels['xyplane']:
                xy_plane_atom = atom.label

            # TODO: maybe check if this label exists in the central group........
            if label == alignment_labels["R"]:
                atom.label = "R" + str(j + 1)
            else:
                atom.label = atom.symbol + str(j + 1)

        # center coordinates on center and rotate it
        fragment.center_coordinates(center_atom)
        fragment = perform_rotations(fragment, [y_axis_atom, xy_plane_atom])

        fragment.invert_if_neccessary()

        write_fragment_to_csv(writer, fragment)

    outputfile.close()


def write_fragment_to_csv(writer, fragment):
    """ This function saves the information of the fragment to a CSV file. """

    [writer.writerow([fragment.id, fragment.from_entry, atom.label, atom.symbol, atom.in_central_group, atom.x,
                      atom.y, atom.z]) for atom in fragment.atoms.values()]


def read_coord_file(filename):
    """ Reads the file and saves its lines as a list. """

    with open(filename) as inputfile:
        lines = inputfile.readlines()

    return lines


def load_fragments_from_coords(lines):
    """ Reads part of the coordinate file and returns an entire fragment. """

    fragments = []
    fragment = None

    print("Reading fragments from .cor file")
    for line in tqdm(lines):

        if "FRAG" in line and fragment is None:
            information = line.split("**")
            fragment = Fragment(fragment_id=information[2].strip(), from_entry=information[0].strip())
        elif "FRAG" in line:
            fragments.append(fragment)

            # if we found the header of the next fragment
            information = line.split("**")
            fragment = Fragment(fragment_id=information[2].strip(), from_entry=information[0].strip())

        else:
            information = line.split()
            x, y, z = information[1].split("("), information[2].split("("), information[3].split("(")

            atom = Atom(label=information[0].strip("%"), coordinates=[float(x[0]), float(y[0]), float(z[0])])

            atom = check_if_label_exists(atom, fragment)

            fragment.add_atom(atom)

    fragments.append(fragment)

    # this return is here for the last fragment
    return fragments


if __name__ == "__main__":
    main()
