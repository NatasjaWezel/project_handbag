# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This is a script that I wrote for my master thesis
# It loads the coordinates of the fragments exported from a conquest query and
# aligns the central groups by using rotation matrices and other linear algebra.
# It then saves the new coordinates in a .csv file.
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import sys

import time

import numpy as np

import pandas as pd
from tqdm import tqdm

from helpers.headers import CENTRAL_GROUPS_CSV


def main():

    if len(sys.argv) != 2:
        print("Usage: python load_from_coords.py <path/to/inputfile>")
        sys.exit(1)

    t0 = time.time()

    filename = sys.argv[1]

    # read part of the datafile first for preparation
    no_atoms, no_atoms_central, label_list = read_coord_file(filename)

    print("Pandas is reading csv...")
    # use amount of atoms for skiprow function to read into df immediately
    data = pd.read_csv(filename, sep='\\s+', usecols=[0, 1, 2, 3], names=['_id', 'x', 'y', 'z'],
                       header=None,
                       skiprows=lambda x: logic(x, no_atoms))
    print("Done")

    data['symbol'] = data.apply(lambda row: get_atom_symbol(row), axis=1)

    amount_rows = len(data)

    # calc amount of fragments
    fragments = int(amount_rows / no_atoms)

    # give each fragment a unique id
    fragment_ids = range(0, fragments)
    fragment_ids = np.repeat(fragment_ids, no_atoms)
    data['fragment_id'] = fragment_ids

    # give each atom in each fragment their label from conquest
    labels = label_list * fragments
    data['label'] = labels

    alignment = alignment_dict(central_group_name=filename.rsplit('\\')[-1].rsplit('.', 1)[0]
                               .rsplit('_aligned', 1)[0].split("_")[0])

    data_xyz_matrix = np.array([np.array(data.x), np.array(data.y), np.array(data.z)]).T

    # TODO: give R a label

    print("Rotating all fragments...")
    for i in tqdm(range(fragments)):
        # translate and rotate first fragment onto the origin as for nice viewing
        begin, end = i * no_atoms, (i + 1) * no_atoms
        part_matrix = data_xyz_matrix[begin:end].copy()

        # give the index of which atom is needed for the calculations
        part_matrix = perform_translation(part_matrix, index_center=label_list.index(alignment['center']))

        data_xyz_matrix[begin:end] = perform_rotations(part_matrix.copy(),
                                                       [label_list.index(alignment['yaxis']),
                                                        label_list.index(alignment['xyplane'])])

        # TODO: invert?

        # calculate and save error
        # TODO

    # put back into df
    data_xyz_matrix_T = data_xyz_matrix.T
    x_vec, y_vec, z_vec = data_xyz_matrix_T[0], data_xyz_matrix_T[1], data_xyz_matrix_T[2]
    data.x, data.y, data.z = x_vec, y_vec, z_vec

    # reindex the data to a more readable format
    data = data.reindex(['fragment_id', '_id', 'symbol', 'label', 'x', 'y', 'z'], axis=1)

    # save as csv
    data.to_csv('test_rotations.csv', index=False)

    t1 = time.time() - t0

    print("Duration: %.2f s." % t1)


def alignment_dict(central_group_name):
    alignment = {}
    df = pd.read_csv(CENTRAL_GROUPS_CSV)
    df = df[df.name == central_group_name]

    alignment['center'] = df.center_label.max()
    alignment['yaxis'] = df.y_axis_label.max()
    alignment['xyplane'] = df.xy_plane_label.max()

    alignment['R'] = df.R.max()

    if alignment["R"] == "-":
        alignment["R"] = None

    return alignment


def get_atom_symbol(row):
    if row._id[1] in '0123456789':
        return row._id[:1]

    return row._id[:2]


def logic(index, amount_atoms):
    if index % (amount_atoms+1) == 0:
        return True
    return False


def get_translation_dict(filename):

    with open(filename) as inputfile:
        lines = [next(inputfile) for x in range(2)]

    translation_dict = {}

    labels = lines[0].split(',')
    atoms = lines[1].split(',')

    for label, atom in zip(labels, atoms):
        if "LAB" in label:
            translation_dict[atom.strip()] = label.strip()

    return translation_dict


def read_coord_file(filename):
    """ Reads the first 100 lines of a file and saves them as a list. """

    csv_filename = filename.rsplit('.', 1)[0] + '.csv'
    translation_dict = get_translation_dict(csv_filename)

    with open(filename) as inputfile:
        lines = [next(inputfile) for x in range(100)]

    label_list = []

    # count number of atoms in a fragment
    no_atoms = 0
    no_atoms_central = 0

    for line in lines[1:]:
        if "FRAG" in line:
            break

        label_list.append(translation_dict.get(line.split(' ')[0], '-'))
        no_atoms += 1

        if translation_dict.get(line.split(' ')[0], '-') != "-":
            no_atoms_central += 1

    return no_atoms, no_atoms_central, label_list


def perform_translation(fragment, index_center):
    translation_vec = np.array([-float(fragment[index_center][0]),
                                -float(fragment[index_center][1]),
                                -float(fragment[index_center][2])])

    n = fragment.shape[0]
    fragment = (fragment.T + np.tile(np.matrix(translation_vec).T, (1, n))).T

    return np.array(fragment)


def perform_rotations(fragment, labels):

    # calculate rotation angle based on atom
    coord_vector = find_coord_vector(ax='z', atom=fragment[labels[0]])
    angle = find_angles(ax='z', point_vector=coord_vector)
    angle = angle * find_rotation_direction(ax='z', atom=fragment[labels[0]])

    # rotate the entire fragment with that angle
    fragment = calculate_rotation(fragment=fragment.copy(), angle=angle, ax='z')

    # Rotation 2
    coord_vector = find_coord_vector(ax='y', atom=fragment[labels[0]])
    angle = find_angles(ax='y', point_vector=coord_vector)
    angle = angle * find_rotation_direction(ax='y', atom=fragment[labels[0]])
    fragment = calculate_rotation(fragment=fragment.copy(), angle=angle, ax='y')

    # Rotation 3
    coord_vector = find_coord_vector(ax='x', atom=fragment[labels[1]])
    angle = find_angles(ax='x', point_vector=coord_vector)
    angle = angle * find_rotation_direction(ax='x', atom=fragment[labels[1]])
    fragment = calculate_rotation(fragment=fragment.copy(), angle=angle, ax='x')

    return fragment


def find_coord_vector(ax, atom):
    coord_vector = np.copy(atom)
    if ax == "x":
        coord_vector[0] = 0
    elif ax == "y":
        coord_vector[1] = 0
    else:
        coord_vector[2] = 0

    return coord_vector


def find_rotation_direction(ax, atom):
    """ Defines the direction of the rotation, clockwise or counter clockwise. """

    if ax == "x" and atom[2] < 0:
        return -1
    elif ax == "y" and atom[2] < 0:
        return -1
    elif ax == "z" and atom[1] < 0:
        return -1

    return 1


def find_angles(point_vector, ax):
    """ Rotates the molecule so that the contact fragment is always in the same position.
        Rotates only the important part of the molecule. """

    assert ax in ["x", "y", "z"], "Ax must be either x, y or z."

    # x, y, z unitary row vectors:
    x, y = np.array([1, 0, 0]), np.array([0, 1, 0])

    norm = np.linalg.norm(point_vector)

    if norm < 1e-10:
        return 0.0
    elif ax == "x":
        # return angle with y axis (beta)
        return np.arccos(np.dot(point_vector, y) / norm)

    # if y or z: return angle with x axis (alpha)
    return np.arccos(np.dot(point_vector, x) / norm)


def calculate_rotation(fragment, angle, ax):
    # rotate all coordinates according to the previously defined rotation angle
    if ax == "x":
        rotation_matrix = rotate_x(angle)
    elif ax == "y":
        rotation_matrix = rotate_y(angle)
    else:
        rotation_matrix = rotate_z(angle)

    # TODO: check why this doesn't work in a single dot product
    for i in range(len(fragment)):
        fragment[i] = np.dot(fragment[i], rotation_matrix)

    # fragment = np.dot(rotation_matrix, fragment.T).T

    return fragment


def rotate_x(angle):
    """ Rotation matrix for rotation around x-axis. """

    return np.array(([1,             0,              0],
                     [0,             np.cos(angle),  -np.sin(angle)],
                     [0,             np.sin(angle),  np.cos(angle)]))


def rotate_y(angle):
    """ Rotation matrix for rotation around y-axis. """

    return np.array(([np.cos(angle),  0,              -np.sin(angle)],
                     [0,              1,              0],
                     [np.sin(angle),  0,              np.cos(angle)]))


def rotate_z(angle):
    """ Rotation matrix for rotation around z-axis. """

    return np.array(([np.cos(angle), -np.sin(angle), 0],
                     [np.sin(angle), np.cos(angle),  0],
                     [0,             0,              1]))


if __name__ == "__main__":
    main()
