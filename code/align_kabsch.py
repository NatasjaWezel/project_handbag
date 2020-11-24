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

    # translate and rotate first fragment onto the origin as for nice viewing
    data[:no_atoms] = perform_translation(data[:no_atoms].copy(), alignment['center'])
    data[:no_atoms] = perform_rotations(data[:no_atoms].copy(), [alignment['yaxis'], alignment['xyplane']])

    # TODO: give R a label

    data_xyz_matrix = np.array([np.array(data.x), np.array(data.y), np.array(data.z)]).T

    A = data_xyz_matrix[0:no_atoms_central]

    print("Applying Kabsch Algorithm...")
    for i in tqdm(range(1, fragments)):
        B_central = data_xyz_matrix[i * no_atoms:no_atoms_central + i * no_atoms]
        B_total = data_xyz_matrix[i * no_atoms: (i + 1) * no_atoms]

        # run kabsch, shift frame in data each time
        B_total_2 = kabsch_align(A, B_central, B_total)

        # TODO: invert?

        # calculate and save error
        # TODO

        # put back into overall matrix
        data_xyz_matrix[i * no_atoms: (i + 1) * no_atoms] = B_total_2

    # put back into df
    data_xyz_matrix_T = data_xyz_matrix.T
    x_vec, y_vec, z_vec = data_xyz_matrix_T[0], data_xyz_matrix_T[1], data_xyz_matrix_T[2]
    data.x, data.y, data.z = x_vec, y_vec, z_vec

    # save as csv
    data.to_csv('test.csv', index=False)

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


def kabsch_align(A, B, B2):
    assert len(A) == len(B)

    N = A.shape[0]  # total points

    # center the points
    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)
    AA = A - np.tile(centroid_A, (N, 1))
    BB = B - np.tile(centroid_B, (N, 1))

    # np.dot is matrix multiplication for array
    H = np.dot(np.transpose(BB), AA)

    # decompose into singular values
    U, S, Vt = np.linalg.svd(H)

    rotation_matrix = np.dot(Vt.T, U.T)

    # special reflection case
    if np.linalg.det(rotation_matrix) < 0:
        Vt[2, :] *= -1
        rotation_matrix = np.dot(Vt.T, U.T)

    translation_vector = -np.dot(rotation_matrix, centroid_B.T) + centroid_A.T

    n = B2.shape[0]

    B2 = (np.dot(rotation_matrix, B2.T)) + np.tile(np.matrix(translation_vector).T, (1, n))

    return B2.T


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


def perform_translation(fragment, label):
    move_x = -float(fragment[fragment.label == label]['x'])
    move_y = -float(fragment[fragment.label == label]['y'])
    move_z = -float(fragment[fragment.label == label]['z'])

    fragment.loc[fragment.label == label, 'x'] = 0
    fragment.loc[fragment.label == label, 'y'] = 0
    fragment.loc[fragment.label == label, 'z'] = 0

    fragment.loc[fragment.label != label, 'x'] = fragment.loc[fragment.label != label, 'x'] + move_x
    fragment.loc[fragment.label != label, 'y'] = fragment.loc[fragment.label != label, 'y'] + move_y
    fragment.loc[fragment.label != label, 'z'] = fragment.loc[fragment.label != label, 'z'] + move_z

    return fragment


def perform_rotations(fragment, labels):
    # calculate rotation angle based on atom
    coord_vector = find_coord_vector(ax='z', atom=fragment[fragment['label'] == labels[0]])
    angle = find_angles(ax='z', coord_vector=coord_vector)
    angle = angle * find_rotation_direction(ax='z', atom=fragment[fragment['label'] == labels[0]])

    # rotate the entire fragment with that angle
    fragment = calculate_rotation(fragment=fragment.copy(), angle=angle, ax='z')

    # Rotation 2
    coord_vector = find_coord_vector(ax='y', atom=fragment[fragment['label'] == labels[0]])
    angle = find_angles(ax='y', coord_vector=coord_vector)
    angle = angle * find_rotation_direction(ax='y', atom=fragment[fragment['label'] == labels[0]])
    fragment = calculate_rotation(fragment=fragment.copy(), angle=angle, ax='y')

    # Rotation 3
    coord_vector = find_coord_vector(ax='x', atom=fragment[fragment['label'] == labels[1]])
    angle = find_angles(ax='x', coord_vector=coord_vector)
    angle = angle * find_rotation_direction(ax='x', atom=fragment[fragment['label'] == labels[1]])
    fragment = calculate_rotation(fragment=fragment.copy(), angle=angle, ax='x')

    return fragment


def find_coord_vector(ax, atom):
    if ax == "x":
        return [0, float(atom.y), float(atom.z)]
    elif ax == "y":
        return [float(atom.x), 0, float(atom.z)]

    return [float(atom.x), float(atom.y), 0]


def find_rotation_direction(ax, atom):
    """ Defines the direction of the rotation, clockwise or counter clockwise. """

    if ax == "x" and float(atom.z) < 0:
        return -1
    elif ax == "y" and float(atom.z) < 0:
        return -1
    elif ax == "z" and float(atom.y) < 0:
        return -1
    else:
        return 1


def find_angles(coord_vector, ax):
    """ Rotates the molecule so that the contact fragment is always in the same position.
        Rotates only the important part of the molecule. """

    # x, y, z unitary row vectors:
    x, y = np.array([1, 0, 0]), np.array([0, 1, 0])

    point_vector = np.array(coord_vector)

    norm = np.linalg.norm(point_vector)

    if norm < 1e-10:
        return 0.0
    elif ax == "x":
        # return angle with y axis (beta)
        return np.arccos(np.dot(point_vector, y) / norm)
    elif ax == "y" or ax == "z":
        # return angle with x axis (alpha)
        return np.arccos(np.dot(point_vector, x) / norm)
    else:
        assert ax in ["x", "y", "z"], "Ax must be either x, y or z."
    assert ax in ["x", "y", "z"], "Ax must be either x, y or z."

    # x, y, z unitary row vectors:
    x, y = np.array([1, 0, 0]), np.array([0, 1, 0])

    point_vector = np.array(coord_vector)

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
    for i, atom in fragment.iterrows():
        coord_vector = np.array([atom.x, atom.y, atom.z])

        if ax == "x":
            coord_vector = np.dot(coord_vector, rotate_x(angle))
        elif ax == "y":
            coord_vector = np.dot(coord_vector, rotate_y(angle))
        else:
            coord_vector = np.dot(coord_vector, rotate_z(angle))

        fragment.loc[fragment.index == i, 'x'] = coord_vector[0]
        fragment.loc[fragment.index == i, 'y'] = coord_vector[1]
        fragment.loc[fragment.index == i, 'z'] = coord_vector[2]

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
