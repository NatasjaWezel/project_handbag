import numpy as np
import pandas as pd


def read_raw_data(filename, no_atoms):
    """ Reads the raw datafiles making use of a logical function. One part reads the fragments their coordinates
        and the other part reads the structure names. """

    print("Pandas is reading csv..." + filename)

    # use amount of atoms for skiprow function to read into df immediately
    data = pd.read_csv(filename, sep='\\s+', usecols=[0, 1, 2, 3], names=['_id', 'x', 'y', 'z'],
                       header=None,
                       skiprows=lambda x: logic(x, no_atoms))

    # get the actual symbols of each atom
    data['symbol'] = data.apply(lambda row: get_atom_symbol(row), axis=1)

    structures = pd.read_csv(filename, sep='*', usecols=[0], names=['structure_id'],
                             header=None,
                             skiprows=lambda x: not logic(x, no_atoms))

    structures['structure_id'] = structures['structure_id'].str.rstrip()

    structures['rmse'] = 0
    structures['mirrored'] = False

    print("Done")

    return data, structures


def get_atom_symbol(row):
    """ Regex the row id to get the actual name of the element. """

    if row._id[1] in '0123456789':
        return row._id[:1]

    return row._id[:2]


def logic(index, amount_atoms):
    """ Logic function to skip each header line - or read only the headers for the structure names - when reading
        in the data. """
    if index % (amount_atoms + 1) == 0:
        return True
    return False


def kabsch_align(A, B, B2, n):
    """ Performs the kabsch algorithm on two central groups, then translates and multiplies the
        entire fragment with the calculated translation vector and rotation matrix. """\

    assert len(A) == len(B), "Fragment 1 and fragment to align do not have the same length"

    # center the points
    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)
    AA = A - np.tile(centroid_A, (n, 1))
    BB = B - np.tile(centroid_B, (n, 1))

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


def perform_translation(fragment, index_center):
    """ Lays the atom on index_center on the origin and moves the rest of the atoms as well. """

    translation_vec = np.array([-float(fragment[index_center][0]),
                                -float(fragment[index_center][1]),
                                -float(fragment[index_center][2])])

    n = fragment.shape[0]
    fragment = (fragment.T + np.tile(np.matrix(translation_vec).T, (1, n))).T

    return np.array(fragment)


def perform_rotations(fragment, labels):
    """ Performs three rotations to lie three of the atoms in the xy plane, one of those
        on the x-axis.
        First rotation: puts first atom on xy-plane if it already was on a plane,
        and above the x-axis if it wasn't by rotating around the z-axis.
        Second rotation: puts first atom on x-axis by rotating around y-axis.
        Third rotation: puts second atom in x-y plane by rotating around the x-axis. """

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


def calc_rmse(A, B, n):
    """ Calculate the RMSE of two matrices. """

    err = A - B
    err = np.multiply(err, err)
    err = np.sum(err)

    return np.sqrt(err / n)
