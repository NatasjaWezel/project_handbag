import numpy as np


def perform_rotations(fragment, atoms):
    """ Performs three rotations to lie three of the atoms in the xy plane, one of those
        on the x-axis.
        First rotation: puts first atom on xy-plane if it already was on a plane,
        and above the x-axis if it wasn't by rotating around the z-axis.
        Second rotation: puts first atom on x-axis by rotating around y-axis.
        Third rotation: puts second atom in x-y plane by rotating around the x-axis. """

    atoms = [fragment.atoms[atoms[0]], fragment.atoms[atoms[0]], fragment.atoms[atoms[1]]]
    axis = ["z", "y", "x"]

    for i, atom in enumerate(atoms):

        coord_vector = find_coord_vector(ax=axis[i], atom=atom)
        angle = find_angles(ax=axis[i], coord_vector=coord_vector)
        angle = angle * find_rotation_direction(ax=axis[i], atom=atom)
        fragment = calculate_rotation(fragment=fragment, angle=angle, ax=axis[i])

    return fragment


def find_coord_vector(ax, atom):
    """ Returns the coordinate vectore from the origin to the point. If the atom is not in
        any plane, returns the coordinate vector projected onto a plane. """

    coord_vector = [atom.x, atom.y, atom.z]

    # if not in any plane project it onto the xy plane for the first rotation
    if ax == "x":
        coord_vector = [0, atom.y, atom.z]
    elif ax == "y":
        coord_vector = [atom.x, 0, atom.z]
    elif ax == "z":
        coord_vector = [atom.x, atom.y, 0]

    return coord_vector


def find_rotation_direction(ax, atom):
    """ Defines the direction of the rotation, clockwise or counter clockwise. """

    if ax == "x" and atom.z < 0:
        return -1
    elif ax == "y" and atom.z < 0:
        return -1
    elif ax == "z" and atom.y < 0:
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


def calculate_rotation(fragment, angle, ax):
    # rotate all coordinates according to the previously defined rotation matrices
    for atom in fragment.atoms.values():
        coord_vector = np.array([atom.x, atom.y, atom.z])

        if ax == "x":
            coord_vector = np.dot(coord_vector, rotate_x(angle))
        elif ax == "y":
            coord_vector = np.dot(coord_vector, rotate_y(angle))
        else:
            coord_vector = np.dot(coord_vector, rotate_z(angle))

        atom.x, atom.y, atom.z = coord_vector[0], coord_vector[1], coord_vector[2]

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
